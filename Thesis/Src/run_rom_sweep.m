function run_rom_sweep(Struct, contact, shock, cfg, run_dir)
%RUN_ROM_SWEEP Reduced-order runs over every requested combination.
%
%   RUN_ROM_SWEEP(Struct, contact, shock, cfg, run_dir)
%
% Sweeps quality factor, penalty multiplier, number of linear modes, method
% and number of CC modes, and writes one .mat per combination.
%
% See also BUILD_ROM, INTERFACE_REDUCTION, TRANSIENTSOLVERODE,
% TRANSIENTSOLVERMASSLESS.

    rom_list = {'MT', 'MC', 'CB', 'Rubin', 'MCB', 'MN'};
    rom_list = rom_list(cellfun(@(m) cfg.run.(m), rom_list) == 1);
    if isempty(rom_list), return; end

    % Methods whose basis can take a secondary interface reduction. 
    ir_eligible = {'CB', 'Rubin'};

    fprintf('\n=========================================\n');
    fprintf('                 ROM\n');
    fprintf('=========================================\n');

    max_phi = max(cfg.array_linModes);
    labels  = fieldnames(contact.Interfaces)';

    Kbb_guyan = []; Mbb_guyan = [];
    if cfg.interface_reduction.enabled
        fprintf('Computing the Guyan interface pencil (%d constraint modes)...\n', ...
            contact.n_bnd);
        tic_g = tic;
        [Kbb_guyan, Mbb_guyan] = guyan_interface_pencil(Struct, contact.dofs);
        fprintf('  done in %.1f s\n', toc(tic_g));
    end

    [Q_pairs, f_anchor] = damping_spec(cfg);
    for i_Q = 1:size(Q_pairs, 2)
        Q = Q_pairs(1, i_Q);            % also the tag carried by the file name

        % compute_rayleigh_damping does two things: it updates Struct.C, used
        % by the penalty ROMs, and it returns alpha and beta, from which the
        % massless ROMs build an equivalent diagonal modal damping. Using the
        % same pair keeps the damping identical across all methods.
        [~, alpha_ray, beta_ray] = Struct.compute_rayleigh_damping( ...
            Q_pairs(1, i_Q), Q_pairs(2, i_Q), f_anchor);
        rayleigh = struct('alpha', alpha_ray, 'beta', beta_ray);

        w_hi = 2*pi * Struct.frequencies(max_phi);
        z_hi = 0.5*(alpha_ray/w_hi + beta_ray*w_hi);
        fprintf('  zeta(mode %d, %.4g Hz) = %.4f | requested zeta at the anchors = %.4f, %.4f\n', ...
            max_phi, Struct.frequencies(max_phi), z_hi, ...
            1/(2*Q_pairs(1, i_Q)), 1/(2*Q_pairs(2, i_Q)));
        if z_hi > 1
            warning('ROMSWEEP:Overdamped', ...
                'Rayleigh makes the high modes OVERDAMPED (zeta_%d = %.2f).', ...
                max_phi, z_hi);
        end

        for k_mult = cfg.array_k_mult
            k_contact = contact.k_base * k_mult;

            for phi = cfg.array_linModes
                for im = 1:numel(rom_list)
                    model = rom_list{im};

                    % The massless models use exact set-valued contact. They run at the first value only.
                    is_massless = any(strcmp(model, {'MCB', 'MN'}));
                    if is_massless && k_mult ~= cfg.array_k_mult(1), continue; end

                    fprintf('\n--- %s | Phi %d | Q %d | k_mult %g ---\n', ...
                        model, phi, Q, k_mult);

                    tic_offline = tic;
                    [rom, Pc, Mr, Kr, Cr] = build_rom(Struct, model, phi, ...
                        contact, k_contact, rayleigh);
                    offline_base = toc(tic_offline);

                    % Spatial part of the forcing projected once.
                    F_r_full = Pc' * shock.F_spatial;

                    % n_cc = 0 is the baseline without interface reduction,
                    % always run so the comparison term is present.
                    if cfg.interface_reduction.enabled && any(strcmp(model, ir_eligible))
                        cc_list = [0, cfg.array_ccModes];
                        if isfield(cfg.interface_reduction, 'skip_baseline') && ...
                                cfg.interface_reduction.skip_baseline
                            cc_list = cfg.array_ccModes;
                        end
                    else
                        cc_list = 0;
                    end

                    for n_cc = cc_list
                        ir = apply_interface_reduction(Mr, Kr, Cr, F_r_full, ...
                            rom, model, n_cc, cfg, contact, Pc, ...
                            Kbb_guyan, Mbb_guyan, offline_base);

                        fprintf('  [%s] reduced size %d\n', ir.tag, size(ir.K, 1));

                        z0       = zeros(size(ir.K, 1), 1);
                        F_handle = @(t) ir.F * shock.profile(t);

                        tic;
                        cf = [];
                        if is_massless
                            [t, q_rom, lambda, info] = integrate_massless( ...
                                ir, rom, contact, Pc, z0, F_handle, cfg);
                        else
                            [t, q_rom, cf] = integrate_penalty(ir, rom, contact, ...
                                Pc, n_cc, z0, F_handle, shock, cfg, k_contact);
                            lambda = []; info = [];
                        end
                        cpu_time = toc;

                        % Back through T_CC first and only then through Pc:
                        % applying blkdiag(Phi_CC, I) to the time history is
                        % far cheaper than forming Pc*T_CC.
                        if n_cc == 0
                            q_full = q_rom;
                        else
                            q_full = ir.info.T * q_rom;
                        end
                        u_rom = Pc * q_full;
                        y_contact = extract_contact_response(Struct, ...
                            contact.Interfaces, labels, u_rom);
                        % The force from the solver, not k*(normal - gap): with
                        % the static correction the two differ (see
                        % attach_contact_force).
                        [y_contact, contact_force] = attach_contact_force( ...
                            y_contact, contact, cf, t);

                        % The whole field only on request: at 22302 DOFs and
                        % 1001 samples it is 180 MB per ROM, against a few MB
                        % for y_contact, and only the modal analysis in
                        % test_05_spectralAnalysis needs it. Same layout as the
                        % FOM u_full (free DOFs x output instants), so the two
                        % are projected on the modes by the same code.
                        if ~(isfield(cfg, 'save_rom_full_field') && cfg.save_rom_full_field)
                            u_rom = [];
                        elseif ismethod(Struct, 'to_global_constrained')
                            u_rom = Struct.to_global_constrained(u_rom);   % global components
                        end
                        save_rom_result(run_dir, contact, model, ir, phi, ...
                            n_cc, Q, k_mult, t, y_contact, cpu_time, ...
                            lambda, info, u_rom, contact_force);
                    end
                end
            end
        end
    end
end

% =====================================================================
function ir = apply_interface_reduction(Mr, Kr, Cr, F_r_full, rom, model, ...
                                        n_cc, cfg, contact, Pc, ...
                                        Kbb_guyan, Mbb_guyan, t_off)
%APPLY_INTERFACE_REDUCTION Secondary modal reduction of the interface block.
% Returns the untouched matrices when n_cc is zero, which is the baseline.

    ir = struct('M', Mr, 'K', Kr, 'C', Cr, 'F', F_r_full, 'Phi_CC', [], ...
                'info', [], 'mode', 'none', 'tag', model, 'offline', t_off);
    ir.tag = tag_integrator(ir.tag, cfg);
    if n_cc == 0, return; end

    tic_ir = tic;

    % The Guyan pencil is known in PHYSICAL interface coordinates, while
    % interface_reduction works in the ROM's own interface coordinates. With
    % u_phys = A*u_rom (A being the contact rows of the boundary block of Pc,
    % the identity for CB and the scaling for Rubin), the congruence A'*(.)*A
    % moves the pencil across and leaves the eigenvalues unchanged, so Phi_CC
    % comes back directly in ROM coordinates.
    A = Pc(contact.dofs, 1:rom.n_bnd);
    pencil = struct('K', A' * Kbb_guyan * A, 'M', A' * Mbb_guyan * A);

    do_static = isfield(cfg.interface_reduction, 'static_correction') && ...
                cfg.interface_reduction.static_correction;

    % Refinement of the CC modes (Ahn et al. 2022). Independent of the contact
    % compliance above: that one changes the contact law, this one changes the
    % basis. They compose, and either can be used alone.
    do_refine = isfield(cfg.interface_reduction, 'refine') && ...
                cfg.interface_reduction.refine;

    % Optional: give every contact face the same number of CC modes, instead of
    % letting the frequency pooling decide the split. n_cc is distributed as
    % evenly as possible over the faces (exactly equal when it is a multiple of
    % their number), so "k modes per interface" is array_ccModes = k*n_faces.
    cc_alloc = [];
    equal_pf = isfield(cfg.interface_reduction, 'equal_per_face') && ...
               cfg.interface_reduction.equal_per_face;
    if equal_pf && strcmpi(cfg.interface_reduction.mode, 'per_interface')
        nf = numel(contact.blocks);
        base = floor(n_cc / nf);
        rem  = n_cc - base*nf;                 % spread the remainder one-each
        cc_alloc = base * ones(1, nf);
        cc_alloc(1:rem) = cc_alloc(1:rem) + 1; % faces are ordered as in contact.blocks
    end

    [ir.M, ir.K, ir.C, ir.Phi_CC, ir.info] = interface_reduction( ...
        Mr, Kr, Cr, rom.n_bnd, n_cc, cfg.interface_reduction.mode, ...
        contact.blocks, pencil, do_static, cc_alloc, do_refine);

    ir.mode = ir.info.mode;
    switch ir.mode
        case 'global',        ir.tag = [model 'IRG'];
        case 'per_interface', ir.tag = [model 'IRP'];
        otherwise,            ir.tag = [model 'IR'];
    end
    % A correction changes the model, so it must change the file name too:
    % without this a corrected and an uncorrected run write the SAME file and
    % the second silently overwrites the first. The suffixes are letters only,
    % which is what the post-processing method parser accepts, so each variant
    % appears as its own method next to the plain one. They can both be on:
    %   SC  series contact compliance (corrects the contact law)
    %   RF  refined CC modes (corrects the basis)
    if ~isempty(ir.info)
        if isfield(ir.info, 'R_res') && ~isempty(ir.info.R_res)
            ir.tag = [ir.tag 'SC'];
        end
        if isfield(ir.info, 'refined') && ir.info.refined
            ir.tag = [ir.tag 'RF'];
        end
    end
    ir.tag = tag_integrator(ir.tag, cfg);

    % T' applied to the projected forcing. Taken from the reduction, not rebuilt
    % from Phi_CC: the refined transformation carries an interface/substructure
    % coupling block that blkdiag(Phi_CC, I) does not have.
    ir.F = ir.info.T' * F_r_full;
    ir.offline = t_off + toc(tic_ir);
end

% =====================================================================
function [t, q, cf] = integrate_penalty(ir, rom, contact, Pc, n_cc, z0, ...
                                        F_handle, shock, cfg, k_contact)
%INTEGRATE_PENALTY Penalty contact with ode15s, one path for every method.

    N_run = contact.N * Pc;

    % Compliance of the interface content removed by the truncation, mapped onto
    % the contact constraints: one row and column per constraint. Only the CMS
    % models have an interface partition (MT and MC carry no n_bnd), and only
    % they can be interface-reduced, so everything below is gated on n_cc.
    Scomp = [];
    if n_cc > 0
        % Operator mapping the ROM's INTERFACE coordinates onto the constraints,
        % captured before the CC reduction is applied to it. It is the right map
        % for the compliance: R_res comes out of the interface pencil, so it
        % lives in the ROM's interface coordinates, which for a scaled basis
        % such as Rubin's are NOT the physical ones. Using the physical contact
        % operator here silently produced a compliance twelve orders of
        % magnitude too small.
        Nb_rom = N_run(:, 1:rom.n_bnd);

        % Right-multiplication by the interface transformation. Taken from the
        % reduction rather than assumed block diagonal: after refinement it is
        % not, so N_run(:,1:n_bnd)*Phi_CC would drop the coupling block.
        N_run = N_run * ir.info.T;

        if isfield(ir.info, 'R_res') && ~isempty(ir.info.R_res)
            Scomp = compliance_on_constraints(Nb_rom, ir.info.R_res);
        end
    end
    contact_solve('reset');

    % Same integrator as the FOM, from the same cfg field. Deliberately not a
    % separate setting: if the reference and the reduced model were integrated
    % differently, the measured "reduction error" would also contain the
    % difference between the two schemes, and there would be no way to tell
    % the two apart afterwards.
    integ = 'ode15s';
    if isfield(cfg, 'integrator') && ~isempty(cfg.integrator)
        integ = lower(cfg.integrator);
    end
    switch integ
        case 'ode15s'
            solver = TransientSolverOde(ir.M, ir.K, ir.C);
            [t, q, cf] = solver.solve(cfg.tmax, cfg.dt, z0, z0, F_handle, ...
                'ContactCompliance', Scomp, ...
                'ContactOperator', N_run, ...
                'ContactGap',      contact.gaps, ...
                'ContactPenalty',  k_contact, ...
                'Label',           ir.tag, ...
                'Eref',            shock.Eref, ...
                'RelTol',          cfg.RelTol, ...
                'OutputTimes',     shock.t_out);
        case {'newmark', 'genalpha'}
            rho = 0.7;
            if isfield(cfg,'rho_inf') && ~isempty(cfg.rho_inf), rho = cfg.rho_inf; end
            al = 0;
            if isfield(cfg,'newmark_alpha') && ~isempty(cfg.newmark_alpha)
                al = cfg.newmark_alpha;
            end
            solver = TransientSolverNewmark(ir.M, ir.K, ir.C);
            [t, q, cf] = solver.solve(cfg.tmax, cfg.dt, z0, z0, F_handle, ...
                'ContactCompliance', Scomp, ...
                'ContactOperator', N_run, ...
                'ContactGap',      contact.gaps, ...
                'ContactPenalty',  k_contact, ...
                'Label',           ir.tag, ...
                'Scheme',          integ, ...
                'RhoInf',          rho, ...
                'Alpha',           al, ...
                'RelTol',          newton_tol(cfg), ...
                'ForceStride',     force_stride(cfg, cfg.dt), ...
                'OutputTimes',     shock.t_out);
        otherwise
            error('ROMSWEEP:BadIntegrator', ...
                'cfg.integrator must be ''ode15s'', ''newmark'' or ''genalpha''.');
    end

    if ~isempty(Scomp)
        s = contact_solve('stats');
        if s.calls_contact > 0
            fprintf(['  [SC] active-set sweeps over %d contact calls (%d total): ' ...
                     'mean %.3f | max %d | budget hits %d | hist(1..4) = %s\n'], ...
                s.calls_contact, s.calls, s.sweeps_contact / s.calls_contact, ...
                s.max_sweeps, s.budget_hits, mat2str(s.hist(1:4)));
            % The cost of the correction is the factorisations: each one is
            % O(n_a^3), a cache hit is two triangular solves (contact_factor).
            fprintf('  [SC] factorisations %d (mean n_a %.0f) | cache hits %d\n', ...
                s.factorisations, s.n_a_mean, s.factor_hits);
        else
            fprintf('  [SC] compliance active but no contact occurred in this run\n');
        end
    end
end

% =====================================================================
function S = compliance_on_constraints(Nb, R)
%COMPLIANCE_ON_CONSTRAINTS S = Nb*R*Nb', without the products when Nb is diagonal-like.
%
% On every contact model so far each constraint reads ONE interface coordinate:
% N(:, dofs) is a signed permutation, and N*Pc(:, 1:n_bnd) keeps one entry per
% row - unit for CB, 1/d_i for the scaled Rubin basis. Then
%
%       S(i,j) = s_i * R(p_i, p_j) * s_j
%
% which is a gather instead of two dense n^3 products: 0.4 s against 21 s at the
% fine mesh interface (n_c = 11240), and the same bits when the pattern is
% exact (verified). Rubin's Pc carries inv(T1_bb)*T1_bb, i.e. the identity up to
% roundoff, so its off-pattern entries are rounding noise of order 1e-16; they
% are dropped below a relative 1e-10 and the general product is used whenever
% anything larger is found.
    TOL = 1e-10;
    [n_c, n_b] = size(Nb);
    [top2, idx] = maxk(abs(Nb), 2, 2);           % [n_c x 2], largest two per row
    p = idx(:, 1);
    s = full(Nb(sub2ind([n_c n_b], (1:n_c)', p)));
    off = top2(:, min(2, n_b));
    if n_b > 1 && any(off > TOL * abs(s))
        fprintf('  [SC] contact operator is not one entry per row (max off/diag %.1e): dense product\n', ...
            max(off ./ abs(s)));
        S = full(Nb * R * Nb');
    else
        S = s .* R(p, p) .* s.';
    end
    S = (S + S') / 2;
end


% =====================================================================
function save_rom_result(run_dir, contact, model, ir, phi, n_cc, Q, k_mult, ...
                         t, y_contact, cpu_time, lambda, info, u_full, contact_force)
%SAVE_ROM_RESULT One file per combination, named after it.
% The variable names are the contract with the post-processing.

    Interfaces   = contact.Interfaces;
    model_tag    = ir.tag;
    ir_mode      = ir.mode;
    ir_info      = ir.info;
    % R_res is n_bnd x n_bnd and dense - 1 GB on the fine mesh (11240 interface
    % DOFs) against a few hundred MB for the rest of the file. It is only needed
    % while integrating, and nothing reads it back from disk, so it is not
    % saved; a flag records that the static correction was on.
    if isstruct(ir_info) && isfield(ir_info, 'R_res')
        ir_info.static_correction = ~isempty(ir_info.R_res);
        ir_info = rmfield(ir_info, 'R_res');
    end
    offline_time = ir.offline;
    n_modes      = phi;

    if n_cc == 0
        fname = sprintf('ROM_%s_Phi%03d_Q%g_K%g.mat', ir.tag, phi, Q, k_mult);
    else
        fname = sprintf('ROM_%s_Phi%03d_CC%03d_Q%g_K%g.mat', ...
            ir.tag, phi, n_cc, Q, k_mult);
    end

    save(fullfile(run_dir, fname), ...
        't', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time', ...
        'model', 'model_tag', 'n_modes', 'Q', 'k_mult', 'lambda', 'info', ...
        'n_cc', 'ir_mode', 'ir_info');
    if nargin >= 14 && ~isempty(u_full)
        save(fullfile(run_dir, fname), 'u_full', '-append');
    end
    % Sparse, only the active constraints: a few MB even at the step size.
    if nargin >= 15 && ~isempty(contact_force)
        save(fullfile(run_dir, fname), 'contact_force', '-append');
    end
end

% =====================================================================
function s = force_stride(cfg, h)
%FORCE_STRIDE Integration steps between contact-force samples.
% cfg.contact_force_dt sets the force sampling; empty or absent means the
% output grid. It must be a multiple of h and divide the output spacing, so
% that the force lands exactly on the output instants.
    s = [];
    if ~isfield(cfg, 'contact_force_dt') || isempty(cfg.contact_force_dt), return; end
    s = round(cfg.contact_force_dt / h);
    dt_out = cfg.output_stride * cfg.dt;
    if s < 1 || abs(s*h - cfg.contact_force_dt) > 1e-9*cfg.contact_force_dt || ...
            abs(dt_out/(s*h) - round(dt_out/(s*h))) > 1e-9
        error('ROMSWEEP:BadForceDt', ...
            ['cfg.contact_force_dt = %.4g s must be a multiple of h = %.4g s ' ...
             'and divide the output spacing %.4g s.'], cfg.contact_force_dt, h, dt_out);
    end
end

% =====================================================================
function tag = tag_integrator(tag, cfg)
%TAG_INTEGRATOR Append the integrator to a model tag.

    if isfield(cfg, 'integrator') && ~isempty(cfg.integrator)
        switch lower(cfg.integrator)
            case 'genalpha', tag = [tag 'GA'];
            case 'newmark',  tag = [tag 'NM'];
        end
    end
end
