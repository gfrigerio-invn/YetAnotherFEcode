function run_fom(Struct, contact, shock, cfg, run_dir)
%RUN_FOM Reference full-order runs, one per quality factor and penalty.
%
%   RUN_FOM(Struct, contact, shock, cfg, run_dir)

    fprintf('\n=========================================\n');
    fprintf('                 FOM\n');
    fprintf('=========================================\n');

    Mc = Struct.AssemblyObj.constrain_matrix(Struct.M);
    Kc = Struct.AssemblyObj.constrain_matrix(Struct.K);
    labels = fieldnames(contact.Interfaces)';
    
    % Damping definition
    f_anchor = cfg.Q_freq(:)';
    Q_pair   = cfg.array_QFactor(:);
    
    if numel(f_anchor) ~= 2
        error('RUNFOM:BadQFreq', ...
            'cfg.Q_freq must contain exactly two frequencies.');
    end
    
    if numel(Q_pair) ~= 2
        error('RUNFOM:BadQFactor', ...
            'cfg.array_QFactor must contain one Q value for each anchor frequency.');
    end
    
    Struct.compute_rayleigh_damping( ...
        Q_pair(1), Q_pair(2), f_anchor);
    
        Cc = Struct.AssemblyObj.constrain_matrix(Struct.C);
    
    for k_mult = cfg.array_k_mult
        k_contact = contact.k_base * k_mult;
        fprintf('\n[FOM] Q = %d | k_mult = %g\n', Q, k_mult);
    
        tic;
        integ = lower(cfg.integrator);
            switch integ
                case 'ode15s'
                    solver = TransientSolverOde(Mc, Kc, Cc);
                    [t, q, cf] = solver.solve(cfg.tmax, cfg.dt, shock.q0, shock.qd0, shock.handle, ...
                        'ContactOperator', contact.N, ...
                        'ContactGap',      contact.gaps, ...
                        'ContactPenalty',  k_contact, ...
                        'Label',           'FOM', ...
                        'Eref',            shock.Eref, ...
                        'RelTol',          cfg.RelTolFOM, ...
                        'OutputTimes',     shock.t_out);
                case {'newmark', 'genalpha'}
                    rho = 0.7;
                    if isfield(cfg, 'rho_inf') && ~isempty(cfg.rho_inf), rho = cfg.rho_inf; end
                    al = 0;
                    if isfield(cfg, 'newmark_alpha') && ~isempty(cfg.newmark_alpha)
                        al = cfg.newmark_alpha;
                    end
                    solver = TransientSolverNewmark(Mc, Kc, Cc);
                    [t, q, cf] = solver.solve(cfg.tmax, cfg.dt, shock.q0, shock.qd0, shock.handle, ...
                        'ContactOperator', contact.N, ...
                        'ContactGap',      contact.gaps, ...
                        'ContactPenalty',  k_contact, ...
                        'Label',           'FOM', ...
                        'Scheme',          integ, ...
                        'RhoInf',          rho, ...
                        'Alpha',           al, ...
                        'TimeStep',        h_int, ...
                        'RelTol',          newton_tol(cfg), ...
                        'ForceStride',     fom_force_stride(cfg, h_int), ...
                        'OutputTimes',     shock.t_out);
                otherwise
                    error('RUNFOM:BadIntegrator', ...
                        ['cfg.integrator must be ''ode15s'', ''newmark'' or ' ...
                         '''genalpha'', got ''%s''.'], integ);
            end
            cpu_time = toc;

            y_contact = extract_contact_response(Struct, contact.Interfaces, labels, q);
            % Same fields as the ROMs (.force, .pen), so the post-processing
            % reads one quantity for every model. Without the static correction
            % .pen equals max(normal - gap, 0), which is a check in itself.
            [y_contact, contact_force] = attach_contact_force(y_contact, contact, cf, t);

        % SAVE FIRST, diagnose after.
            Interfaces   = contact.Interfaces;
        model        = 'FOM';
            n_modes      = size(Mc, 1);
            offline_time = 0;
            % The machine cpu_time was measured on. A reference timed on one
            % host and ROMs timed on another give a speedup that is a ratio
            % between two computers; the post-processing checks this.
            host         = run_host();
            fname = fullfile(run_dir, sprintf('%s_Q%g_K%g.mat', model, Q, k_mult));
            save(fname, 't', 'y_contact', 'Interfaces', 'cpu_time', ...
                'offline_time', 'model', 'n_modes', 'Q', 'k_mult', 'host', 'integrator');
            if ~isempty(contact_force)
                save(fname, 'contact_force', '-append');
            end

            % The full displacement field, on the free DOFs and on the output
            % grid, is what animate_contact_3d needs to deform the whole mesh
            % rather than the contact nodes alone. It is a few MB against
            % hours of computation, so it is worth keeping, but it stays
            % optional because the sweep does not need it.
            if isfield(cfg, 'save_full_field') && cfg.save_full_field
                % Saved in GLOBAL components even when the contact nodes carry
                % local frames, so post-processing and animations never need to
                % know about them.
                u_full = q;
                if ismethod(Struct, 'to_global_constrained')
                    u_full = Struct.to_global_constrained(u_full); %#ok<NASGU>
                end
                save(fname, 'u_full', '-append');
                d = dir(fname);
                fprintf('  full displacement field stored (%.1f MB)\n', d.bytes/1e6);
            end
            fprintf('Saved %s\n', fname);

            fprintf('\n--- Contact activity ---\n');
            try
                activity = contact_activity(y_contact, contact.Interfaces, labels, t);
                save(fname, 'activity', '-append');
            catch ME
                warning('RUNFOM:ActivityFailed', ...
                    ['contact_activity failed (%s). The response is saved, ' ...
                     'so it can be rerun on the file.'], ME.message);
            end
        end
    end

% =====================================================================
function s = fom_force_stride(cfg, h)
%FOM_FORCE_STRIDE Steps between contact-force samples, from cfg.contact_force_dt.
% Empty means the output grid. Same rule as in run_rom_sweep.
    s = [];
    if ~isfield(cfg, 'contact_force_dt') || isempty(cfg.contact_force_dt), return; end
    s = round(cfg.contact_force_dt / h);
    dt_out = cfg.output_stride * cfg.dt;
    if s < 1 || abs(s*h - cfg.contact_force_dt) > 1e-9*cfg.contact_force_dt || ...
            abs(dt_out/(s*h) - round(dt_out/(s*h))) > 1e-9
        error('RUNFOM:BadForceDt', ...
            ['cfg.contact_force_dt = %.4g s must be a multiple of h = %.4g s ' ...
             'and divide the output spacing %.4g s.'], cfg.contact_force_dt, h, dt_out);
    end
end
