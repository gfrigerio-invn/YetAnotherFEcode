function run_fom(Struct, contact, shock, cfg, run_dir)
%RUN_FOM Reference full-order runs, one per quality factor and penalty.
%
%   RUN_FOM(Struct, contact, shock, cfg, run_dir)
%
% Also runs contact_activity on each result and stores it alongside: knowing
% whether the contact is flat or on a corner is what decides whether a
% reduction study is worth starting, so it should never be an afterthought.
%
% See also TRANSIENTSOLVERODE, CONTACT_ACTIVITY, EXTRACT_CONTACT_RESPONSE.

    fprintf('\n=========================================\n');
    fprintf('                 FOM\n');
    fprintf('=========================================\n');

    Mc = Struct.AssemblyObj.constrain_matrix(Struct.M);
    Kc = Struct.AssemblyObj.constrain_matrix(Struct.K);
    labels = fieldnames(contact.Interfaces)';

    [Q_pairs, f_anchor] = damping_spec(cfg);
    for i_Q = 1:size(Q_pairs, 2)
        Q = Q_pairs(1, i_Q);            % also the tag carried by the file name
        Struct.compute_rayleigh_damping(Q_pairs(1, i_Q), Q_pairs(2, i_Q), f_anchor);
        Cc = Struct.AssemblyObj.constrain_matrix(Struct.C);

        for k_mult = cfg.array_k_mult
            k_contact = contact.k_base * k_mult;
            fprintf('\n[FOM] Q = %d | k_mult = %g\n', Q, k_mult);

            tic;
            % Integrator to run the FOM. Default: ode15s
            integ = 'ode15s';
            if isfield(cfg, 'integrator') && ~isempty(cfg.integrator)
                integ = lower(cfg.integrator);
            end
            switch integ
                case 'ode15s'
                    solver = TransientSolverOde(Mc, Kc, Cc);
                    [t, q] = solver.solve(cfg.tmax, cfg.dt, shock.q0, shock.qd0, shock.handle, ...
                        'ContactOperator', contact.N, ...
                        'ContactGap',      contact.gaps, ...
                        'ContactPenalty',  k_contact, ...
                        'Label',           'FOM', ...
                        'Eref',            shock.Eref, ...
                        'RelTol',          cfg.RelTolFOM, ...
                        'OutputTimes',     shock.t_out);
                case {'newmark', 'genalpha'}
                    h_int = cfg.dt;
                    if isfield(cfg, 'h_int') && ~isempty(cfg.h_int), h_int = cfg.h_int; end
                    rho = 0.7;
                    if isfield(cfg, 'rho_inf') && ~isempty(cfg.rho_inf), rho = cfg.rho_inf; end
                    al = 0;
                    if isfield(cfg, 'newmark_alpha') && ~isempty(cfg.newmark_alpha)
                        al = cfg.newmark_alpha;
                    end
                    solver = TransientSolverNewmark(Mc, Kc, Cc);
                    [t, q] = solver.solve(cfg.tmax, cfg.dt, shock.q0, shock.qd0, shock.handle, ...
                        'ContactOperator', contact.N, ...
                        'ContactGap',      contact.gaps, ...
                        'ContactPenalty',  k_contact, ...
                        'Label',           'FOM', ...
                        'Scheme',          integ, ...
                        'RhoInf',          rho, ...
                        'Alpha',           al, ...
                        'TimeStep',        h_int, ...
                        'OutputTimes',     shock.t_out);
                otherwise
                    error('RUNFOM:BadIntegrator', ...
                        ['cfg.integrator must be ''ode15s'', ''newmark'' or ' ...
                        '''genalpha'', got ''%s''.'], integ);
            end

            cpu_time = toc;

            y_contact = extract_contact_response(Struct, contact.Interfaces, labels, q);

            % SAVE FIRST, diagnose after.
            Interfaces   = contact.Interfaces;
            model        = 'FOM';
            n_modes      = size(Mc, 1);
            offline_time = 0;
            fname = fullfile(run_dir, sprintf('FOM_Q%g_K%g.mat', Q, k_mult));
            save(fname, 't', 'y_contact', 'Interfaces', 'cpu_time', ...
                'offline_time', 'model', 'n_modes', 'Q', 'k_mult');

            % The full displacement field, on the free DOFs and on the output
            % grid, is what animate_contact_3d needs to deform the whole mesh
            % rather than the contact nodes alone. It is a few MB against
            % hours of computation, so it is worth keeping, but it stays
            % optional because the sweep does not need it.
            if isfield(cfg, 'save_full_field') && cfg.save_full_field
                u_full = q; %#ok<NASGU>
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
end
