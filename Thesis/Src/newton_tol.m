function tol = newton_tol(cfg)
%NEWTON_TOL Newton residual tolerance of the second-order schemes.
%
%   tol = NEWTON_TOL(cfg)
%
% Read in one place so that run_fom and run_rom_sweep cannot drift apart: the
% reference and the reduced models must solve each step to the same residual,
% or part of the measured "reduction error" would be a difference in solver
% settings. Before this, the ROMs were handed cfg.RelTol - the ode15s
% tolerance, 1e-4 in some runs - while the FOM used the solver default.
%
% Missing field -> 1e-6, the TransientSolverNewmark default, so every existing
% run_config keeps the FOM it was produced with.
%
% See also TRANSIENTSOLVERNEWMARK, RUN_FOM, RUN_ROM_SWEEP.

    tol = 1e-6;
    if isfield(cfg, 'newton_tol') && ~isempty(cfg.newton_tol)
        tol = cfg.newton_tol;
    end
    if ~(isscalar(tol) && isfinite(tol) && tol > 0)
        error('NEWTONTOL:Bad', 'cfg.newton_tol must be a positive scalar.');
    end
end
