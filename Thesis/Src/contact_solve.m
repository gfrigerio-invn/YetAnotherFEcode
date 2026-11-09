function [fa, act, cache] = contact_solve(pen, k_pen, Scomp, cache)
%CONTACT_SOLVE Contact force and the active set consistent with it.
%
%   [fa, act] = CONTACT_SOLVE(pen, k_pen, Scomp)
%   [fa, act, cache] = CONTACT_SOLVE(pen, k_pen, Scomp, cache)
%
% INPUT
%   pen    [n_c x 1] penetration N*q - gap, positive where a constraint is
%          violated
%   k_pen  [n_c x 1] penalty stiffness per constraint
%   Scomp  [n_c x n_c] contact compliance, or [] for the plain penalty law
%   cache  optional, the factor cache of this contact problem (see
%          contact_factor). Without it every call factorises from scratch,
%          with the same result.
%
% OUTPUT
%   fa     [nnz(act) x 1] force on the ACTIVE constraints only
%   act    [n_c x 1] logical active set, matching fa in length
%   cache  updated cache, to be passed back on the next call
%
% Shared by every integrator so that the force, the tangent and the residual
% can never disagree about which constraints are touching.
    % With the correction the constraints are COUPLED:
    %
    %       p_true = p - S*f
    %
    % and activity has to be decided on p_true, not on p. A node
    % entering the set changes the matrix being inverted, so every
    % activation makes the force on all the other nodes jump.
    persistent stats
    if isempty(stats), stats = fresh_stats(); end
    if nargin == 1 && ischar(pen)
        act = [];  cache = [];
        switch lower(pen)
            case 'reset'
                stats = fresh_stats();  contact_factor('reset');  fa = [];
            case 'stats'
                fa = stats;
                fs = contact_factor('stats');
                fa.factorisations = fs.factorisations;
                fa.factor_hits    = fs.hits;
                fa.n_a_mean       = fs.n_a_sum / max(fs.factorisations, 1);
            otherwise, error('contact_solve:BadQuery', 'Unknown query ''%s''.', pen);
        end
        return
    end
    if nargin < 4, cache = []; end

    % The computation starts HERE

    act = pen > 0;
    if isempty(Scomp) % In case of NO STATIC CORRECTION
        fa = k_pen(act) .* pen(act);
        return
    end
    % In case of STATIC CORRECTION
    BUDGET = 40;
    n_sweep   = 0;
    converged = false;
    for sweep = 1:BUDGET
        n_sweep = sweep;
        if any(act)
            [F, cache] = contact_factor(act, k_pen, Scomp, cache);
            fa = contact_factor('solve', F, pen(act));
        else
            fa = zeros(0, 1);
        end

        f_full = zeros(numel(pen), 1);
        f_full(act) = fa;
        p_true = pen - Scomp * f_full;

        act_new = act;
        act_new(act  & f_full <= 0) = false;   % pulled, so not touching
        act_new(~act & p_true  > 0) = true;    % pushed in by a neighbour
        if isequal(act_new, act), converged = true; break; end % if set is converged
        act = act_new;
    end

    if ~converged
        % Did not settle in the sweep budget (a degenerate active set
        % can cycle). 
        if any(act)
            [F, cache] = contact_factor(act, k_pen, Scomp, cache);
            fa = contact_factor('solve', F, pen(act));
            keep = fa > 0;
            if ~all(keep)
                idx = find(act);
                act(idx(~keep)) = false;
                fa = fa(keep);
            end
        else
            fa = zeros(0, 1);
        end
    end

    in_contact = any(act) || any(pen > 0);
    stats.calls = stats.calls + 1;
    if in_contact
        stats.calls_contact  = stats.calls_contact + 1;
        stats.sweeps_contact = stats.sweeps_contact + n_sweep;
        stats.max_sweeps     = max(stats.max_sweeps, n_sweep);
        stats.hist(min(n_sweep, BUDGET)) = stats.hist(min(n_sweep, BUDGET)) + 1;
        stats.budget_hits    = stats.budget_hits + ~converged;
    end
end

function s = fresh_stats()
    s = struct('calls', 0, 'calls_contact', 0, 'sweeps_contact', 0, ...
               'max_sweeps', 0, 'budget_hits', 0, 'hist', zeros(1, 40));
end
