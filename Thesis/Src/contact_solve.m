function [fa, act] = contact_solve(pen, k_pen, Scomp)
%CONTACT_SOLVE Contact force and the active set consistent with it.
%
%   [fa, act] = CONTACT_SOLVE(pen, k_pen, Scomp)
%
% INPUT
%   pen    [n_c x 1] penetration N*q - gap, positive where a constraint is
%          violated
%   k_pen  [n_c x 1] penalty stiffness per constraint
%   Scomp  [n_c x n_c] contact compliance, or [] for the plain penalty law
%
% OUTPUT
%   fa     [nnz(act) x 1] force on the ACTIVE constraints only
%   act    [n_c x 1] logical active set, matching fa in length
%
% Shared by every integrator.
    act = pen > 0;
    if isempty(Scomp)
        fa = k_pen(act) .* pen(act);
        return
    end

    for sweep = 1:40
        if any(act)
            fa = (diag(1 ./ k_pen(act)) + Scomp(act, act)) \ pen(act);
        else
            fa = zeros(0, 1);
        end

        f_full = zeros(numel(pen), 1);
        f_full(act) = fa;
        p_true = pen - Scomp * f_full;

        act_new = act;
        act_new(act  & f_full <= 0) = false;   % pulled, so not touching
        act_new(~act & p_true  > 0) = true;    % pushed in by a neighbour
        if isequal(act_new, act), return; end
        act = act_new;
    end

    % Did not settle in the sweep budget (a degenerate active set
    % can cycle). Return a force CONSISTENT with the final act -
    % solving on it and dropping any pulling node - so that fa and
    % act always match in length. Omitting this is what let a stale
    % fa reach the caller and crash the state-space product.
    if any(act)
        fa = (diag(1 ./ k_pen(act)) + Scomp(act, act)) \ pen(act);
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
