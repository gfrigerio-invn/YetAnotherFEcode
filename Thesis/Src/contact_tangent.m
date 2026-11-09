function [F_pen, K_eff, act, fa, cache] = contact_tangent(q, K, Nop, Nop_t, gap, k_pen, Scomp, cache)
%CONTACT_TANGENT Penalty contact force and the tangent consistent with it.
%
%   [F_pen, K_eff, act, fa] = CONTACT_TANGENT(q, K, Nop, Nop_t, gap, k_pen, Scomp)
%   [F_pen, K_eff, act, fa, cache] = CONTACT_TANGENT(..., Scomp, cache)
%
% INPUT
%   q      [n x 1] displacement
%   K      [n x n] structural stiffness
%   Nop    [n_c x n] contact operator, penetration is Nop*q - gap
%   Nop_t  Nop' precomputed, since it is applied at every evaluation
%   gap    [n_c x 1] positive gaps
%   k_pen  [n_c x 1] penalty stiffness per constraint
%   Scomp  [n_c x n_c] contact compliance, or [] for the plain penalty law
%   cache  optional contact cache (see contact_factor). It must belong to THIS
%          K, Nop, k_pen and Scomp: the caller keeps one per integration.
%
% OUTPUT
%   F_pen  [n x 1] contact force on the structure
%   K_eff  [n x n] K plus the contact contribution, i.e. d(K*q + F_pen)/dq
%   act    [n_c x 1] logical active set
%   fa     force on the active constraints
%   cache  updated cache, to be passed back on the next call
%
% The force and the tangent come from ONE call to contact_solve, so they are
% always the exact derivative of one another. 

    if nargin < 8 || isempty(cache), cache = contact_factor('new'); end
    n = size(K, 1);
    [fa, act, cache] = contact_solve(Nop*q - gap, k_pen, Scomp, cache);

    if ~any(act)
        F_pen = zeros(n, 1);
        K_eff = K;
        return
    end

    F_pen = Nop_t(:, act) * fa;
    if isequal(cache.keff_act, act)
        K_eff = cache.K_eff;
        return
    end

    if isempty(Scomp)
        K_eff = K + Nop_t(:, act) * (k_pen(act) .* Nop(act, :));
    else
        [F, cache] = contact_factor(act, k_pen, Scomp, cache);
        Keff_a = contact_factor('solve', F, Nop(act, :));
        K_eff  = K + Nop_t(:, act) * Keff_a;
    end
    cache.keff_act = act;
    cache.K_eff    = K_eff;
end
