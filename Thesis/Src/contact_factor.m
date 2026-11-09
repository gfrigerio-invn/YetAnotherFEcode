function [F, cache] = contact_factor(act, k_pen, Scomp, cache)
%CONTACT_FACTOR Factor of the compliant contact matrix on an active set.
%
%   [F, cache] = CONTACT_FACTOR(act, k_pen, Scomp, cache)
%   x = contact_factor('solve', F, b)        % A_aa \ b
%   cache = contact_factor('new')            % empty cache
%   contact_factor('reset')                  % clears the counters
%   s = contact_factor('stats')              % factorisations and cache hits
%
% INPUT
%   act    [n_c x 1] logical active set, at least one true
%   k_pen  [n_c x 1] penalty stiffness per constraint
%   Scomp  [n_c x n_c] contact compliance, exactly symmetric
%   cache  struct from a previous call, or [] to start an empty one
%
% OUTPUT
%   F      what 'solve' needs to apply A_aa^-1, with
%          A_aa = diag(1/k_pen(act)) + Scomp(act,act): its upper Cholesky
%          factor F.R, or for a small set the matrix itself, F.A
%   cache  updated cache, to be passed back on the next call
%
% WHY. With the compliance correction every contact evaluation solves with
% A_aa, and A_aa depends on the active set ONLY - not on q, not on t. Within a
% Newton iteration, and between neighbouring steps, the set is almost always
% the same, yet the matrix used to be rebuilt and refactorised from scratch at
% every sweep of contact_solve and once more in contact_tangent. At the fine
% mesh interface (n_c = 11240) a factorisation costs 0.31 s with half of the
% constraints active and 1.7 s with all of them, against 5.5 ms and 24 ms for
% the two triangular solves that a stored factor needs.
%
% EXACTNESS. The cache changes the cost and nothing else: the results are
% bit-identical to the A_aa \ b this replaces. That holds because chol followed
% by the two linsolve calls of 'solve' reproduces mldivide's own Cholesky path
% bit for bit - measured from n = 17 up to 3000, with 1 to 1064 right-hand
% sides. Below that it does NOT: up to n = 16 mldivide takes a different route
% and the last bits differ in practically every case, which on a full run grew
% to a 1.7e-7 difference in the response. Small sets therefore keep the matrix
% and solve with mldivide exactly as before; they cost microseconds, so the
% cache loses nothing there. The threshold is set at 32 to leave a margin.
%
% Two sets are kept, not one. When the first guess of contact_solve (pen > 0)
% is inconsistent, the call factorises that guess and then the corrected set;
% while a node is being pushed in by its neighbours this repeats on every
% call, and keeping only the last set would miss the first one each time.
%
% The cache belongs to ONE contact problem (one Scomp, one k_pen). Callers keep
% it in the workspace of their integration, so it dies with it and is never
% shared between models.


    N_SMALL = 32;

    persistent stats
    if isempty(stats), stats = fresh_stats(); end
    if ischar(act)
        switch lower(act)
            case 'solve'
                % Called as contact_factor('solve', F, b), so the arguments
                % are renamed. The two triangular solves are the same LAPACK
                % calls as mldivide's Cholesky path, hence the same bits.
                Fs = k_pen;  b = Scomp;
                if isempty(Fs.R)
                    F = Fs.A \ b;
                else
                    o_t = struct('UT', true, 'TRANSA', true);
                    o_n = struct('UT', true);
                    F = linsolve(Fs.R, linsolve(Fs.R, b, o_t), o_n);
                end
            case 'new',   F = empty_cache();
            case 'reset', stats = fresh_stats(); F = [];
            case 'stats', F = stats;
            otherwise, error('contact_factor:BadQuery', 'Unknown query ''%s''.', act);
        end
        return
    end

    if nargin < 4 || isempty(cache), cache = empty_cache(); end

    for i = 1:numel(cache.sets)
        if isequal(cache.sets{i}, act)
            F = cache.F{i};
            stats.hits = stats.hits + 1;
            return
        end
    end

    A = diag(1 ./ k_pen(act)) + Scomp(act, act);
    n_a = size(A, 1);
    if n_a < N_SMALL
        F = struct('R', [], 'A', A);
    else
        [R, not_pd] = chol(A);
        if not_pd > 0
            error('contact_factor:NotPD', ...
                ['diag(1/k) + S is not positive definite on the active set ' ...
                 '(%d constraints). S is a Gram matrix by construction, so this ' ...
                 'means the compliance passed in is not the one interface_reduction ' ...
                 'returned.'], n_a);
        end
        F = struct('R', R, 'A', []);
    end
    stats.factorisations = stats.factorisations + 1;
    stats.n_a_sum = stats.n_a_sum + n_a;

    % Newest first; the oldest of the two is dropped.
    cache.sets = [{act}, cache.sets(1:min(end, 1))];
    cache.F    = [{F},   cache.F(1:min(end, 1))];
end

function c = empty_cache()
% Factors of the last two active sets, newest first, and the last tangent
% (contact_tangent). The one place that defines the layout.
    c = struct('sets', {{}}, 'F', {{}}, 'keff_act', [], 'K_eff', []);
end

function s = fresh_stats()
    s = struct('factorisations', 0, 'hits', 0, 'n_a_sum', 0);
end
