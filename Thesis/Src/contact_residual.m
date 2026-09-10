function fh = contact_residual(M, K, C, Nop, gap, k_pen, F_ext, Scomp)
%CONTACT_RESIDUAL Residual closure for the YaFEc second-order integrators.
%
%   fh = CONTACT_RESIDUAL(M, K, C, Nop, gap, k_pen, F_ext, Scomp)
%   [r, drdqdd, drdqd, drdq, c0] = fh(q, qd, qdd, t)
%
% Builds the handle ImplicitNewmark and GeneralizedAlpha expect (see
% src/TimeIntegration/ResidualFunctions/residual_nonlinear.m for the
% convention) for a system with unilateral penalty contact:
%
%   r      = M*qdd + C*qd + K*q + F_contact(q) - F_ext(t)
%   drdqdd = M
%   drdqd  = C
%   drdq   = K + d F_contact / dq
%
%
% INPUT
%   M, K, C   [n x n] system matrices, already constrained
%   Nop       [n_c x n] contact operator, penetration is Nop*q - gap
%   gap       [n_c x 1] positive gaps
%   k_pen     scalar or [n_c x 1] penalty stiffness
%   F_ext     function handle of t returning [n x 1]
%   Scomp     [n_c x n_c] contact compliance, or [] for the plain penalty law
%
% See also CONTACT_TANGENT, CONTACT_SOLVE, TRANSIENTSOLVERODE.

    if isscalar(k_pen), k_pen = k_pen * ones(numel(gap), 1); end
    k_pen = k_pen(:);
    gap   = gap(:);
    Nop_t = Nop';                      % applied at every evaluation

    fh = @residual;

    function [r, drdqdd, drdqd, drdq, c0] = residual(q, qd, qdd, t)
        [F_pen, K_eff] = contact_tangent(q, K, Nop, Nop_t, gap, k_pen, Scomp);

        F_inertial = M * qdd;
        F_damping  = C * qd;
        F_elastic  = K * q;
        F_external = F_ext(t);

        r      = F_inertial + F_damping + F_elastic + F_pen - F_external;
        drdqdd = M;
        drdqd  = C;
        drdq   = K_eff;

        c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + ...
             norm(F_pen) + norm(F_external);
        if c0 == 0, c0 = 1; end        % at rest with no load the test is vacuous
    end
end
