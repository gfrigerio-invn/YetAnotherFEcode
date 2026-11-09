function [y_contact, contact_force] = attach_contact_force(y_contact, contact, cf, t)
%ATTACH_CONTACT_FORCE Add the solver's contact force to the per-face output.
%
%   [y_contact, contact_force] = ATTACH_CONTACT_FORCE(y_contact, contact, cf, t)
%
% INPUT
%   y_contact  per-face response from extract_contact_response
%   contact    from build_interfaces (contact.info.rows maps faces to rows)
%   cf         third output of TransientSolverNewmark / TransientSolverOde:
%              .t [1 x n_f], .F [n_c x n_f] sparse, .k [n_c x 1]
%   t          output grid of y_contact
%
% OUTPUT
%   y_contact.(face).force  [n_nodes x n_t] normal contact force per node on
%                           the output grid, >= 0 (compressive)
%   y_contact.(face).pen    [n_nodes x n_t] force ./ k, the penetration the
%                           springs actually see
%   contact_force           cf plus .rows and .labels, to be saved as it is:
%                           with a ForceStride finer than the output it is the
%                           only record of the impacts in between
%
% WHY .pen AND NOT normal - gap. With the static correction the springs see
% p - S*f, not the penetration read from the reduced coordinates: the
% displacement in y_contact misses the -S*f term, and at few CC modes it
% overstates the penetration several times over (33x the impulse of Rubin at
% 22 CC on Sim4). Without the correction the two coincide, so .pen is the one
% quantity valid for every model.
%
% Rows of N follow the retained nodes of each face in order (contact_operator),
% which is the node order of y_contact, so row j of a face is node j.
%
% See also TRANSIENTSOLVERNEWMARK, EXTRACT_CONTACT_RESPONSE, CONTACT_OPERATOR.

    if isempty(cf), contact_force = []; return; end

    % Output instants inside the force grid. Exact matches by construction
    % (both are multiples of the step), so no interpolation.
    [ok, idx] = ismembertol(t, cf.t, 1e-6, 'DataScale', max(abs(cf.t(end)), eps));
    if ~all(ok)
        error('attach_contact_force:GridMismatch', ...
            ['%d output instants are not force samples: the force stride ' ...
             'must divide the output spacing.'], nnz(~ok));
    end

    rows = contact.info.rows;
    for i = 1:numel(contact.labels)
        lbl = contact.labels{i};
        if ~isfield(y_contact, lbl), continue; end
        Fi = full(cf.F(rows{i}, idx));
        if size(Fi, 1) ~= size(y_contact.(lbl).normal, 1)
            error('attach_contact_force:RowMismatch', ...
                'Face %s: %d constraint rows but %d nodes in y_contact.', ...
                lbl, size(Fi, 1), size(y_contact.(lbl).normal, 1));
        end
        y_contact.(lbl).force = Fi;
        y_contact.(lbl).pen   = Fi ./ cf.k(rows{i});
    end

    contact_force        = cf;
    contact_force.rows   = rows;
    contact_force.labels = contact.labels;
end
