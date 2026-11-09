function [N, g, info] = contact_operator(Struct, ifaces)
%CONTACT_OPERATOR Build the linear contact operator of a set of rigid walls.
%
%   [N, g, info] = CONTACT_OPERATOR(Struct, ifaces)
%
% Every contact node contributes one unilateral constraint, written as one row
% of N.
%
%       p = N*u - g                                              (penetration)
%
% is positive where the node has crossed its wall
%
% Reduced models: N*Pc is.
%
% INTERFACE SPECIFICATION
%   ifaces is a struct array, one entry per wall, with fields:
%
%     set      name of the node set of the contacting surface.
%     normal   [1 x n_dim] direction along which penetration is measured.
%     gap      distance from the undeformed node to the wall, measured along
%              normal. Scalar (same for every node) or [n_nodes x 1].
%     plane    ALTERNATIVE to gap, for a wall that is NOT parallel to the
%              surface: struct('point', p0) with p0 a point of the wall plane.
%              The gap is then computed node by node,
%                   g_i = normal * (p0 - x_i)',
%     dofs     'normal' (default) keeps one DOF per node, the one along the
%              normal, and requires the normal to be aligned with a DOF of the
%              node: an axis, or any direction once local contact frames are
%              set (FeStructure.set_contact_frames), which is the way to treat
%              an oblique wall with one DOF per node.
%              'all' keeps every DOF of the node, which is what an oblique
%              normal needs, since the contact force then has components on
%              all of them.
%
% OUTPUT
%   N      [n_c x n_free] sparse contact operator, one row per contact node
%   g      [n_c x 1] positive gaps
%   info


    ifaces = normalize_interface_spec(ifaces);
    n_if   = numel(ifaces);
    nD     = Struct.MeshObj.nDOFPerNode;
    n_free = size(Struct.AssemblyObj.constrain_matrix(Struct.K), 1);

    rows = []; cols = []; vals = [];
    g = [];
    info = struct('nodes', {cell(1,n_if)}, 'bnd_dofs', [], ...
                  'blocks', zeros(1,n_if), 'rows', {cell(1,n_if)}, ...
                  'labels', {cell(1,n_if)}, 'spec', ifaces);
    bnd = [];
    r0 = 0;

    for i = 1:n_if
        s = ifaces(i);
        nrm = s.normal(:)' / norm(s.normal);
        ids = Struct.get_nodes(s.set);
        x   = Struct.nodes(ids, :);

        % --- gaps, uniform or from a wall plane --------------------------
        if ~isempty(s.plane)
            gi = (s.plane.point(:)' - x) * nrm(:);
        else
            gi = s.gap(:);
            if isscalar(gi), gi = gi * ones(numel(ids), 1); end
            if numel(gi) ~= numel(ids)
                error('contact_operator:GapSize', ...
                    ['Interface ''%s'': gap has %d entries but the set has %d ' ...
                     'nodes.'], s.set, numel(gi), numel(ids));
            end
        end
        if any(gi <= 0)
            error('contact_operator:NodeInsideWall', ...
                ['Interface ''%s'': %d nodes have a non-positive gap, i.e. they ' ...
                 'start on or beyond the wall. Check the sign of the normal: it ' ...
                 'must point from the body towards the wall.'], ...
                s.set, nnz(gi <= 0));
        end

        % --- which DOFs of a node take part -------------------------------
        if ~any(strcmpi(s.dofs, {'normal', 'all'}))
            error('contact_operator:BadDofs', ...
                'Interface ''%s'': dofs must be ''normal'' or ''all''.', s.set);
        end
        has_frames = isprop(Struct, 'frame_index') && ~isempty(Struct.frame_index);

        % --- assemble the rows, skipping clamped nodes --------------------
        kept = false(numel(ids), 1);
        n_blk0 = numel(bnd);
        for j = 1:numel(ids)
            n_loc = nrm(:);
            if has_frames && Struct.frame_index(ids(j)) > 0
                n_loc = Struct.frame_R(:, :, Struct.frame_index(ids(j)))' * n_loc;
            end
            if strcmpi(s.dofs, 'normal')
                [~, ax] = max(abs(n_loc));
                if abs(abs(n_loc(ax)) - 1) >= 1e-10
                    error('contact_operator:ObliqueNormal', ...
                        ['Interface ''%s'': dofs = ''normal'' keeps only the DOF ' ...
                         'along the normal, which requires the normal to be ' ...
                         'aligned with a DOF of the node. This one is [%s] in the ' ...
                         'node frame. Set local contact frames ' ...
                         '(cfg.local_contact_frames = true) or use dofs = ''all''.'], ...
                         s.set, num2str(n_loc', ' %.4g'));
                end
                kk = ax;
            else
                kk = (1:nD)';
            end
            nk = n_loc(kk); nk = nk(:);
            cd = Struct.AssemblyObj.free2constrained_index((ids(j)-1)*nD + kk');
            cd = cd(:);

            % A row is created only for the DOFs that are free AND carry a
            % component of the normal. A node whose normal-bearing DOFs are
            % all clamped can never come into contact.
            use = cd > 0 & nk ~= 0;
            if ~any(use), continue; end
            kept(j) = true;
            r0 = r0 + 1;
            rows = [rows; repmat(r0, nnz(use), 1)];  %#ok<AGROW>
            cols = [cols; cd(use)];                  %#ok<AGROW>
            vals = [vals; nk(use)];                  %#ok<AGROW>

            % The interface partition keeps every free DOF selected by dofs,
            % including those with a zero normal component: the contact force
            % is transmitted to them through the stiffness of the substructure.
            bnd = [bnd; cd(cd > 0)];                 %#ok<AGROW>
        end

        n_skip = nnz(~kept);
        if n_skip > 0
            fprintf(['  Interface %-4s: %d of %d nodes are clamped and cannot ' ...
                'come into contact, dropped.\n'], s.set, n_skip, numel(ids));
        end
        g = [g; gi(kept)];                           %#ok<AGROW>

        info.nodes{i}  = ids(kept);
        info.labels{i} = s.set;
        info.blocks(i) = numel(bnd) - n_blk0;
        info.rows{i}   = (r0 - nnz(kept) + 1):r0;
    end

    N = sparse(rows, cols, vals, r0, n_free);

    [u_bnd, ia] = unique(bnd, 'stable');
    if numel(u_bnd) ~= numel(bnd)
        error('contact_operator:SharedDof', ...
            ['%d DOFs belong to more than one interface. Overlapping contact ' ...
             'sets would make the interface partition ambiguous; remove the ' ...
             'shared nodes from one of the sets.'], numel(bnd) - numel(u_bnd));
    end
    info.bnd_dofs = bnd(sort(ia));
end
