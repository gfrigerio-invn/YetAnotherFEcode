classdef FeStructure < handle
    % FESTRUCTURE Import a finite element mesh and build the YAFEC model.
    %
    % Works in 2D and 3D: the number of DOFs per node is deduced from the nodal
    % coordinates, and the element constructor branches accordingly (the 2D
    % elements of YAFEC take the out-of-plane thickness, the 3D ones do not).
    %
    % MESH SOURCES, selected by file extension:
    %   .inp   Abaqus input deck, read through abqmesh. Node sets are parsed
    %          here rather than by abqmesh, because abqmesh only matches
    %          assembly-level sets and the compact '*Nset, ..., generate' form
    %          has to be expanded (reading its three numbers as three node IDs
    %          silently produces wrong sets).
    %   .mat   MATLAB file holding the nodal coordinates and the connectivity.
    %          Accepted variable names: nodes/Nodes/coord and
    %          elements/Elements/elem/conn.
    %
    % NODE SETS come from two places, and can be mixed:
    %   - sets declared inside the .inp file. Names starting with
    %     'ContactInterface' become sets labelled by their suffix
    %     ('ContactInterface_T' -> 'T', bare 'ContactInterface' -> 'C');
    %     'BottomFixed' and 'TopFixed' become boundary conditions.
    %   - geometric selectors declared in set_specs before build(). A selector
    %     is a struct with a 'box' field, [n_dim x 2] of [min max] per
    %     direction, with -inf/inf for the unbounded ones; a struct array is
    %     the union of its boxes. This is mesh-format independent and survives
    %     remeshing, unlike node IDs.
    %
    % Which sets are clamped is declared in bc_sets.
    %
    % Typical use with an Abaqus deck (sets already in the file):
    %   S = FeStructure();
    %   S.mesh_file = 'model.inp';  S.element_type = 'TRI3';
    %   S.build();
    %
    % Typical use with a bare mesh and geometric selection:
    %   S = FeStructure();
    %   S.mesh_file = 'model.mat';  S.element_type = 'WED15';
    %   S.set_specs.anchor = struct('box', [-inf inf; -inf inf; -inf -474e-6]);
    %   S.set_specs.stopZ  = struct('box', [ x1  x2 ;  y1  y2 ;  ztop-tol inf]);
    %   S.bc_sets = {'anchor'};
    %   S.build();

    properties
        % --- Input ---
        mesh_file    = '';
        element_type = '';          % TRI3 TRI6 QUAD4 QUAD8 | WED15 HEX8 HEX20 TET4 TET10

        % --- Material ---
        thickness = 10e-6;          % out-of-plane thickness [m], 2D only
        E   = 165e9;                % Young's modulus [Pa]
        rho = 2330;                 % density [kg/m^3]
        nu  = 0.25;                 % Poisson's ratio

        % --- Node set definitions, applied at build() ---
        set_specs = struct();       % label -> geometric selector
        bc_sets   = {};             % labels to clamp in every direction

        % --- Mesh ---
        nodes
        elements
        n_dim                       % 2 or 3, from the nodal coordinates
        bc_nodes

        % --- Node sets ---
        node_sets  = struct();      % label -> node IDs
        set_labels = {};            % sorted list of the labels available

        % --- YAFEC objects ---
        MeshObj
        AssemblyObj
        K
        M
        C

        % --- Analysis results ---
        frequencies
        mode_shapes

        % --- Local contact frames (set_contact_frames) ---
        % Trot maps the local nodal DOFs to the global ones on the FULL
        % (unconstrained) numbering, u_global = Trot * u_local; empty when no
        % frame has been set, in which case every DOF is global. K, M and C are
        % stored in the local DOFs once it is set.
        Trot        = []
        frame_nodes = []            % node IDs that carry a local frame
        frame_R     = []            % [nD x nD x n] rotation of each, columns = local axes
        frame_index = []            % [n_nodes x 1] index into frame_R, 0 = no frame
        k_diag_ref  = []            % max(diag(Kc)) BEFORE the rotation, for the penalty scale
        Trot_c      = []            % Trot restricted to the free DOFs (cached)
    end

    methods
        function obj = FeStructure(varargin)
            if nargin > 0 && mod(nargin,2) == 0
                for i = 1:2:nargin
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end

        function build(obj)
            % BUILD Import the mesh, resolve the node sets, assemble the model.
            if isempty(obj.mesh_file)
                error('FeStructure:NoFile', 'Set mesh_file before calling build().');
            end
            obj.import_mesh();
            obj.apply_set_specs();
            obj.resolve_bc_nodes();
            obj.setup_yafec_model();
        end

        % =================================================================
        %  MESH IMPORT
        % =================================================================
        function import_mesh(obj)
            [~, ~, ext] = fileparts(obj.mesh_file);
            switch lower(ext)
                case '.inp', obj.import_inp();
                case '.mat', obj.import_mat();
                otherwise
                    error('FeStructure:BadExtension', ...
                        'Unsupported mesh file extension ''%s''. Use .inp or .mat.', ext);
            end

            obj.n_dim = size(obj.nodes, 2);
            if ~ismember(obj.n_dim, [2 3])
                error('FeStructure:BadDimension', ...
                    'Nodal coordinates have %d columns; expected 2 or 3.', obj.n_dim);
            end

            fprintf('--- Mesh imported from %s ---\n', obj.mesh_file);
            fprintf(' %dD | %d nodes | %d elements (%s)\n', ...
                obj.n_dim, size(obj.nodes,1), size(obj.elements,1), upper(obj.element_type));
        end

        function import_mat(obj)
            % IMPORT_MAT Nodal coordinates and connectivity from a .mat file.
            S = load(obj.mesh_file);
            obj.nodes    = pick_field(S, {'nodes','Nodes','coord','Coord','XYZ'}, ...
                                      'nodal coordinates');
            obj.elements = pick_field(S, {'elements','Elements','elem','Elem','conn'}, ...
                                      'element connectivity');

            if ~isnumeric(obj.nodes) || ~ismatrix(obj.nodes)
                error('FeStructure:BadNodes', 'The nodal coordinates must be a numeric matrix.');
            end
            obj.nodes    = double(obj.nodes);
            obj.elements = double(obj.elements);   % le mesh TDK arrivano int64
            bad = obj.elements < 1 | obj.elements > size(obj.nodes,1);
            if any(bad(:))
                error('FeStructure:BadConnectivity', ...
                    'The connectivity references %d node IDs outside 1..%d.', ...
                    nnz(bad), size(obj.nodes,1));
            end

            % Sets may travel with the mesh; anything else is defined by selectors.
            obj.node_sets = struct(); obj.set_labels = {}; obj.bc_nodes = [];
            if isfield(S, 'node_sets') && isstruct(S.node_sets)
                f = fieldnames(S.node_sets);
                for i = 1:numel(f)
                    obj.add_node_set(f{i}, S.node_sets.(f{i}));
                end
                fprintf(' %d node set(s) found inside the .mat\n', numel(f));
            end
        end

        function import_inp(obj)
            % IMPORT_INP Abaqus deck: mesh through abqmesh, node sets parsed here.
            meshinfo = abqmesh(obj.mesh_file);
            obj.nodes    = meshinfo.nodes;
            obj.elements = meshinfo.elem{1};

            obj.node_sets = struct(); obj.set_labels = {}; obj.bc_nodes = [];

            fid = fopen(obj.mesh_file, 'r');
            if fid == -1
                error('FeStructure:CannotOpen', ...
                    'Cannot open %s to read the node sets.', obj.mesh_file);
            end
            cleaner = onCleanup(@() fclose(fid));

            label = ''; is_bc = false; is_gen = false;
            while ~feof(fid)
                line = strtrim(fgetl(fid));
                if isempty(line) || startsWith(line, '**'), continue; end

                if startsWith(line, '*Nset', 'IgnoreCase', true)
                    [label, is_bc, is_gen] = parse_nset_header(line);
                elseif startsWith(line, '*')
                    label = ''; is_bc = false; is_gen = false;
                elseif is_bc || ~isempty(label)
                    nums = sscanf(line, '%f,');
                    if is_gen, nums = expand_generate(nums); end
                    if is_bc
                        obj.bc_nodes = [obj.bc_nodes; nums];
                    else
                        obj.add_node_set(label, nums);
                    end
                end
            end
            if ~isempty(obj.bc_nodes), obj.bc_nodes = unique(obj.bc_nodes); end
        end

        % =================================================================
        %  NODE SETS
        % =================================================================
        function add_node_set(obj, label, ids)
            % ADD_NODE_SET Append node IDs to a set, creating it if new.
            label = matlab.lang.makeValidName(label);
            if ~isfield(obj.node_sets, label)
                obj.node_sets.(label) = [];
                obj.set_labels{end+1} = label;
                obj.set_labels = sort(obj.set_labels);
            end
            obj.node_sets.(label) = unique([obj.node_sets.(label); ids(:)]);
        end

        function ids = select_nodes(obj, spec)
            % SELECT_NODES Node IDs matching a geometric selector.
            %   spec(k).box  [n_dim x 2] of [min max] per direction, -inf/inf
            %                for the unbounded ones. A struct array is the
            %                union of its boxes.
            if isempty(obj.nodes)
                error('FeStructure:NoMesh', 'Import the mesh before selecting nodes.');
            end
            keep = false(size(obj.nodes,1), 1);
            for k = 1:numel(spec)
                if ~isfield(spec(k), 'box')
                    error('FeStructure:BadSpec', 'A selector needs a ''box'' field.');
                end
                b = spec(k).box;
                if ~isequal(size(b), [obj.n_dim 2])
                    error('FeStructure:BadBox', ...
                        'box must be %dx2 for a %dD mesh, got %s.', ...
                        obj.n_dim, obj.n_dim, mat2str(size(b)));
                end
                in = true(size(obj.nodes,1), 1);
                for d = 1:obj.n_dim
                    in = in & obj.nodes(:,d) >= b(d,1) & obj.nodes(:,d) <= b(d,2);
                end
                keep = keep | in;
            end
            ids = find(keep);
        end

        function apply_set_specs(obj)
            % APPLY_SET_SPECS Turn the declared geometric selectors into node sets.
            if isempty(fieldnames(obj.set_specs)), return; end
            f = fieldnames(obj.set_specs);
            for i = 1:numel(f)
                ids = obj.select_nodes(obj.set_specs.(f{i}));
                if isempty(ids)
                    warning('FeStructure:EmptySelector', ...
                        'Selector ''%s'' matched no node. Check the box against the mesh extent.', f{i});
                end
                obj.add_node_set(f{i}, ids);
            end
        end

        function resolve_bc_nodes(obj)
            % RESOLVE_BC_NODES Decide which nodes are clamped.
            %
            % Two mutually exclusive sources, so that the outcome never depends
            % on what the mesh file happens to call its sets:
            %   bc_sets empty     the boundary conditions come from the file
            %                     (an Abaqus set whose name contains BC), which
            %                     import_inp has already put into bc_nodes;
            %   bc_sets non-empty those sets are authoritative and replace the
            %                     ones from the file. A .mat mesh carries no
            %                     boundary conditions at all, so this is the
            %                     only route available there.
            if isempty(obj.bc_sets), return; end

            n_file = numel(obj.bc_nodes);
            ids = [];
            for i = 1:numel(obj.bc_sets)
                ids = [ids; obj.get_nodes(obj.bc_sets{i})]; %#ok<AGROW>
            end
            obj.bc_nodes = unique(ids);

            if n_file > 0
                fprintf(['Boundary conditions taken from bc_sets (%d nodes); ' ...
                    'the %d nodes of the file BC set are ignored.\n'], ...
                    numel(obj.bc_nodes), n_file);
            end
        end

        function ids = get_nodes(obj, label)
            % GET_NODES Node IDs of a set.
            label = obj.check_label(label);
            ids = obj.node_sets.(label);
        end

        function dofs = get_set_dofs(obj, label, dirs)
            % GET_SET_DOFS Constrained DOFs of a set, in the given directions.
            %   dirs  vector of directions, 1 = X, 2 = Y, 3 = Z. Defaults to all.
            %
            % The DOFs come out ordered node by node, and within a node in the
            % order given by dirs. Returns [] if every DOF is clamped.
            if nargin < 3 || isempty(dirs), dirs = 1:obj.n_dim; end
            if any(dirs < 1 | dirs > obj.n_dim)
                error('FeStructure:BadDir', ...
                    'Directions must lie in 1..%d for a %dD mesh.', obj.n_dim, obj.n_dim);
            end

            ids = obj.get_nodes(label);
            if isempty(ids), dofs = []; return; end

            nD = obj.MeshObj.nDOFPerNode;
            g = (ids(:)-1)*nD + dirs(:)';      % [n_nodes x numel(dirs)]
            g = reshape(g.', [], 1);           % node-major ordering
            dofs = obj.AssemblyObj.free2constrained_index(g);
            dofs = dofs(:);
            dofs = dofs(dofs > 0);
        end

        function describe_node_sets(obj)
            % DESCRIBE_NODE_SETS Geometric summary of every set.
            % Meant to catch a misread or mis-selected set immediately: the
            % nodes of a face must be aligned, i.e. have one coordinate nearly
            % constant.
            if isempty(obj.set_labels)
                fprintf('No node set defined.\n'); return;
            end
            ax = 'XYZ';
            fprintf('\n--- Node sets ---\n');
            fprintf('%-10s %6s  ', 'Label', 'Nodes');
            for d = 1:obj.n_dim, fprintf('%-25s', [ax(d) ' range [m]']); end
            fprintf('\n');
            for i = 1:numel(obj.set_labels)
                l = obj.set_labels{i};
                n = obj.node_sets.(l);
                if isempty(n), continue; end
                fprintf('%-10s %6d  ', l, numel(n));
                for d = 1:obj.n_dim
                    fprintf('[%9.3e %9.3e]   ', min(obj.nodes(n,d)), max(obj.nodes(n,d)));
                end
                fprintf('\n');
            end
            if ~isempty(obj.bc_nodes)
                fprintf('%-10s %6d  (clamped in all directions)\n', 'BC', numel(obj.bc_nodes));
            end
            fprintf('\n');
        end

        % =================================================================
        %  MODEL
        % =================================================================
        function setup_yafec_model(obj)
            % SETUP_YAFEC_MODEL Mesh, boundary conditions, assembly, matrices.
            myMaterial = KirchoffMaterial();
            set(myMaterial, 'YOUNGS_MODULUS', obj.E, 'DENSITY', obj.rho, ...
                'POISSONS_RATIO', obj.nu);

            % The 2D elements take the out-of-plane thickness, the 3D ones do
            % not and use their own default Gauss rule when called with the
            % material alone.
            switch upper(obj.element_type)
                case 'TRI3',  myMaterial.PLANE_STRESS = true; ctor = @()Tri3Element(obj.thickness, myMaterial);
                case 'TRI6',  myMaterial.PLANE_STRESS = true; ctor = @()Tri6Element(obj.thickness, myMaterial);
                case 'QUAD4', myMaterial.PLANE_STRESS = true; ctor = @()Quad4Element(obj.thickness, myMaterial);
                case 'QUAD8', myMaterial.PLANE_STRESS = true; ctor = @()Quad8Element(obj.thickness, myMaterial);
                case 'WED15', ctor = @()Wed15Element(myMaterial);
                case 'HEX8',  ctor = @()Hex8Element(myMaterial);
                case 'HEX20', ctor = @()Hex20Element(myMaterial);
                case 'TET4',  ctor = @()Tet4Element(myMaterial);
                case 'TET10', ctor = @()Tet10Element(myMaterial);
                otherwise
                    error('FeStructure:BadElement', ...
                        'Unsupported element type ''%s''.', obj.element_type);
            end

            obj.MeshObj = Mesh(obj.nodes);
            obj.MeshObj.create_elements_table(obj.elements, ctor);

            if ~isempty(obj.bc_nodes)
                obj.MeshObj.set_essential_boundary_condition( ...
                    obj.bc_nodes, 1:obj.MeshObj.nDOFPerNode, 0);
            else
                warning('FeStructure:NoBC', ...
                    ['No node is clamped: the model has rigid body modes. ' ...
                     'Several of the reduction methods assume a constrained ' ...
                     'substructure and will fail on a singular stiffness.']);
            end

            obj.AssemblyObj = Assembly(obj.MeshObj);
            obj.M = obj.AssemblyObj.mass_matrix();
            u0 = zeros(obj.MeshObj.nDOFs, 1);
            [obj.K, ~] = obj.AssemblyObj.tangent_stiffness_and_force(u0);

            obj.fix_element_orientation();

            obj.AssemblyObj.DATA.K = obj.K;
            obj.AssemblyObj.DATA.M = obj.M;

            fprintf(' %d DOFs total, %d free\n', obj.MeshObj.nDOFs, ...
                size(obj.AssemblyObj.constrain_matrix(obj.K), 1));
        end

        function fix_element_orientation(obj)
            % FIX_ELEMENT_ORIENTATION Detect and correct a globally inverted mesh.
            %
            % The diagonal of a consistent mass matrix is strictly positive:
            % each entry integrates rho*N_i^2 over the element. If the
            % connectivity orders the nodes so that the Jacobian determinant
            % is negative, the assembly integrates with det(J) < 0 and BOTH K
            % and M come out as the exact negative of the correct matrices.
            % Only the volume factor carries the wrong sign; the shape
            % function gradients, which depend on inv(J), stay right. Negating
            % the assembled matrices is therefore an exact correction, not a
            % patch.
            %
            % This failure is SILENT in a modal analysis, because
            % K*v = lambda*M*v is invariant when both sides change sign: the
            % frequencies come out correct and nothing looks wrong. Everything
            % else is broken. In particular max(diag(K)), which the runs use to
            % scale the contact penalty, comes out negative, and a negative
            % penalty pulls the node THROUGH the wall instead of pushing it
            % back. On the TDK accelerometer mesh this was the case on all
            % 1955 elements.
            d = diag(obj.M);
            d = d(d ~= 0);              % unreferenced nodes contribute nothing
            if isempty(d) || all(d > 0), return; end

            if all(d < 0)
                obj.K = -obj.K;
                obj.M = -obj.M;
                fprintf([' [orientation] Negative Jacobian on every element: K and M ' ...
                    'have been negated.\n               Frequencies were unaffected, ' ...
                    'everything else was silently wrong.\n']);
                return
            end

            error('FeStructure:MixedOrientation', ...
                ['%d of %d mass diagonal entries are negative: the mesh has ' ...
                 'INCONSISTENT element orientation, and a global sign change ' ...
                 'cannot fix it. The connectivity has to be repaired element ' ...
                 'by element.'], nnz(d < 0), numel(d));
        end

        function set_contact_frames(obj, ifaces)
            % SET_CONTACT_FRAMES Express the DOFs of every contact node in a
            % local frame whose first axis is the contact normal.
            %
            %   S.set_contact_frames(cfg.interfaces)
            %
            % WHY. A unilateral contact without friction loads a node only along
            % the wall normal, so one DOF per node is all the interface of a CMS
            % basis needs - but only if that DOF exists. With a wall normal to a
            % coordinate axis it does; with an oblique wall the normal
            % displacement mixes all the Cartesian DOFs of the node, and the
            % interface would have to keep all of them (dofs = 'all'), tripling
            % n_bnd, the Guyan pencil and R_res for nothing. Rotating the node
            % DOFs into (normal, tangent, tangent) makes the normal a DOF again.
            %
            % WHAT CHANGES. K, M (and C if already set) become T'*(.)*T, with T
            % block diagonal: one rotation per contact node, identity elsewhere.
            % Everything that reads S.K / S.M afterwards - ROM bases, Guyan
            % pencil, eigenvalues, Rayleigh damping, the FOM - works in local
            % DOFs without knowing it. What has to be converted is only what
            % crosses the boundary of the model: the shock forcing on the way in
            % (shock_forcing), the displacements on the way out
            % (extract_contact_response, the saved u_full), and the contact
            % operator (contact_operator reads the frame of each node).
            %
            % FRAME CHOICE. The normal goes in the slot of its largest Cartesian
            % component and the other two slots are completed by Gram-Schmidt on
            % the remaining axes. For a wall normal to an axis the rotation is
            % then a pure sign change of one DOF, which is exact in floating
            % point: switching the frames on changes no result. That is also the
            % non-regression test.
            %
            % Call once, after build() and before building contacts, forcing or
            % ROMs.
            if ~isempty(obj.Trot)
                error('FeStructure:FramesSet', 'Local contact frames are already set.');
            end
            if isempty(obj.K) || isempty(obj.M)
                error('FeStructure:NoMatrices', 'Run build() before set_contact_frames.');
            end
            spec   = normalize_interface_spec(ifaces);
            nD     = obj.MeshObj.nDOFPerNode;
            n_full = size(obj.K, 1);
            n_node = size(obj.nodes, 1);

            idx = zeros(n_node, 1);  Rs = zeros(nD, nD, 0);  ids_all = zeros(0, 1);
            n_clamped = 0;
            for i = 1:numel(spec)
                nrm = zeros(nD, 1);
                nn  = spec(i).normal(:);
                nrm(1:min(nD, numel(nn))) = nn(1:min(nD, numel(nn)));
                nrm = nrm / norm(nrm);
                R   = frame_from_normal(nrm);
                ids = obj.get_nodes(spec(i).set);
                for j = 1:numel(ids)
                    id = ids(j);
                    if idx(id) > 0
                        if norm(Rs(:,:,idx(id)) - R, 'fro') > 1e-12
                            error('FeStructure:SharedFrame', ...
                                ['Node %d belongs to two contact sets with different ' ...
                                 'normals; one node can carry only one frame.'], id);
                        end
                        continue
                    end
                    dofs = (id - 1)*nD + (1:nD);
                    free = obj.AssemblyObj.free2constrained_index(dofs) > 0;
                    if ~any(free)
                        n_clamped = n_clamped + 1;  continue   % dropped by contact_operator
                    end
                    if ~all(free) && norm(R - diag(diag(R)), 'fro') > 0
                        error('FeStructure:PartialClamp', ...
                            ['Node %d is clamped in some directions only and its ' ...
                             'contact normal is oblique: a rotation would mix fixed ' ...
                             'and free DOFs.'], id);
                    end
                    Rs(:, :, end+1) = R;  %#ok<AGROW>
                    idx(id) = size(Rs, 3);
                    ids_all(end+1, 1) = id; %#ok<AGROW>
                end
            end

            % Block-diagonal T from triplets: identity on every DOF that is not
            % rotated, one nD x nD block per framed node.
            framed_dof = false(n_full, 1);
            I = []; J = []; V = [];
            for j = 1:numel(ids_all)
                dofs = (ids_all(j) - 1)*nD + (1:nD);
                framed_dof(dofs) = true;
                [r, c] = ndgrid(dofs, dofs);
                R = Rs(:, :, idx(ids_all(j)));
                I = [I; r(:)]; J = [J; c(:)]; V = [V; R(:)]; %#ok<AGROW>
            end
            plain = find(~framed_dof);
            T = sparse([I; plain], [J; plain], [V; ones(numel(plain), 1)], n_full, n_full);
            T = T .* (abs(T) > 0);                     % drop explicit zeros of the blocks

            obj.k_diag_ref = full(max(diag(obj.AssemblyObj.constrain_matrix(obj.K))));
            % No symmetrisation afterwards: the assembled matrices are symmetric
            % only to roundoff, and (A+A')/2 would change their last bits even
            % where T is a pure sign change - which is exactly the case that
            % must reproduce the run without frames bit for bit.
            obj.K = T' * obj.K * T;
            obj.M = T' * obj.M * T;
            obj.AssemblyObj.DATA.K = obj.K;
            obj.AssemblyObj.DATA.M = obj.M;
            if ~isempty(obj.C)
                obj.C = T' * obj.C * T;
                obj.AssemblyObj.DATA.C = obj.C;
            end
            obj.Trot        = T;
            obj.Trot_c      = obj.AssemblyObj.constrain_matrix(T);
            obj.frame_nodes = ids_all;
            obj.frame_R     = Rs;
            obj.frame_index = idx;

            n_obl = 0;
            for j = 1:size(Rs, 3), n_obl = n_obl + (norm(Rs(:,:,j) - diag(diag(Rs(:,:,j))), 'fro') > 0); end
            fprintf(['--- Local contact frames: %d nodes (%d oblique, %d aligned with ' ...
                     'an axis), %d clamped nodes skipped ---\n'], numel(ids_all), n_obl, ...
                     numel(ids_all) - n_obl, n_clamped);
        end

        function u = to_global(obj, u_full)
            % TO_GLOBAL Local nodal DOFs -> global Cartesian, full numbering.
            if isempty(obj.Trot), u = u_full; else, u = obj.Trot * u_full; end
        end

        function u = to_global_constrained(obj, u_c)
            % TO_GLOBAL_CONSTRAINED Same, on the free (constrained) numbering.
            if isempty(obj.Trot), u = u_c; else, u = obj.Trot_c * u_c; end
        end

        function v = to_local(obj, v_full)
            % TO_LOCAL Global Cartesian -> local nodal DOFs, full numbering
            % (T is orthogonal, so its inverse is its transpose).
            if isempty(obj.Trot), v = v_full; else, v = obj.Trot' * v_full; end
        end

        function compute_eigenmodes(obj, n_modes)
            % COMPUTE_EIGENMODES Normal modes of the constrained model.
            if isempty(obj.K) || isempty(obj.M)
                error('FeStructure:NoMatrices', ...
                    'Model matrices are missing. Run build() first.');
            end
            Kc = obj.AssemblyObj.constrain_matrix(obj.K);
            Mc = obj.AssemblyObj.constrain_matrix(obj.M);
            [V0, om] = eigs(Kc, Mc, n_modes, 'SM');
            [obj.frequencies, ind] = sort(sqrt(diag(om)) / (2*pi));
            V0 = V0(:, ind);
            % Always stored in GLOBAL components, so plots and anything that
            % reads node positions see physical directions; the eigenvalues are
            % unaffected by the local frames (orthogonal change of DOFs).
            obj.mode_shapes = obj.to_global(obj.AssemblyObj.unconstrain_vector(V0));

            % Normalized on the largest NODAL displacement, as the reference
            % 3D script does, rather than on the largest single component: the
            % two differ by up to sqrt(n_dim) and only the first is a geometric
            % quantity. Done on the full vector, whose layout is node-major by
            % construction; the constrained one need not be, since a node may
            % have only some of its DOFs clamped.
            nD = obj.MeshObj.nDOFPerNode;
            for ii = 1:n_modes
                v = reshape(obj.mode_shapes(:, ii), nD, []);
                obj.mode_shapes(:, ii) = obj.mode_shapes(:, ii) / ...
                    max(sqrt(sum(v.^2, 1)));
            end

            fprintf('--- Modal analysis complete (%d modes) ---\n', n_modes);
            n_show = min(n_modes, 10);
            for ii = 1:n_show
                fprintf(' Mode %d: %.3f Hz\n', ii, obj.frequencies(ii));
            end
            if n_modes > n_show
                fprintf(' ... Mode %d: %.3f Hz\n', n_modes, obj.frequencies(n_modes));
            end
        end

        function [C, alpha_ray, beta_ray] = compute_rayleigh_damping(obj, Q1, Q2, f_anchor)
            % COMPUTE_RAYLEIGH_DAMPING Rayleigh damping from two quality factors.
            %
            %   compute_rayleigh_damping(Q1, Q2)
            %   compute_rayleigh_damping(Q1, Q2, [fa fb])
            %
            % Q1 is the quality factor AT fa and Q2 the one AT fb, both in Hz.
            % Omitting f_anchor anchors on the first two natural frequencies,
            % which is what every earlier run did and is kept as the default so
            % those results stay reproducible.
            %
            % WHY THE ANCHORS MATTER. C = alpha*M + beta*K gives
            %
            %       1/Q(f) = alpha/(2*pi*f) + 2*pi*f*beta
            %
            % so the two anchors fix alpha and beta and EVERYTHING ELSE IS
            % EXTRAPOLATION. Anchoring on modes 1 and 2 of the 3D model puts
            % both of them at 7.08 and 7.83 kHz, ten per cent apart: asking for
            % Q = 1000 there silently produces Q = 15 at 1 MHz and Q = 3 at
            % 5 MHz, which is 17% of critical on exactly the frequencies the
            % contact excites. Nobody chose that, it fell out of the fit.
            % Passing f_anchor lets the two coefficients be pinned where the
            % physics is, e.g. [30e3 5e6] with Q = [300 5000].
            %
            % The anchors are FREQUENCIES and never mode indices. Rayleigh gives
            %
            %       Phi' C Phi = alpha*I + beta*diag(w_n^2)
            %
            % hence zeta_n = alpha/(2*w_n) + beta*w_n/2, in which the mode index
            % has cancelled: the damping depends on the frequency alone. A mode
            % index would only be a way of NAMING a frequency through the
            % spectrum, and a poor one up high, where the number of modes below
            % f grows like f^3 and locating a 5 MHz anchor would mean computing
            % tens of thousands of them to extract a single value of zeta.
            %
            % The map (zeta1, zeta2) -> (alpha, beta) is linear and invertible,
            % so any pair is reachable from any two distinct anchors. Choosing
            % the anchors changes nothing about what CAN be represented, only
            % about which numbers you have to write to get there.
            %
            % Also returns alpha and beta, which the massless ROMs use to build
            % an equivalent diagonal modal damping.
            if isempty(obj.frequencies) || length(obj.frequencies) < 2
                obj.compute_eigenmodes(2);
            end
            if nargin < 4 || isempty(f_anchor)
                f_anchor = obj.frequencies(1:2);
                anchor_src = 'modes 1-2';
            else
                if numel(f_anchor) ~= 2
                    error('FeStructure:BadAnchor', ...
                        'f_anchor must hold exactly two frequencies [Hz], got %d.', ...
                        numel(f_anchor));
                end
                anchor_src = 'user';
            end
            if any(~isfinite(f_anchor)) || any(f_anchor <= 0)
                error('FeStructure:BadAnchorValue', ...
                    ['Anchor frequencies must be positive and finite, got %s Hz.\n' ...
                     'A zero anchor usually means the model is UNCONSTRAINED and ' ...
                     'modes 1-2 are rigid body modes: apply the boundary ' ...
                     'conditions before asking for the damping, or pass ' ...
                     'f_anchor explicitly.'], mat2str(f_anchor(:)', 4));
            end
            if abs(f_anchor(2) - f_anchor(1)) < eps(max(abs(f_anchor)))
                error('FeStructure:DegenerateAnchor', ...
                    ['The two anchors coincide (%.6g Hz): alpha and beta are ' ...
                     'not separable.'], f_anchor(1));
            end

            w1 = f_anchor(1)*2*pi;   w2 = f_anchor(2)*2*pi;
            z1 = 1/(2*Q1);           z2 = 1/(2*Q2);

            alpha_ray = (2*w1*w2*(z1*w2 - z2*w1)) / (w2^2 - w1^2);
            beta_ray  = (2*(z2*w2 - z1*w1)) / (w2^2 - w1^2);

            obj.C = alpha_ray*obj.M + beta_ray*obj.K;
            obj.AssemblyObj.DATA.C = obj.C;
            C = obj.C;

            fprintf('\n--- Rayleigh damping ---\n');
            fprintf('Anchors (%s): f = %.4g Hz (Q = %g), f = %.4g Hz (Q = %g)\n', ...
                anchor_src, f_anchor(1), Q1, f_anchor(2), Q2);
            fprintf('alpha = %e | beta = %e\n', alpha_ray, beta_ray);
            if alpha_ray < 0 || beta_ray < 0
                warning('FeStructure:NegativeRayleigh', ...
                    ['alpha or beta came out NEGATIVE: the requested pair is ' ...
                     'not realisable as a passive Rayleigh damping.']);
            end

            % The whole point of the anchors is that the rest of the band is
            % extrapolated, so the extrapolation gets printed instead of being
            % discovered later in a transient that will not settle. The span
            % has to reach well past the anchors: pinned on modes 1-2 at 7 kHz
            % the interesting damage happens at 1-5 MHz, three decades higher.
            f_top  = max([50*f_anchor(2), 1e3*f_anchor(1), obj.frequencies(end)]);
            f_show = unique([f_anchor(:)', ...
                10.^(ceil(log10(f_anchor(1))) : floor(log10(f_top)))]);
            fprintf('Resulting quality factor over the band:\n');
            for ff = f_show
                ww = 2*pi*ff;
                if any(abs(ff - f_anchor) < 1e-9*ff), mk = '  <- anchor'; else, mk = ''; end
                fprintf('   %9.4g Hz  ->  Q = %9.4g   (zeta = %.3e)%s\n', ...
                    ff, 1/(alpha_ray/ww + beta_ray*ww), ...
                    0.5*(alpha_ray/ww + beta_ray*ww), mk);
            end
        end

        % =================================================================
        %  PLOTS
        % =================================================================
        function plot_node_sets(obj)
    % PLOT_NODE_SETS Mesh con i nodi appartenenti a ciascun set.

    if isempty(obj.set_labels)
        error('FeStructure:NoSets', ...
              'No node set to plot.');
    end

    figure('Name', 'Node sets', 'Color', 'w');

    PlotMesh(obj.nodes, obj.plot_connectivity(), 0);
    hold on;

    mk = {'o','s','d','^','v','>','<','p'};

    % Handle e label validi solamente per i set effettivamente disegnati
    h = gobjects(0,1);
    leg_labels = cell(0,1);

    for i = 1:numel(obj.set_labels)

        raw_label = obj.set_labels{i};
        n = obj.node_sets.(raw_label);

        % Ignora eventuali set vuoti
        if isempty(n)
            continue;
        end

        c = obj.nodes(n, :);
        marker = mk{mod(i-1, numel(mk)) + 1};

        % Crea il grafico e conserva l'handle restituito
        if obj.n_dim == 3
            hi = plot3(c(:,1), c(:,2), c(:,3), marker, ...
                'MarkerSize', 7, ...
                'LineWidth', 1.3, ...
                'LineStyle', 'none');
        else
            hi = plot(c(:,1), c(:,2), marker, ...
                'MarkerSize', 7, ...
                'LineWidth', 1.3, ...
                'LineStyle', 'none');
        end

        % Aggiunge solo handle validi
        h(end+1,1) = hi;

        % Formattazione dell'etichetta della legenda
        %
        % Esempi:
        % stop_A_top -> stop A top
        % stop_B_bot -> stop B bot
        % stop_Atop  -> stop A top
        %
        token = regexp(raw_label, ...
            '^stop_([A-Da-d])_?(top|bot)$', ...
            'tokens', 'once');

        if ~isempty(token)
            legend_label = sprintf('stop %s %s', ...
                upper(token{1}), lower(token{2}));
        elseif strcmpi(raw_label, 'anchor')
            legend_label = 'anchor';
        else
            % Comportamento generico per eventuali altri nomi
            legend_label = strrep(raw_label, '_', ' ');
        end

        leg_labels{end+1,1} = sprintf('%s (%d)', ...
            legend_label, numel(n));
    end

    % Crea la legenda solo se esistono elementi validi
    if ~isempty(h)
        lgd = legend(h, leg_labels, ...
            'Location', 'bestoutside', ...
            'Interpreter', 'none');

        % Evita che nuovi plot vengano aggiunti automaticamente
        lgd.AutoUpdate = 'off';
    end

    title('Contact Interfaces and Fixed nodes');
    axis equal;
    grid on;

    if obj.n_dim == 3
        view(3);
    end
end
function plot_first_four_modes(obj, scale_factor)
    % PLOT_FIRST_FOUR_MODES Plot dei primi quattro modi in una figura.

    n_modes = 4;

    if isempty(obj.mode_shapes) || size(obj.mode_shapes, 2) < n_modes
        obj.compute_eigenmodes(n_modes);
    end

    if nargin < 2 || isempty(scale_factor)
        scale_factor = 0.2 * ...
            max(max(obj.nodes, [], 1) - min(obj.nodes, [], 1));
    end

    conn = obj.plot_connectivity();
    nD = obj.MeshObj.nDOFPerNode;
    cmap = jet(256);

    figure('Name', 'Primi quattro modi', 'Color', 'w');

    for mode_idx = 1:n_modes
        ax = subplot(2, 2, mode_idx);

        v = reshape(obj.mode_shapes(:, mode_idx), nD, []).';

        PlotFieldonDeformedMesh(obj.nodes, conn, v, ...
            'factor', scale_factor);

        % Impostazioni applicate dopo PlotFieldonDeformedMesh
        set(ax, 'Colormap', cmap);
        clim(ax, [0 1]);

        cb = colorbar(ax);
        cb.Limits = [0 1];

        title(ax, sprintf('Modo %d - %.4g Hz', ...
            mode_idx, obj.frequencies(mode_idx)), ...
            'FontWeight', 'bold');

        axis(ax, 'equal');
        grid(ax, 'on');

        set(findobj(ax, '-property', 'Marker'), ...
            'Marker', 'none');

        if obj.n_dim == 3
            view(ax, 3);
        end
    end

    sgtitle('Primi quattro modi della struttura');
end
        function plot_mode(obj, mode_idx, scale_factor)
            % PLOT_MODE Mode shape mode_idx over the undeformed mesh.
            if isempty(obj.mode_shapes)
                error('FeStructure:NoModes', 'Run compute_eigenmodes() first.');
            end
            % compute_eigenmodes normalizes each mode to unit largest NODAL
            % displacement, so the shape carries no length scale: a factor of 1
            % would draw a one-metre deformation on a micrometre-sized model.
            % Default to a fraction of the model size, which is independent of
            % the units and works in 2D and 3D alike.
            if nargin < 3 || isempty(scale_factor)
                scale_factor = 0.2 * max(max(obj.nodes, [], 1) - min(obj.nodes, [], 1));
            end
            v = reshape(obj.mode_shapes(:, mode_idx), obj.MeshObj.nDOFPerNode, []).';
            figure('Name', sprintf('Mode %d', mode_idx), 'Color', 'w');
            PlotMesh(obj.nodes, obj.plot_connectivity(), 0); hold on;
            PlotFieldonDeformedMesh(obj.nodes, obj.plot_connectivity(), v, ...
                'factor', scale_factor);
            colormap jet; colorbar;
            set(findobj(gca,'-property','Marker'), 'Marker', 'none');
            title(sprintf('\\Phi_{%d} - %.4g Hz', mode_idx, obj.frequencies(mode_idx)));
            axis equal; grid on;
            if obj.n_dim == 3, view(3); end
        end

        function conn = plot_connectivity(obj)
            % PLOT_CONNECTIVITY Node ordering used to draw each element type.
            % In 3D the full connectivity is handed to YAFEC's PlotMesh, which
            % extracts the skin itself. Public because the external plotting
            % helpers need the same ordering to draw on top of the mesh.
            switch upper(obj.element_type)
                case 'TRI3',  conn = obj.elements(:, 1:3);
                case 'TRI6',  conn = obj.elements(:, [1 4 2 5 3 6]);
                case 'QUAD4', conn = obj.elements(:, 1:4);
                case 'QUAD8', conn = obj.elements(:, [1 5 2 6 3 7 4 8]);
                otherwise,    conn = obj.elements;
            end
        end
    end

    % =====================================================================
    methods (Access = private)
        function label = check_label(obj, label)
            if ~ischar(label) && ~isstring(label)
                error('FeStructure:BadLabel', 'The set label must be a string.');
            end
            label = char(label);
            if ~isfield(obj.node_sets, label)
                error('FeStructure:NoSuchSet', ...
                    'Node set ''%s'' is not defined. Available: %s', ...
                    label, strjoin(obj.set_labels, ', '));
            end
        end
    end
end

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================
function v = pick_field(S, names, what)
% Return the first field of S present in names.
for i = 1:numel(names)
    if isfield(S, names{i}), v = S.(names{i}); return; end
end
error('FeStructure:MissingField', ...
    'The .mat file has no %s. Looked for: %s. Found: %s.', ...
    what, strjoin(names, ', '), strjoin(fieldnames(S)', ', '));
end

function [label, is_bc, is_gen] = parse_nset_header(line)
% Extract from '*Nset, nset=NAME, ..., generate' the set name, whether it is a
% boundary set, and whether it uses the compact 'generate' form.
label  = '';
is_bc  = false;
is_gen = ~isempty(regexpi(line, ',\s*generate\s*(,|$)', 'once'));

tok = regexpi(line, 'nset\s*=\s*([^,]+)', 'tokens', 'once');
if isempty(tok), return; end
name = strtrim(tok{1});

if strcmpi(name, 'BottomFixed') || strcmpi(name, 'TopFixed')
    is_bc = true; return;
end

prefix = 'ContactInterface';
if strncmpi(name, prefix, numel(prefix))
    suffix = regexprep(name(numel(prefix)+1:end), '^[_\-]', '');
    if isempty(suffix), label = 'C'; else, label = matlab.lang.makeValidName(suffix); end
else
    label = matlab.lang.makeValidName(name);
end
end

function nums = expand_generate(raw)
% Expand the (first, last, increment) triples of the compact 'generate' form.
raw = raw(:).';
if isempty(raw) || mod(numel(raw), 3) ~= 0
    warning('FeStructure:BadGenerate', ...
        '''generate'' line with %d values (expected a multiple of 3): read as-is.', numel(raw));
    nums = raw(:); return;
end
nums = [];
for i = 1:3:numel(raw)
    step = raw(i+2); if step == 0, step = 1; end
    nums = [nums, raw(i):step:raw(i+1)]; %#ok<AGROW>
end
nums = nums(:);
end

% =====================================================================
function R = frame_from_normal(n)
% Orthonormal frame whose column k is the unit normal n, k being the slot of its
% largest component; the other columns are the remaining coordinate axes made
% orthogonal by Gram-Schmidt, in their natural order. For n = +-e_k this returns
% the identity with a sign on column k, so an aligned wall only flips the sign
% of one DOF.
    nD = numel(n);
    [~, k] = max(abs(n));
    R = eye(nD);
    R(:, k) = n;
    for j = [1:k-1, k+1:nD]
        v = R(:, j);
        for i = [k, 1:j-1]
            if i == j, continue; end
            if i ~= k && i > j, continue; end
            v = v - (R(:, i)' * v) * R(:, i);
        end
        R(:, j) = v / norm(v);
    end
end
