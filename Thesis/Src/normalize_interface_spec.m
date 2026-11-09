function s = normalize_interface_spec(spec)
%NORMALIZE_INTERFACE_SPEC Bring both accepted formats to one struct array.
% The legacy format is the cell array {label, direction, signed_gap; ...},
% where the sign of the gap told which side of the DOF the wall was on. It
% maps to a normal along that axis with the sign of the gap, and a positive
% gap equal to its magnitude.

    if iscell(spec)
        n = size(spec, 1);
        s = repmat(struct('set', '', 'normal', [], 'gap', [], ...
                          'plane', [], 'dofs', 'normal'), 1, n);
        for i = 1:n
            d = spec{i,2}; gp = spec{i,3};
            if gp == 0
                error('contact_operator:ZeroGap', ...
                    'Interface ''%s'': a zero gap gives no side for the wall.', spec{i,1});
            end
            nrm = zeros(1, max(d, 2));
            nrm(d) = sign(gp);
            s(i).set = spec{i,1};
            s(i).normal = nrm;
            s(i).gap = abs(gp);
        end
        return
    end

    if ~isstruct(spec)
        error('contact_operator:BadSpec', ...
            'The interfaces must be a struct array or the legacy cell array.');
    end

    s = spec(:)';
    for i = 1:numel(s)
        if ~isfield(s, 'dofs') || isempty(s(i).dofs), s(i).dofs = 'normal'; end
        if ~isfield(s, 'plane'), s(i).plane = []; end
        if ~isfield(s, 'gap'),   s(i).gap = []; end
        if isempty(s(i).gap) && isempty(s(i).plane)
            error('contact_operator:NoGap', ...
                'Interface ''%s'': give either gap or plane.', s(i).set);
        end
        if ~isempty(s(i).gap) && ~isempty(s(i).plane)
            error('contact_operator:TwoGaps', ...
                'Interface ''%s'': gap and plane are alternatives, not both.', s(i).set);
        end
        if isempty(s(i).normal)
            error('contact_operator:NoNormal', ...
                'Interface ''%s'': the normal is missing.', s(i).set);
        end
    end
end
