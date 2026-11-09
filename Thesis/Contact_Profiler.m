%% BENCHMARK_CONTACT_DISTRIBUTION
% Benchmark contact matrix products using a measured active-set
% distribution.
%
% Required variables:
%
%   active_counts : number of active constraints at each time sample
%   n_contact     : total number of contact constraints
%
% Optional variable:
%
%   Scomp        : actual dense contact compliance matrix
%
% If Scomp is not available, a random dense matrix is generated.

clearvars -except active_counts t Scomp n_contact;
clc;

%% Configuration

if ~exist('active_counts', 'var')
    error('The variable active_counts is missing.');
end

active_counts = active_counts(:);

if ~exist('n_contact', 'var')
    error(['Set n_contact to the total number of contact constraints, ' ...
           'for example: n_contact = numel(contact.gaps);']);
end

if any(active_counts < 0) || ...
        any(active_counts > n_contact) || ...
        any(active_counts ~= round(active_counts))
    error('active_counts contains invalid values.');
end

% Number of bins used to approximate the measured distribution.
% Increase this value for a more detailed benchmark.
n_bins = 20;

% Block size for the block-wise product
chunk_size = 512;

% Set this to true only if Scomp is not already available.
generate_test_matrix = ~exist('Scomp', 'var');

if generate_test_matrix
    memory_GB = 8 * n_contact^2 / 1e9;

    fprintf(['Generating a random dense Scomp matrix: %d x %d ' ...
             '(approximately %.2f GB)\n'], ...
        n_contact, n_contact, memory_GB);

    Scomp = rand(n_contact, n_contact);
else
    if ~isequal(size(Scomp), [n_contact, n_contact])
        error('Scomp has dimensions [%d x %d], expected [%d x %d].', ...
            size(Scomp, 1), size(Scomp, 2), n_contact, n_contact);
    end

    fprintf('Using the supplied Scomp matrix.\n');
end

%% Build representative active-set cases

% Assign each sample to an equally spaced active-count bin.
bin_id = floor(active_counts / n_contact * n_bins) + 1;
bin_id = min(bin_id, n_bins);

cases = struct( ...
    'n_active', {}, ...
    'n_samples', {}, ...
    'weight', {});

for ib = 1:n_bins

    samples = active_counts(bin_id == ib);

    if isempty(samples)
        continue;
    end

    % Use the median active count as representative of this bin.
    n_active = round(median(samples));

    cases(end+1).n_active = n_active; %#ok<SAGROW>
    cases(end).n_samples = numel(samples);
    cases(end).weight = numel(samples) / numel(active_counts);
end

fprintf('\nMeasured active-set distribution:\n');
fprintf('  Number of samples: %d\n', numel(active_counts));
fprintf('  Minimum active constraints: %d\n', min(active_counts));
fprintf('  Maximum active constraints: %d\n', max(active_counts));
fprintf('  Mean active constraints: %.2f\n', mean(active_counts));

%% Benchmark each representative active-set size

results = [];

pen = randn(n_contact, 1);

for ic = 1:numel(cases)

    n_active = cases(ic).n_active;

    % Generate a representative active set.
    idx_active = sort(randperm(n_contact, n_active)).';

    idx_inactive = true(n_contact, 1);
    idx_inactive(idx_active) = false;
    idx_inactive = find(idx_inactive);

    % Contact forces on active constraints.
    fa = abs(randn(n_active, 1));

    fprintf('\nCase %d/%d: %d active constraints (%.2f%%, weight %.2f%%)\n', ...
        ic, numel(cases), n_active, ...
        100 * n_active / n_contact, ...
        100 * cases(ic).weight);

    % Original implementation: dense full force vector.
    time_dense = timeit(@() product_dense( ...
        Scomp, pen, idx_active, fa));

    % Full sparse force vector, with the all-active case skipped.
    time_sparse = timeit(@() product_sparse( ...
        Scomp, pen, idx_active, idx_inactive, fa));

    % Only inactive rows and active columns.
    time_indexed = timeit(@() product_inactive( ...
        Scomp, pen, idx_active, idx_inactive, fa));

    % Same operation, computed in blocks.
    time_chunked = timeit(@() product_inactive_chunked( ...
        Scomp, pen, idx_active, idx_inactive, fa, chunk_size));

    % Numerical consistency checks.
    y_dense = product_dense( ...
        Scomp, pen, idx_active, fa);

    y_sparse = product_sparse( ...
        Scomp, pen, idx_active, idx_inactive, fa);

    y_indexed = product_inactive( ...
        Scomp, pen, idx_active, idx_inactive, fa);

    y_chunked = product_inactive_chunked( ...
        Scomp, pen, idx_active, idx_inactive, fa, chunk_size);

    % Compare only the inactive constraints because the optimized
    % implementations return only those values.
    if isempty(idx_inactive)
        err_sparse = 0;
        err_indexed = 0;
        err_chunked = 0;
    else
        reference = y_dense(idx_inactive);

        err_sparse = norm(y_sparse(idx_inactive) - reference) / ...
            max(norm(reference), eps);

        err_indexed = norm(y_indexed - reference) / ...
            max(norm(reference), eps);

        err_chunked = norm(y_chunked - reference) / ...
            max(norm(reference), eps);
    end

    results = [results; ...
        n_active, ...
        cases(ic).n_samples, ...
        cases(ic).weight, ...
        1e3 * time_dense, ...
        1e3 * time_sparse, ...
        1e3 * time_indexed, ...
        1e3 * time_chunked, ...
        time_dense / max(time_sparse, eps), ...
        time_dense / max(time_indexed, eps), ...
        time_dense / max(time_chunked, eps), ...
        err_sparse, ...
        err_indexed, ...
        err_chunked]; %#ok<AGROW>
end

%% Display detailed results

results_table = array2table(results, ...
    'VariableNames', { ...
        'n_active', ...
        'n_samples', ...
        'weight', ...
        'dense_ms', ...
        'sparse_ms', ...
        'inactive_indexed_ms', ...
        'inactive_chunked_ms', ...
        'speedup_sparse', ...
        'speedup_indexed', ...
        'speedup_chunked', ...
        'err_sparse', ...
        'err_indexed', ...
        'err_chunked'});

fprintf('\nDetailed benchmark results:\n');
disp(results_table);

%% Compute distribution-weighted average costs

weighted_dense = sum(results_table.weight .* results_table.dense_ms);
weighted_sparse = sum(results_table.weight .* results_table.sparse_ms);
weighted_indexed = sum(results_table.weight .* ...
    results_table.inactive_indexed_ms);
weighted_chunked = sum(results_table.weight .* ...
    results_table.inactive_chunked_ms);

fprintf('\nDistribution-weighted average cost:\n');
fprintf('  Dense full vector       : %.6f ms\n', weighted_dense);
fprintf('  Sparse full vector      : %.6f ms\n', weighted_sparse);
fprintf('  Inactive indexed        : %.6f ms\n', weighted_indexed);
fprintf('  Inactive block-wise     : %.6f ms\n', weighted_chunked);

fprintf('\nDistribution-weighted speedups relative to dense:\n');
fprintf('  Sparse full vector      : %.2fx\n', ...
    weighted_dense / weighted_sparse);

fprintf('  Inactive indexed        : %.2fx\n', ...
    weighted_dense / weighted_indexed);

fprintf('  Inactive block-wise     : %.2fx\n', ...
    weighted_dense / weighted_chunked);

%% Plot the measured distribution

figure('Color', 'w');
histogram(active_counts, ...
    'BinMethod', 'sturges', ...
    'FaceColor', [0.2 0.45 0.85], ...
    'EdgeColor', 'none');

xlabel('Number of active constraints');
ylabel('Number of samples');
title('Measured active-set distribution');
grid on;


%% ========================================================================
function y = product_dense(Scomp, pen, idx_active, fa)
%PRODUCT_DENSE Original implementation using a dense full force vector.

    f_full = zeros(size(pen));
    f_full(idx_active) = fa;

    y = pen - Scomp * f_full;
end


function y = product_sparse( ...
        Scomp, pen, idx_active, idx_inactive, fa)
%PRODUCT_SPARSE Use a sparse full force vector.
%
% The product is skipped when there are no inactive constraints.

    if isempty(idx_inactive)
        y = zeros(0, 1);
        return;
    end

    if isempty(idx_active)
        y = pen;
        return;
    end

    n_contact = numel(pen);

    f_sparse = sparse( ...
        idx_active, ...
        ones(numel(idx_active), 1), ...
        fa, ...
        n_contact, ...
        1);

    y = pen - Scomp * f_sparse;
end


function y = product_inactive( ...
        Scomp, pen, idx_active, idx_inactive, fa)
%PRODUCT_INACTIVE
% Compute corrected penetration only for inactive constraints.

    if isempty(idx_inactive)
        y = zeros(0, 1);
        return;
    end

    if isempty(idx_active)
        y = pen;
        return;
    end

    y = pen(idx_inactive) - ...
        Scomp(idx_inactive, idx_active) * fa;
end


function y = product_inactive_chunked( ...
        Scomp, pen, idx_active, idx_inactive, fa, chunk_size)
%PRODUCT_INACTIVE_CHUNKED
% Compute corrected penetration for inactive constraints in blocks.

    if isempty(idx_inactive)
        y = zeros(0, 1);
        return;
    end

    if isempty(idx_active)
        y = pen;
        return;
    end

    y = pen(idx_inactive);

    for first = 1:chunk_size:numel(idx_active)

        last = min(first + chunk_size - 1, numel(idx_active));
        jj = first:last;

        y = y - ...
            Scomp(idx_inactive, idx_active(jj)) * fa(jj);
    end
end