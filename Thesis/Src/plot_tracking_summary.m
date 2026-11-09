function fig = plot_tracking_summary(results_dir, varargin)
%PLOT_TRACKING_SUMMARY Tracking-time bars from tracking_time.csv, readable labels.
%
%   fig = PLOT_TRACKING_SUMMARY(results_dir)
%   fig = PLOT_TRACKING_SUMMARY(results_dir, 'Exclude', ["Rubin Phi100", ...])
%
% 'Exclude' drops bars by their CLEANED label, as it appears on the axis.
%
% Same chart as the post-processing's Summary_TrackingTime, redrawn from the
% saved tracking_time.csv with:
%   - the integrator tag dropped (NM / GA are the same for every bar of a run,
%     so they only lengthen the labels), and the 'G' of the global IR mode,
%     which is the only IR mode in these sweeps;
%   - the methods grouped (MT, Rubin, Rubin IR, Rubin IR+SC) and, inside each
%     group, ordered by number of modes - phi first, then n_cc - instead of the
%     alphabetical order of the file names, which put Phi1000 before Phi200.
%
% Writes Summary_TrackingTime_sorted.png/.fig next to the original.
%
% See also ROM_TRACKING_TIME.

    p = inputParser;
    addParameter(p, 'Exclude', strings(0));
    parse(p, varargin{:});
    excl = string(p.Results.Exclude);

    T = readtable(fullfile(results_dir, 'tracking_time.csv'), 'TextType', 'string');

    % "RubinIRGSCNM Phi1000 NM CC16" -> method RubinIRGSC, phi 1000, n_cc 16
    n = height(T);
    meth = strings(n, 1);  phi = zeros(n, 1);  ncc = zeros(n, 1);
    for i = 1:n
        tk = regexp(char(T.rom(i)), '^(\w+?)(?:NM|GA)?\s+Phi(\d+)', 'tokens', 'once');
        if isempty(tk)
            error('plot_tracking_summary:BadName', 'Cannot parse ROM name "%s".', T.rom(i));
        end
        meth(i) = regexprep(tk{1}, 'IRG', 'IR');       % drop the G of global
        phi(i)  = str2double(tk{2});
        % n_cc parsed on its own: MATLAB drops a token captured inside an
        % optional trailing group, so folding it into the pattern above loses it.
        tc = regexp(char(T.rom(i)), 'CC(\d+)', 'tokens', 'once');
        if ~isempty(tc), ncc(i) = str2double(tc{1}); end
    end

    % group order, then modes
    grp_order = ["MT", "MC", "CB", "Rubin", "CBIR", "CBIRSC", "RubinIR", "RubinIRSC", "MCB", "MN"];
    g = zeros(n, 1);
    for i = 1:n
        k = find(grp_order == meth(i), 1);
        if isempty(k), k = numel(grp_order) + 1; end
        g(i) = k;
    end
    [~, ord] = sortrows([g phi ncc]);

    lab = meth(ord) + " Phi" + string(phi(ord));
    has_cc = ncc(ord) > 0;
    lab(has_cc) = lab(has_cc) + " CC" + string(ncc(ord(has_cc)));

    if ~isempty(excl)
        unknown = excl(~ismember(excl, lab));
        if ~isempty(unknown)
            error('plot_tracking_summary:UnknownLabel', ...
                'Nothing to exclude called %s. Available: %s', ...
                strjoin('"' + unknown + '"', ', '), strjoin(lab, ', '));
        end
        keep = ~ismember(lab, excl);
        ord = ord(keep);  lab = lab(keep);
        n = numel(ord);
    end

    Tb  = [T.T_0p1_us(ord), T.T_1_us(ord), T.T_10_us(ord)];
    Cs  = logical([T.cens_0p1(ord), T.cens_1(ord), T.cens_10(ord)]);
    thr = [0.1 1 10];

    fig = figure('Name', 'Tracking time', 'Color', 'w', 'Position', [100 100 900 520]);
    hold on; grid on; box on
    bh = bar(Tb, 'grouped');
    for i_thr = 1:3
        cs = Cs(:, i_thr);
        if any(cs)
            plot(bh(i_thr).XEndPoints(cs), Tb(cs, i_thr), '^k', ...
                'MarkerFaceColor', 'w', 'MarkerSize', 7, 'HandleVisibility', 'off');
        end
        bh(i_thr).DisplayName = sprintf('within %g%%', thr(i_thr));
    end
    % thin separators between method groups
    gs = g(ord);
    for b = find(diff(gs) ~= 0)'
        xline(b + 0.5, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    set(gca, 'XTick', 1:n, 'XTickLabel', lab, 'TickLabelInterpreter', 'none', ...
        'XTickLabelRotation', 30);
    ylabel('time inside the error band [\mus]');
    title(['How long each ROM tracks the FOM    ' ...
        '(open marker = never left the band, so the bar is a lower bound)']);
    legend('Location', 'northwest');
    ax = gca;  ax.Toolbar.Visible = 'off';

    exportgraphics(fig, fullfile(results_dir, 'Summary_TrackingTime_sorted.png'), 'Resolution', 300);
    savefig(fig, fullfile(results_dir, 'Summary_TrackingTime_sorted.fig'));
end
