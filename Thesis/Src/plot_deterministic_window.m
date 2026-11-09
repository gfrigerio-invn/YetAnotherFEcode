function res = plot_deterministic_window(folder, phi)
%PLOT_DETERMINISTIC_WINDOW IR vs IR + static correction, before the chaos sets in.
%
%   res = PLOT_DETERMINISTIC_WINDOW(folder, phi)
%
% folder  results folder, full path or a name under test_05_Accelerometer3D/results.
%         Must hold the FOM, ROM_RubinNM_Phi<phi>, and ROM_RubinIRGNM /
%         ROM_RubinIRGSCNM at n_cc = 16, 32, 64, 128 for that phi.
% phi     number of linear modes of the Rubin ROMs to compare.
%
% Vibro-impact is chaotic, so an error integrated over a long run stops ranking
% the models once the trajectories decorrelate: past that point a truncated ROM
% can "beat" its own parent by chance. The comparison is therefore made on a
% window fixed by the PARENT model - the full Rubin at the same phi, of which
% every IR variant is a subspace - namely the time it stays within 1% of the
% FOM (rom_tracking_time, same definition and HoldFraction as the
% post-processing). While the parent tracks the FOM, any error of an IR variant
% is due to the interface reduction; after it, the chaos amplifies every
% difference and the error measures where two uncorrelated trajectories are.
%
% rom_divergence's T_sat is NOT used: it is where the error reaches HALF its
% saturation level, the midpoint of the divergence rather than its onset. On the
% phi = 1000, 1e6 g run the parent already carries a 14% GRE over [0, T_sat],
% and the ordering of the variants is already scrambled there.
%
% Writes, into the results folder:
%   Deterministic_1_result.png      GRE on the window vs n_cc, both variants
%   Deterministic_2_horizon.png     pointwise error, with the window marked
%   Deterministic_3_window_sweep.png  GRE over [0, T] for every T: no window chosen
%   Deterministic_4_bars.png        GRE on the window as bars, one pair per n_cc
%   Deterministic_window.txt        the numbers, and a check that the GRE_win
%                                   of the post-processing is reproduced
%
% See also ROM_TRACKING_TIME, ROM_DIVERGENCE.

    if ~isfolder(folder)
        folder = fullfile(fileparts(mfilename('fullpath')), '..', ...
                          'test_05_Accelerometer3D', 'results', folder);
    end
    D   = folder;
    logf = fullfile(D, 'Deterministic_window.txt');
    fid = fopen(logf, 'w');  fclose(fid);

    ff = dir(fullfile(D, 'FOM*.mat'));
    if isempty(ff), error('plot_deterministic_window:NoFOM', 'No FOM file in %s.', D); end
    F     = load(fullfile(D, ff(1).name), 't', 'y_contact');
    faces = fieldnames(F.y_contact);
    catn  = @(yc) cell2mat(cellfun(@(k) yc.(k).normal, faces(:), 'Uni', 0));
    Yf = catn(F.y_contact);  t = F.t(:)';  nT = numel(t);
    qtag = regexprep(ff(1).name, '^FOM[A-Z]*_(Q.*)\.mat$', '$1');

    ncc = [16 32 64 128];
    ld  = @(m, c) load(fullfile(D, sprintf('ROM_%s_Phi%d%s_%s.mat', m, phi, c, qtag)), 'y_contact');
    S0  = ld('RubinNM', '');  P = catn(S0.y_contact);
    Yir = cell(1, 4);  Ysc = cell(1, 4);
    for k = 1:4
        c = sprintf('_CC%03d', ncc(k));
        S1 = ld('RubinIRGNM', c);    Yir{k} = catn(S1.y_contact);
        S2 = ld('RubinIRGSCNM', c);  Ysc{k} = catn(S2.y_contact);
    end

    % Same GRE as the post-processing; checked on its GRE_win (first 25%).
    gre = @(Y, idx) norm(Yf(:, idx) - Y(:, idx), 'fro') / norm(Yf(:, idx), 'fro') * 100;
    iw  = 1:max(2, round(0.25*nT));
    app(logf, 'GRE_win (first 25%), to compare with GRE_Results.txt:');
    for k = 1:4
        app(logf, sprintf('  CC%3d  IR %.4f%%  IR+SC %.4f%%', ncc(k), gre(Yir{k}, iw), gre(Ysc{k}, iw)));
    end
    app(logf, sprintf('  parent %.4f%%', gre(P, iw)));

    % --- the window: onset of the parent's divergence ---
    scale = max(sqrt(sum(Yf.^2, 1))) + eps;
    dpt   = @(Y) sqrt(sum((Y - Yf).^2, 1)) / scale * 100;    % = post-processing GRE(t)
    evalc('trk = rom_tracking_time(t, dpt(P), 1, ''HoldFraction'', 0.01);');
    Th   = trk.T(1);
    cens = trk.censored(1);
    iTh  = find(t <= Th);
    if cens
        hz_txt = 'the full model never leaves the 1% band in this run: the whole run is deterministic';
    else
        hz_txt = sprintf('the full model leaves the 1%% band at %.2f \\mus', 1e6*Th);
    end
    app(logf, sprintf('\nwindow: [0, %.2f us]%s | run %.1f us', 1e6*Th, ...
        repmat(' (censored: parent never left the band)', 1, cens), 1e6*t(end)));

    g_ir = zeros(1, 4);  g_sc = zeros(1, 4);
    app(logf, sprintf('\nGRE on the window:'));
    for k = 1:4
        g_ir(k) = gre(Yir{k}, iTh);  g_sc(k) = gre(Ysc{k}, iTh);
        app(logf, sprintf('  CC%3d  IR %8.4f%%  IR+SC %8.4f%%  gain %.2fx', ...
            ncc(k), g_ir(k), g_sc(k), g_ir(k)/g_sc(k)));
    end
    g_p = gre(P, iTh);
    app(logf, sprintf('  parent %8.4f%%', g_p));

    % --- GRE over [0, T] for every T ---
    Tgrid = round(linspace(max(5, round(0.02*nT)), nT, 120));
    cum   = @(Y) arrayfun(@(j) gre(Y, 1:j), Tgrid);
    Cir = zeros(4, numel(Tgrid));  Csc = Cir;
    for k = 1:4, Cir(k, :) = cum(Yir{k});  Csc(k, :) = cum(Ysc{k}); end
    Cp = cum(P);

    % Before the first contact the variants are bit-identical (the compliance
    % acts only on an active constraint), so the ordering is checked from the
    % first window where they differ by more than 0.1%.
    inside = t(Tgrid) <= Th;
    app(logf, sprintf('\nIR+SC below IR, from the first difference:'));
    first_rev = nan(1, 4);
    for k = 1:4
        j0 = find(abs(Csc(k, :) - Cir(k, :)) > 1e-3 * Cir(k, :), 1);
        if isempty(j0), app(logf, sprintf('  CC%3d: never differ', ncc(k))); continue; end
        act = false(1, numel(Tgrid));  act(j0:end) = true;
        ok  = all(Csc(k, act & inside) < Cir(k, act & inside));
        jf  = find(act & Csc(k, :) >= Cir(k, :), 1);
        if isempty(jf), fl = 'never';
        else, first_rev(k) = t(Tgrid(jf)); fl = sprintf('%.2f us (%.2f x window)', 1e6*first_rev(k), first_rev(k)/Th); end
        app(logf, sprintf('  CC%3d: below over the whole window: %d | first reversal: %s', ncc(k), ok, fl));
    end

    res = struct('window', Th, 'censored', cens, 'ncc', ncc, 'gre_ir', g_ir, ...
                 'gre_sc', g_sc, 'gre_parent', g_p, 'first_reversal', first_rev);

    cR = [0.80 0.25 0.20];  cB = [0.10 0.40 0.75];  col = lines(4);
    tus = 1e6*t;  Tus = 1e6*t(Tgrid);  Thu = 1e6*Th;

    % ===== 1: the result =====
    f1 = figure('Color', 'w', 'Position', [80 80 760 520]);
    ax = axes(f1);  hold(ax, 'on');  box(ax, 'on');
    plot(ax, ncc, g_ir, '--o', 'Color', cR, 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    plot(ax, ncc, g_sc, '-s', 'Color', cB, 'LineWidth', 2.2, 'MarkerSize', 8, 'MarkerFaceColor', cB);
    yline(ax, g_p, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.6, ...
          'Label', sprintf('full Rubin, \\phi = %d (no IR): %.2f%%', phi, g_p), ...
          'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
    for k = 1:4
        text(ax, ncc(k), g_sc(k), sprintf('  x%.1f', g_ir(k)/g_sc(k)), 'VerticalAlignment', 'top', ...
             'Color', cB, 'FontWeight', 'bold', 'FontSize', 10);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', ncc);  grid(ax, 'on');
    xlabel(ax, 'n_{cc}  (interface DOFs kept, out of 232)');
    ylabel(ax, 'GRE on the deterministic window [%]');
    legend(ax, {'IR', 'IR + static correction'}, 'Location', 'northeast');
    title(ax, {sprintf('Rubin \\phi = %d, IR global: GRE over [0, %.2f \\mus]', phi, Thu), hz_txt});
    xlim(ax, [12 170]);  ylim(ax, [0.5*min([g_sc g_p]) 1.6*max(g_ir)]);
    ax.Toolbar.Visible = 'off';
    exportgraphics(f1, fullfile(D, 'Deterministic_1_result.png'), 'Resolution', 160);

    % ===== 2: why the window ends there =====
    f2 = figure('Color', 'w', 'Position', [80 80 900 480]);
    ax = axes(f2);  hold(ax, 'on');  box(ax, 'on');
    yl = [1e-3 3e2];
    patch(ax, [Thu tus(end) tus(end) Thu], [yl(1) yl(1) yl(2) yl(2)], [0.93 0.93 0.93], ...
          'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, tus, max(dpt(Yir{1}), 1e-6), '-', 'Color', cR, 'LineWidth', 1.0);
    plot(ax, tus, max(dpt(Ysc{1}), 1e-6), '-', 'Color', cB, 'LineWidth', 1.0);
    plot(ax, tus, max(dpt(P), 1e-6), '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.3);
    xline(ax, Thu, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    text(ax, Thu/2, yl(2)/2, {'full model tracks the FOM:', 'differences are due to the IR'}, ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
    if ~cens
        text(ax, (Thu + tus(end))/2, yl(2)/2, {'full model has left the FOM:', 'differences are amplified by the chaos'}, ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0.45 0.45 0.45]);
    end
    set(ax, 'YScale', 'log');  ylim(ax, yl);  xlim(ax, [0 tus(end)]);  grid(ax, 'on');
    xlabel(ax, 't [\mus]');  ylabel(ax, 'pointwise error vs FOM [%]');
    legend(ax, {'IR, n_{cc} = 16', 'IR + SC, n_{cc} = 16', sprintf('full Rubin, \\phi = %d', phi)}, ...
           'Location', 'southeast');
    title(ax, {'Pointwise error vs FOM, n_{cc} = 16', hz_txt});
    ax.Toolbar.Visible = 'off';
    exportgraphics(f2, fullfile(D, 'Deterministic_2_horizon.png'), 'Resolution', 160);

    % ===== 3: every window length, none chosen =====
    f3 = figure('Color', 'w', 'Position', [80 80 900 520]);
    ax = axes(f3);  hold(ax, 'on');  box(ax, 'on');
    jd   = find(any(abs(Csc - Cir) > 1e-3 * Cir, 1), 1);  if isempty(jd), jd = 1; end
    vals = [Cir(:, jd:end); Csc(:, jd:end); Cp(jd:end)];
    yl   = [0.5*min(vals(vals > 0)), 2*max(vals(:))];
    patch(ax, [Thu Tus(end) Tus(end) Thu], [yl(1) yl(1) yl(2) yl(2)], [0.93 0.93 0.93], ...
          'EdgeColor', 'none', 'HandleVisibility', 'off');
    h = gobjects(1, 9);
    for k = 1:4
        h(k)   = plot(ax, Tus, Cir(k, :), '--', 'Color', col(k, :), 'LineWidth', 1.3);
        h(4+k) = plot(ax, Tus, Csc(k, :), '-',  'Color', col(k, :), 'LineWidth', 2.0);
    end
    h(9) = plot(ax, Tus, Cp, '-', 'Color', 'k', 'LineWidth', 2.4);
    xline(ax, Thu, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    set(ax, 'YScale', 'log');  ylim(ax, yl);  xlim(ax, [0.9*Tus(jd) Tus(end)]);  grid(ax, 'on');
    xlabel(ax, 'window length T [\mus]');  ylabel(ax, 'GRE over [0, T] [%]');
    lg = arrayfun(@(c) sprintf('IR+SC, n_{cc} = %d', c), ncc, 'Uni', 0);
    legend(ax, [h(5:8) h(1) h(9)], [lg, {'dashed: IR, same colours'}, {sprintf('full Rubin, \\phi = %d', phi)}], ...
           'Location', 'northwest', 'NumColumns', 2);
    title(ax, {'GRE over [0, T] for every window length T - no window chosen', hz_txt});
    ax.Toolbar.Visible = 'off';
    exportgraphics(f3, fullfile(D, 'Deterministic_3_window_sweep.png'), 'Resolution', 160);

    % ===== 4: bars, same layout as the tracking-time summary =====
    f4 = figure('Color', 'w', 'Position', [80 80 900 520]);
    plot_window_bars(axes(f4), res, phi, Thu, hz_txt);
    exportgraphics(f4, fullfile(D, 'Deterministic_4_bars.png'), 'Resolution', 160);

    fprintf('%s', fileread(logf));
end

function app(f, s)
    fid = fopen(f, 'a');  fprintf(fid, '%s\n', s);  fclose(fid);
end

function plot_window_bars(ax, res, phi, Thu, hz_txt)
% One pair of bars per n_cc (IR, IR + static correction), the parent on its own
% at the right. Log scale because the errors span two decades; every bar carries
% its value so nothing has to be read off the axis.
    cR = [0.80 0.25 0.20];  cB = [0.10 0.40 0.75];  cK = [0.35 0.35 0.35];
    n  = numel(res.ncc);
    Y  = [res.gre_ir(:), res.gre_sc(:)];
    hold(ax, 'on');  box(ax, 'on');
    b = bar(ax, 1:n, Y, 0.8, 'EdgeColor', 'k');
    b(1).FaceColor = cR;  b(2).FaceColor = cB;
    bp = bar(ax, n + 1.3, res.gre_parent, 0.35, 'FaceColor', cK, 'EdgeColor', 'k');
    yline(ax, res.gre_parent, ':', 'Color', cK, 'LineWidth', 1.2, 'HandleVisibility', 'off');

    for j = 1:2
        for k = 1:n
            text(ax, b(j).XEndPoints(k), Y(k, j), sprintf('%.2f%%', Y(k, j)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
        end
    end
    text(ax, n + 1.3, res.gre_parent, sprintf('%.2f%%', res.gre_parent), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);

    set(ax, 'YScale', 'log', 'XTick', [1:n, n + 1.3], ...
        'XTickLabel', [arrayfun(@(c) sprintf('n_{cc} = %d', c), res.ncc, 'Uni', 0), ...
                       {sprintf('full Rubin \phi=%d', phi)}]);
    ylim(ax, [0.5*min([Y(:); res.gre_parent]), 3*max(Y(:))]);
    xlim(ax, [0.4, n + 1.8]);
    grid(ax, 'on');  ax.XGrid = 'off';
    ylabel(ax, 'GRE on the deterministic window [%]  (log scale)');
    legend(ax, [b(1) b(2) bp], {'IR', 'IR + static correction', 'no IR (reference)'}, ...
           'Location', 'northeast');
    title(ax, {sprintf('Error vs FOM over [0, %.2f \mus]  -  Rubin \phi = %d, IR global', Thu, phi), hz_txt});
    ax.Toolbar.Visible = 'off';
end
