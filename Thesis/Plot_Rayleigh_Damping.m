%% Rayleigh damping da due frequenze e due Q factor
clear;
clc;
close all;

%% Parametri utente

% Inserire le due frequenze di riferimento in Hz
f_ref = [7082.177;     % f1 [Hz]
         2e6];         % f2 [Hz]

% Inserire i due Q factor corrispondenti
Q_ref = [5;             % Q alla frequenza f1
         500];          % Q alla frequenza f2

% Damping ratio aggiuntivo costante.
% Porre dmprat = 0 per utilizzare solamente lo smorzamento di Rayleigh.
dmprat = 0.0;

% Numero di modi da visualizzare
number_of_modes = 100;

% Ipotesi per le frequenze modali.
% Sostituire questa frequenza con il vero passo modale del modello,
% oppure direttamente con il vettore delle frequenze modali.
f_mode_spacing = f_ref(1);


%% Controlli sugli input

f_ref = f_ref(:);
Q_ref = Q_ref(:);

if numel(f_ref) ~= 2 || numel(Q_ref) ~= 2
    error('Inserire esattamente due frequenze e due Q factor.');
end

if any(f_ref <= 0)
    error('Le frequenze devono essere positive.');
end

if any(Q_ref <= 0)
    error('I Q factor devono essere positivi.');
end

if dmprat < 0
    error('dmprat deve essere maggiore o uguale a zero.');
end

% Ordina le frequenze mantenendo associati i relativi Q
[f_ref, index] = sort(f_ref);
Q_ref = Q_ref(index);


%% Calcolo di alpha e beta

omega_ref = 2*pi*f_ref;

% Contributo di inverse Q da attribuire a Rayleigh damping
invQ_rayleigh_ref = 1 ./ Q_ref - 2*dmprat;

if any(invQ_rayleigh_ref < 0)
    error(['I Q specificati sono inferiori al limite imposto da dmprat. ', ...
           'Ridurre dmprat oppure modificare i Q di riferimento.']);
end

% Sistema lineare:
%
% 1/Q_i - 2*dmprat = alpha/omega_i + beta*omega_i
%
A = [1 ./ omega_ref, omega_ref];

if rcond(A) < 1e-14
    error(['Le due frequenze sono troppo vicine per determinare ', ...
           'in modo stabile alpha e beta.']);
end

coefficients = A \ invQ_rayleigh_ref;

alpha = coefficients(1);
beta  = coefficients(2);

% Con smorzamento fisico di Rayleigh alpha e beta devono essere non negativi
tolerance_alpha = 1e-12 * max(1, abs(alpha));
tolerance_beta  = 1e-12 * max(1, abs(beta));

if alpha < -tolerance_alpha || beta < -tolerance_beta
    error(['Le due coppie frequenza-Q non possono essere realizzate ', ...
           'con coefficienti di Rayleigh alpha e beta non negativi.']);
end

% Elimina eventuali piccoli errori numerici negativi
alpha = max(alpha, 0);
beta  = max(beta, 0);


%% Verifica dei Q alle due frequenze di riferimento

invQ_ref_calculated = ...
    alpha ./ omega_ref + ...
    beta  .* omega_ref + ...
    2*dmprat;

Q_ref_calculated = 1 ./ invQ_ref_calculated;


%% Frequenze modali per il grafico

mode_number = (1:number_of_modes)';
f_modes = mode_number * f_mode_spacing;
omega_modes = 2*pi*f_modes;


%% Intervallo di frequenze per il grafico continuo

f_min = min(f_ref) / 1000;
f_max = max([f_ref; f_modes]) * 1.05;

f = logspace(log10(f_min), log10(f_max), 1200);
omega = 2*pi*f;


%% Q factor in funzione della frequenza

% Contributo alpha
invQ_alpha = alpha ./ omega;
Q_alpha = 1 ./ invQ_alpha;

% Contributo beta
invQ_beta = beta .* omega;

if beta > 0
    Q_beta = 1 ./ invQ_beta;
else
    Q_beta = nan(size(f));
end

% Contributo damping ratio costante
invQ_dmprat = 2*dmprat * ones(size(f));

if dmprat > 0
    Q_dmprat = 1 ./ invQ_dmprat;
else
    Q_dmprat = nan(size(f));
end

% Q totale
invQ_total = invQ_alpha + invQ_beta + invQ_dmprat;
Q_total = 1 ./ invQ_total;


%% Q factor dei modi

invQ_alpha_modes = alpha ./ omega_modes;
Q_alpha_modes = 1 ./ invQ_alpha_modes;

invQ_beta_modes = beta .* omega_modes;

if beta > 0
    Q_beta_modes = 1 ./ invQ_beta_modes;
else
    Q_beta_modes = nan(size(f_modes));
end

invQ_dmprat_modes = 2*dmprat * ones(size(f_modes));

if dmprat > 0
    Q_dmprat_modes = 1 ./ invQ_dmprat_modes;
else
    Q_dmprat_modes = nan(size(f_modes));
end

invQ_total_modes = ...
    invQ_alpha_modes + ...
    invQ_beta_modes + ...
    invQ_dmprat_modes;

Q_modes = 1 ./ invQ_total_modes;


%% Stampa dei risultati

fprintf('\n');
fprintf('=============================================\n');
fprintf('Rayleigh damping da due punti frequenza-Q\n');
fprintf('=============================================\n');

fprintf('Frequenza 1             = %.9g Hz\n', f_ref(1));
fprintf('Q richiesto 1           = %.9g\n', Q_ref(1));
fprintf('Q calcolato 1           = %.9g\n', Q_ref_calculated(1));

fprintf('\n');

fprintf('Frequenza 2             = %.9g Hz\n', f_ref(2));
fprintf('Q richiesto 2           = %.9g\n', Q_ref(2));
fprintf('Q calcolato 2           = %.9g\n', Q_ref_calculated(2));

fprintf('\n');

fprintf('dmprat                  = %.9g\n', dmprat);
fprintf('alpha                   = %.9g 1/s\n', alpha);
fprintf('beta                    = %.9g s\n', beta);

fprintf('\n');

if beta > 0
    fprintf('Beta damping            = ENABLED\n');
else
    fprintf('Beta damping            = DISABLED\n');
end

fprintf('\n');

fprintf('Q modo 1                = %.9g\n', Q_modes(1));
fprintf('Q modo %d              = %.9g\n', ...
    number_of_modes, Q_modes(end));

fprintf('=============================================\n');
fprintf('\n');


%% Tabella di alcuni modi

selected_modes = [1; 2; 3; 4; 5; 10; 20; 50; 100];
selected_modes = selected_modes(selected_modes <= number_of_modes);

result_table = table( ...
    selected_modes, ...
    f_modes(selected_modes)/1e3, ...
    Q_alpha_modes(selected_modes), ...
    Q_beta_modes(selected_modes), ...
    Q_dmprat_modes(selected_modes), ...
    Q_modes(selected_modes), ...
    'VariableNames', { ...
    'Mode', ...
    'Frequency_kHz', ...
    'Q_alpha', ...
    'Q_beta', ...
    'Q_dmprat', ...
    'Q_total'});

disp(result_table);


%% Grafico Q in funzione della frequenza

figure('Color', 'w');

loglog(f/1e3, Q_total, 'k-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Q_{total}');
hold on;

loglog(f/1e3, Q_alpha, 'r--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Q_{\alpha}');

if beta > 0
    loglog(f/1e3, Q_beta, 'b--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Q_{\beta}');
end

if dmprat > 0
    loglog(f/1e3, Q_dmprat, 'g--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Q_{dmprat}');
end

% Evidenzia i due punti utilizzati per il fit
loglog(f_ref/1e3, Q_ref, 'mo', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'm', ...
    'DisplayName', 'Punti di riferimento');

grid on;
xlabel('Frequency (kHz)');
ylabel('Quality factor Q');
title('Rayleigh damping calibrato su due punti frequenza-Q');
legend('Location', 'best');


%% Grafico Q dei modi

figure('Color', 'w');

plot(mode_number, Q_modes, 'ko-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 4, ...
    'DisplayName', 'Q_{total}');
hold on;

plot(mode_number, Q_alpha_modes, 'r--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Q_{\alpha}');

if beta > 0
    plot(mode_number, Q_beta_modes, 'b--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Q_{\beta}');
end

if dmprat > 0
    plot(mode_number, Q_dmprat_modes, 'g--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Q_{dmprat}');
end

grid on;
xlabel('Mode number');
ylabel('Quality factor Q');
title('Quality factor dei modi');
legend('Location', 'best');


%% Grafico dei contributi in inverse Q

figure('Color', 'w');

semilogx(f/1e3, invQ_alpha, 'r-', ...
    'LineWidth', 1.5, ...
    'DisplayName', '1/Q_{\alpha}');
hold on;

if beta > 0
    semilogx(f/1e3, invQ_beta, 'b-', ...
        'LineWidth', 1.5, ...
        'DisplayName', '1/Q_{\beta}');
end

if dmprat > 0
    semilogx(f/1e3, invQ_dmprat, 'g-', ...
        'LineWidth', 1.5, ...
        'DisplayName', '1/Q_{dmprat}');
end

semilogx(f/1e3, invQ_total, 'k--', ...
    'LineWidth', 2, ...
    'DisplayName', '1/Q_{total}');

grid on;
xlabel('Frequency (kHz)');
ylabel('Inverse quality factor, 1/Q');
title('Contributi dello smorzamento');
legend('Location', 'best');