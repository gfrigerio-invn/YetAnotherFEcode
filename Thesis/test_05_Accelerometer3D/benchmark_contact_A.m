%% BENCHMARK STANDALONE DELLA COSTRUZIONE DI A
%
% Confronta:
%
%   1) Versione originale:
%
%      A = diag(1 ./ k_pen(act)) + Scomp(act, act);
%
%   2) Versione ottimizzata:
%
%      idx = find(act);
%      A = Scomp(idx, idx);
%      A(diag_idx) = A(diag_idx) + d;
%
%   3) Versione ottimizzata con idx, d e diag_idx precalcolati.
%
% Questo script genera una Scomp densa sintetica perché quella reale
% non è salvata su disco.

clear;
clc;

%% ================================================================
% PARAMETRI DEL BENCHMARK
% =================================================================

% Per un test sicuro iniziare con 3000 o 4000.
% Per simulare il tuo caso usare n_c = 11240, se hai memoria sufficiente.
n_c = 4000;

% Frazione di vincoli attivi
active_fraction = 0.50;

% Numero di ripetizioni per metodo
nRep = 3;

% Eseguire anche la versione originale?
%
% La versione originale può richiedere molta memoria perché crea:
%
%   diag(...)
%   Scomp(act,act)
%   A
%
run_original = true;

% Controllo numerico tra originale e ottimizzata.
% Per matrici grandi può richiedere memoria aggiuntiva.
check_result = true;

% Generazione della matrice densa sintetica.
%
% Se true usa randn, che è più costoso ma genera dati generici.
% Se false usa rand, leggermente più veloce.
use_randn = false;

rng(12345);


%% ================================================================
% DIMENSIONI E MEMORIA STIMATA
% =================================================================

n_a = round(active_fraction * n_c);

if n_a < 1
    error('Il numero di vincoli attivi deve essere almeno 1.');
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('BENCHMARK STANDALONE COSTRUZIONE A\n');
fprintf('============================================================\n');
fprintf('Numero totale vincoli       : %d\n', n_c);
fprintf('Numero vincoli attivi       : %d\n', n_a);
fprintf('Frazione attiva             : %.2f\n', n_a / n_c);
fprintf('Numero ripetizioni          : %d\n', nRep);
fprintf('\n');

mem_Scomp_GB = 8 * double(n_c)^2 / 1024^3;
mem_A_GB     = 8 * double(n_a)^2 / 1024^3;

fprintf('Memoria stimata Scomp       : %.3f GB\n', mem_Scomp_GB);
fprintf('Memoria stimata A attiva    : %.3f GB\n', mem_A_GB);
fprintf('\n');

if n_c >= 10000
    warning(['Stai per creare una matrice densa molto grande. ', ...
             'Verifica la memoria disponibile prima di continuare.']);
end


%% ================================================================
% CREAZIONE DEI DATI SINTETICI
% =================================================================

fprintf('Genero Scomp densa %d x %d...\n', n_c, n_c);

tic;

if use_randn
    Scomp = randn(n_c, n_c);
else
    Scomp = rand(n_c, n_c);
end

t_generation = toc;

fprintf('Tempo generazione Scomp      : %.6f s\n', t_generation);
fprintf('\n');


% Active set casuale
act = false(n_c, 1);

active_idx_random = randperm(n_c, n_a);
act(active_idx_random) = true;

% Penalità positive
k_pen = 0.5 + rand(n_c, 1);


%% ================================================================
% PRECALCOLO DELLA VERSIONE CACHED
% =================================================================

idx = find(act);
d = 1 ./ k_pen(idx);

diag_idx = 1:(n_a + 1):(n_a * n_a);


%% ================================================================
% FUNCTION HANDLE
% =================================================================

% Versione originale
f_original = @() build_A_original(act, k_pen, Scomp);

% Versione ottimizzata: calcola idx e d ogni volta
f_optimized = @() build_A_optimized(act, k_pen, Scomp);

% Versione ottimizzata cached: idx, d e diag_idx già calcolati
f_cached = @() build_A_cached(idx, d, diag_idx, Scomp);


%% ================================================================
% BENCHMARK
% =================================================================

method_names = {};
times = [];


% ----------------------------------------------------------------
% Metodo 1: versione originale
% ----------------------------------------------------------------

if run_original

    fprintf('Misuro versione originale...\n');

    t_original = measure_function(f_original, nRep);

    method_names{end+1,1} = 'originale';
    times(end+1,1) = t_original;

    fprintf('Tempo mediano originale     : %.6f s\n', t_original);
    fprintf('\n');
end


% ----------------------------------------------------------------
% Metodo 2: ottimizzata con find
% ----------------------------------------------------------------

fprintf('Misuro versione ottimizzata con find...\n');

t_optimized = measure_function(f_optimized, nRep);

method_names{end+1,1} = 'ottimizzata_find';
times(end+1,1) = t_optimized;

fprintf('Tempo mediano ottimizzata   : %.6f s\n', t_optimized);
fprintf('\n');


% ----------------------------------------------------------------
% Metodo 3: ottimizzata cached
% ----------------------------------------------------------------

fprintf('Misuro versione ottimizzata cached...\n');

t_cached = measure_function(f_cached, nRep);

method_names{end+1,1} = 'ottimizzata_cached';
times(end+1,1) = t_cached;

fprintf('Tempo mediano cached         : %.6f s\n', t_cached);
fprintf('\n');


%% ================================================================
% RISULTATI
% =================================================================

if run_original
    speedup = times(1) ./ times;
else
    speedup = nan(size(times));
end

results = table( ...
    method_names, ...
    times, ...
    speedup, ...
    'VariableNames', { ...
        'metodo', ...
        'tempo_secondi', ...
        'speedup_vs_originale'});

fprintf('============================================================\n');
fprintf('RISULTATI\n');
fprintf('============================================================\n');

disp(results);


%% ================================================================
% CONTROLLO NUMERICO
% =================================================================

if check_result && run_original

    fprintf('Controllo numerico...\n');

    % Creo una matrice alla volta e poi confronto
    A_original = f_original();
    A_optimized = f_optimized();

    difference = A_original - A_optimized;

    err_abs = max(abs(difference(:)));

    scale = max(abs(A_original(:)));

    if scale == 0
        scale = 1;
    end

    err_rel = err_abs / scale;

    fprintf('Errore assoluto massimo     : %.6e\n', err_abs);
    fprintf('Errore relativo massimo     : %.6e\n', err_rel);

    clear A_original A_optimized difference;
end


%% ================================================================
% FINE
% =================================================================

fprintf('============================================================\n');
fprintf('Benchmark terminato.\n');
fprintf('============================================================\n');


%% =================================================================
% FUNZIONI LOCALI
% =================================================================

function A = build_A_original(act, k_pen, Scomp)
% Versione originale

    A = diag(1 ./ k_pen(act)) + Scomp(act, act);

end


function A = build_A_optimized(act, k_pen, Scomp)
% Versione ottimizzata.
%
% idx e d vengono calcolati all'interno della funzione.

    idx = find(act);

    d = 1 ./ k_pen(idx);

    n_a = numel(idx);

    A = Scomp(idx, idx);

    diag_idx = 1:(n_a + 1):(n_a * n_a);

    A(diag_idx) = A(diag_idx) + d;

end


function A = build_A_cached(idx, d, diag_idx, Scomp)
% Versione ottimizzata con dati precalcolati.

    A = Scomp(idx, idx);

    A(diag_idx) = A(diag_idx) + d;

end


function t = measure_function(fun, nRep)
% Misura il tempo mediano di esecuzione.

    % Warm-up JIT
    A = fun(); %#ok<NASGU>
    clear A;

    times = zeros(nRep, 1);

    for i = 1:nRep

        tic;

        A = fun(); %#ok<NASGU>

        times(i) = toc;

        clear A;
    end

    t = median(times);

end