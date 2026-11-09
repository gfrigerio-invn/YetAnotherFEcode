clear; clc; close all;

% Parametri del primo impulso
cfg1.impulse_g   = 1*1e6;        % ampiezza [g]
cfg1.impulse_dir = [0 0 1];   % direzione
cfg1.t_shock     = 1e-6;       % durata [s]

% Parametri del secondo impulso
cfg2.impulse_g   = 1e6;      % ampiezza [g]
cfg2.impulse_dir = [0 0 1];    % direzione
cfg2.t_shock     = 0.5e-6;       % durata [s]

% Direzioni normalizzate
dir1 = cfg1.impulse_dir / norm(cfg1.impulse_dir);
dir2 = cfg2.impulse_dir / norm(cfg2.impulse_dir);

% Asse temporale
Tmax = max(cfg1.t_shock, cfg2.t_shock);
Tmax = 8e-6;
t = linspace(0, 1.1*Tmax, 2000);

% Impulsi half-sine lungo la direzione z
pulse1_g = zeros(size(t));
pulse2_g = zeros(size(t));

idx1 = (t >= 0) & (t <= cfg1.t_shock);
idx2 = (t >= 0) & (t <= cfg2.t_shock);

pulse1_g(idx1) = cfg1.impulse_g * dir1(3) .* ...
                 sin(pi*t(idx1)/cfg1.t_shock);

pulse2_g(idx2) = cfg2.impulse_g * dir2(3) .* ...
                 sin(pi*t(idx2)/cfg2.t_shock);

% Conversione da g a kilo-g
% 1 kg = 1000 g, dove kg indica kilo-g
pulse1_kg = pulse1_g / 1e3;
pulse2_kg = pulse2_g / 1e3;

% Grafico
figure('Color','w');
plot(t, pulse1_kg, 'b', 'LineWidth', 2); hold on;
%plot(t, pulse2_kg, 'r--', 'LineWidth', 2);

grid on;
xlabel('Time [s]');
ylabel('Acceleration [kilo-g]');
title('Impulse Sim 5');

% legend( ...
%     sprintf('Impulse 4: %.2e g', cfg1.impulse_g), ...
%     sprintf('Impulse 3: %.2e g', cfg2.impulse_g), ...
%     'Location', 'best');

xlim([0 1.1*Tmax]);