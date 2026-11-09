%% =====================================================================
%  SIMULATION MAIN - 3D accelerometer, FOM vs ROM 
%
%  Everything below the configuration is a call: the machinery lives in
%  Thesis/Src. To set up a run, edit section 1.
%
%    FeStructure       mesh, material, node sets, boundary conditions, modes
%    build_interfaces  contact operator and interface bookkeeping
%    shock_forcing     separable impulsive forcing and reference energy
%    run_fom           reference full-order runs, with the contact diagnostic
%    run_rom_sweep     every reduction method over the requested sweep
%
%  Available methods:
%    FOM     full model
%    MT      modal truncation
%    MC      Milman-Chu
%    CB      Craig-Bampton fixed-interface CMS
%    Rubin   free-interface CMS
%
%  NOTE: CB and Rubin can take a secondary interface reduction.
% =====================================================================
clear; close all; clc;

%% --- 0. PATHS ---------------------------------------------------------
test_root = fileparts(mfilename('fullpath'));
addpath(fullfile(test_root, '..', 'Src'));
addpath(fullfile(test_root, 'mesh'));

% YaFEc startup.
if isempty(which('KirchoffMaterial'))
    run(fullfile(test_root, '..', '..', 'startup.m'));
end

cd(test_root);

%% --- 1. CONFIGURATION -------------------------------------------------

% --- Mesh ---
cfg.mesh_file    = 'Accel3D_V1_F.mat';   % see mesh/README.md
cfg.element_type = 'WED15';

% --- Material ---
% UNITS: mesh works in the consistent system um / MPa / kg.
cfg.E   = 168e3;      % [MPa]     
cfg.nu  = 0.23;
cfg.rho = 2.33e-15;   % [kg/um^3] 

% --- Node sets ---
z_top = 30.3;              % upper face of the proof mass
z_bot = 0.0;               % lower face
z_anc = [52.3, -2.18];     % the two anchoring levels
tol   = 1e-3;              % tolerance on the coordinates

tab = { 'A', [-202.875 -172.875], [-115.000 -105.000]
        'B', [ 172.875  202.875], [-115.000 -105.000]
        'C', [-151.350 -121.350], [ 139.825  149.825]
        'D', [ 162.875  192.875], [ 139.825  149.825] };

box = @(x, y, z) [x(1)-tol x(2)+tol; y(1)-tol y(2)+tol; z-tol z+tol];
for it = 1:size(tab, 1)
    cfg.node_sets.(['stop_' tab{it,1} '_top']) = ...
        struct('box', box(tab{it,2}, tab{it,3}, z_top));
    cfg.node_sets.(['stop_' tab{it,1} '_bot']) = ...
        struct('box', box(tab{it,2}, tab{it,3}, z_bot));
end

cfg.node_sets.anchor = [struct('box', [-inf inf; -inf inf; z_anc(1)-tol z_anc(1)+tol]), ...
                        struct('box', [-inf inf; -inf inf; z_anc(2)-tol z_anc(2)+tol])];
cfg.bc_sets = {'anchor'};

% --- Contact interfaces ---
% For a wall that is NOT parallel to the contacting surface, replace 'gap' with 
%     cfg.interfaces(1).plane = struct('point', [x0 y0 z0]);
gap_top = 4.0;            % [um]
gap_bot = 2.0;            % [um]

k = 0;
for it = 1:size(tab, 1)
    for face = {'top', [0 0 1], gap_top; 'bot', [0 0 -1], gap_bot}'
        k = k + 1;
        cfg.interfaces(k).set    = ['stop_' tab{it,1} '_' face{1}];
        cfg.interfaces(k).normal = face{2};
        cfg.interfaces(k).gap    = face{3};
        cfg.interfaces(k).dofs   = 'normal';
    end
end

% --- Methods to run ---
cfg.run.FOM   = 0;
cfg.run.MT    = 0;
cfg.run.MC    = 0;
cfg.run.CB    = 0;
cfg.run.Rubin = 1;
cfg.run.MCB   = 0;  
cfg.run.MN    = 0;

% --- Interface reduction ---
cfg.interface_reduction.enabled = 1;
cfg.interface_reduction.skip_baseline = 1;
cfg.interface_reduction.mode    = 'per_interface';   % 'global' | 'per_interface' 
cfg.interface_reduction.static_correction = true;
cfg.array_ccModes  = [16];

% --- Sweep ---
cfg.array_linModes = [1];   

% --- Damping ---
% Rayleigh, C = alpha*M + beta*K,  1/Q(f) = alpha/(2*pi*f) + 2*pi*f*beta.
cfg.Q_freq         = [7e3, 2e6];    % anchor frequencies [Hz]
cfg.array_QFactor  = [5,  200];   % Q at anchor frequencies
cfg.array_k_mult   = 0.003125;  % contact stiffness as a multiple of max(diag(K)).
                            
% --- Impulsive forcing ---
cfg.impulse_g   = 1e6;        % amplitude [g]
cfg.impulse_dir = [0 0 -1];   % direction in space, normalized afterwards
cfg.g_value     = 9.81e6;     % gravity [um/s^2]
cfg.t_shock     = 1e-6;      % half-sine duration [s]

% --- Integration ---
cfg.tmax      = 2e-6;
cfg.output_stride = 10;       % output every N steps of dt
cfg.integrator    = 'genalpha';       % 'ode15s' | 'genalpha' | 'newmark'
cfg.dt         = 0.5e-9;           % fixed step of the second-order schemes
cfg.rho_inf       = 0.7;            % genalpha only
cfg.newmark_alpha = 0;              % newmark only
cfg.RelTol    = 1e-8;         % ode15s only
cfg.RelTolFOM = 1e-8;         % ode15s only
% Store the full displacement field of the FOM.
cfg.save_full_field = true;

cfg.test_name = sprintf('Shock3D_%dg_%.1es', cfg.impulse_g, cfg.tmax);

%% --- 2. RESULTS DIRECTORY ---------------------------------------------
run_dir = fullfile('results', sprintf('%s_%s', cfg.test_name, ...
    char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'))));
if ~exist(run_dir, 'dir'), mkdir(run_dir); end
fprintf('Results directory: %s\n\n', run_dir);

%% --- 3. MODEL ---------------------------------------------------------
% --- Build Structure ---
Struct = FeStructure();
Struct.mesh_file    = cfg.mesh_file;
Struct.element_type = cfg.element_type;
Struct.E            = cfg.E;
Struct.nu           = cfg.nu;
Struct.rho          = cfg.rho;
Struct.set_specs    = cfg.node_sets;
Struct.bc_sets      = cfg.bc_sets;
Struct.build();
Struct.describe_node_sets();

% --- Modes ---
Struct.compute_eigenmodes(max(cfg.array_linModes));

%% --- 4. CONTACT AND FORCING -------------------------------------------
contact = build_interfaces(Struct, cfg);
shock   = shock_forcing(Struct, cfg);

%% --- 5. RUN -----------------------------------------------------------
save_run_config(run_dir, cfg, Struct, contact, shock);

if cfg.run.FOM
    run_fom(Struct, contact, shock, cfg, run_dir);
end

run_rom_sweep(Struct, contact, shock, cfg, run_dir);

fprintf('\n=========================================\n');
fprintf('  Benchmark complete.\n  Results in: %s\n', run_dir);
fprintf('=========================================\n');
