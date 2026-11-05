%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Made By:                        %
%            Pablo Venicius Costa Guimarães            %
%        Discente de Engenharia de Energia-UFPA        %
%                                                      %
%  VERSÃO REVISADA — BANCA:                            %
%    - Modos de controle por VELOCIDADE do vento      %
%        Cut-in: V0 < 4 m/s                           %
%        MPPT:   4 ≤ V0 < 11 m/s                      %
%        MPTA:   V0 ≥ 11 m/s                          %
%    - Gerador PMSG redimensionado para o novo        %
%      envelope de potência (~1160 W)                 %
%    - MPPT reformulado para rotação fixa             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========== WIND SIMULATION WITH ADVANCED BEMT AND GEARBOX ========== %%
clc; clear; close all;

%% ========== GLOBAL PARAMETERS ========== %%
tol = 1e-5;
maxIter = 10000;
min_r_ratio = 0.1379;
R = 1.37860;
relax = 0.2;
B = 4;
rho = 1.225;
mu = 1.8e-5;
Re_ref = 1e4;
use_Re_scaling = true;
limit_high_Ct = false;

% FIXED ROTOR RPM (design constraint preserved)
RPM_rotor = 260;
Omega_rotor = pi * RPM_rotor / 30;

% GEARBOX PARAMETERS
gear_ratio = 2;
efficiency_gearbox = 0.97;

% ========== PMSG REDIMENSIONADO (novo envelope) ==========
P = 6;
Rs       = 0.1;      % was 0.2
Ld       = 0.005;    % was 0.003
Lq       = 0.005;    % was 0.003
lambda_m = 0.15;     % was 0.1
Vdc      = 120;      % was 48 — 120 V evita enfraquecimento de campo

% ========== PONTOS DE REFERÊNCIA ==========
P_nominal = 1160;      % Novo nominal elétrico (em V0 = 11 m/s) [W]
V0_nominal = 6.63;     % Ponto de projeto aerodinâmico (pico Cp) — mantido em plots

% Thermal model parameters
T_amb = 25;
T_winding = T_amb;
R_th = 1.5;
C_th = 5000;
alpha_cu = 0.00393;

% Enhanced loss model parameters
K_h = 0.005;
K_e = 3e-6;
K_rot = 1e-5;

% Fatigue parameters (mantidos — recalibráveis para o novo gerador)
k_stress = 1.15;
m = 8.5;
S_N = 25e6;
base_mechanical_stress = 2.0e6;
base_thermal_stress = 0.12e6;
stress_min_detectable = 0.05e6;
damage_threshold = 0.008;

% ========== CONTROL MODE THRESHOLDS (REVISED BY BANCA) ==========
V_cut_in_to_mppt = 4.0;
V_mppt_to_mpta   = 11.0;
V_hysteresis     = 0.0;

P_cut_in_threshold = 20;
current_smoothing_factor = 0.2;

fprintf('Calibrated Fatigue Parameters:\n');
fprintf('  Base Mechanical Stress: %.1f MPa\n', base_mechanical_stress/1e6);
fprintf('  Base Thermal Stress: %.1f MPa/°C\n', base_thermal_stress/1e6);
fprintf('  Fatigue Strength: %.0f MPa\n', S_N/1e6);

fprintf('\n=== NEW CONTROL MODE THRESHOLDS ===\n');
fprintf('  Cut-in → MPPT: %.1f m/s\n', V_cut_in_to_mppt);
fprintf('  MPPT → MPTA:   %.1f m/s\n', V_mppt_to_mpta);
fprintf('  Hysteresis:    ±%.1f m/s\n', V_hysteresis);

%% ========== AUTOMATED PARAMETER CALCULATION ========== %%
eta_gen_est = 0.92;

RPM_generator = RPM_rotor * gear_ratio;
Omega_generator = Omega_rotor * gear_ratio;
Omega_nominal = Omega_generator;

T_nominal = (P_nominal / eta_gen_est) / Omega_generator;
I_nominal = (2/3) * T_nominal / (P * lambda_m);

I_max_continuous = I_nominal * 1.2;
I_max_peak = I_nominal * 1.5;

V_ph_nominal = Vdc / (2 * sqrt(3));
V_ph_max = Vdc * 0.9 / sqrt(3);

T_max_winding = T_amb + 130;
T_derating_start = T_amb + 110;

fprintf('Automated Parameter Calculation (NEW generator):\n');
fprintf('  Nominal Power:           %.0f W\n', P_nominal);
fprintf('  Nominal Torque:          %.2f Nm\n', T_nominal);
fprintf('  Nominal Current:         %.2f A\n', I_nominal);
fprintf('  Max Continuous Current:  %.2f A\n', I_max_continuous);
fprintf('  Max Peak Current:        %.2f A\n', I_max_peak);
fprintf('  Max Phase Voltage:       %.2f V\n', V_ph_max);
fprintf('  Max Winding Temperature: %.0f°C\n', T_max_winding);

fprintf('Gearbox Configuration:\n');
fprintf('  Rotor RPM: %.1f\n', RPM_rotor);
fprintf('  Generator RPM: %.1f (Ratio: %.1f:1)\n', RPM_generator, gear_ratio);
fprintf('  Mechanical Efficiency: %.1f%%\n', efficiency_gearbox*100);

% Wind speed range
NS = 49;
Vi = 1;
Vf = 13;
V0_range = linspace(Vi, Vf, NS);
A = pi * R^2;

%% ========== BLADE GEOMETRY INPUT ========== %%
fid = fopen('shape.txt');
if fid == -1
    error('shape.txt file not found!');
end

blade_data = [];
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || line(1) == '%'
        continue;
    end
    data_line = sscanf(line, '%f %f %f');
    if numel(data_line) == 3
        blade_data = [blade_data; data_line'];
    end
end
fclose(fid);

if isempty(blade_data)
    error('No valid data found in shape.txt');
end

r = blade_data(:,1);
chord = blade_data(:,2);
beta_rad = blade_data(:,3);
beta_deg = rad2deg(beta_rad);
beta_deg_base = beta_deg;

assert(all(r >= 0), 'Negative radial positions detected!');
assert(all(chord > -1), 'Non-positive chords detected!');

if length(r) > 1
    mid_points = (r(1:end-1) + r(2:end)) / 2;
    boundaries = [0; mid_points; r(end)];
    dr = diff(boundaries);
else
    boundaries = [0; r(end)];
    dr = diff(boundaries);
end

%% ========== AERODYNAMIC DATA INPUT ========== %%
if exist('NACA653618CLCD.txt', 'file') ~= 2
    error('NACA653618CLCD.txt file not found!');
end

polar_data = load('NACA653618CLCD.txt');
polar.alpha = polar_data(:,1);
polar.CL = polar_data(:,2);
polar.CD = polar_data(:,3);

[CL_max, idx_stall_pos] = max(polar.CL);
stall_alpha_pos = polar.alpha(idx_stall_pos);
[CL_min, idx_stall_neg] = min(polar.CL);
stall_alpha_neg = polar.alpha(idx_stall_neg);
CD_min = min(polar.CD);

F_CL = griddedInterpolant(polar.alpha, polar.CL, 'pchip', 'linear');
F_CD = griddedInterpolant(polar.alpha, polar.CD, 'pchip', 'linear');

viterna_calc = @(alpha_deg, CL_max) [
    (CL_max/2) * sind(2*alpha_deg),
    CL_max * (sind(alpha_deg))^2
];

%% ========== CORRECTED AERODYNAMIC MODELS ==========
shen_tip_correction = @(g, a) (2/pi)*acos(exp(-g .* (1 + (0.1*B)./(1 - 2*a))));

tau_ref = 4 * mean(chord) / 10;
CL_old = zeros(size(r));
CD_old = zeros(size(r));
CL_prev = zeros(size(r));
CD_prev = zeros(size(r));

%% ========== GLOBAL OUTPUT STRUCTURE ==========
global_output = struct();
global_output.V0 = V0_range;
global_output.TSR = zeros(1, NS);
global_output.Cp = zeros(1, NS);
global_output.Cq = zeros(1, NS);
global_output.Ct = zeros(1, NS);
global_output.P_available = zeros(1, NS);
global_output.P_rotor = zeros(1, NS);
global_output.P_mechanical = zeros(1, NS);
global_output.Q_rotor = zeros(1, NS);
global_output.iterations = zeros(1, NS);
global_output.converged = false(1, NS);
global_output.convergence_failed = false(1, NS);
global_output.final_errors = zeros(NS, 6);

num_radial_points = length(r);
global_output.a = zeros(NS, num_radial_points);
global_output.ap = zeros(NS, num_radial_points);
global_output.alpha = zeros(NS, num_radial_points);
global_output.phi = zeros(NS, num_radial_points);
global_output.F = zeros(NS, num_radial_points);
global_output.CL = zeros(NS, num_radial_points);
global_output.CD = zeros(NS, num_radial_points);
global_output.efficiency = zeros(NS, num_radial_points);
global_output.Re_local = zeros(NS, num_radial_points);

global_output.gear_ratio = gear_ratio;
global_output.efficiency_gearbox = efficiency_gearbox;
global_output.RPM_generator = RPM_generator;
global_output.Omega_generator = Omega_generator;
global_output.Q_generator = zeros(1, NS);
global_output.P_generator = zeros(1, NS);
global_output.damage = zeros(1, NS);
global_output.stress_range = zeros(1, NS);

global_output.P_electrical = zeros(1, NS);
global_output.I_phase = zeros(1, NS);
global_output.V_phase = zeros(1, NS);
global_output.T_em = zeros(1, NS);
global_output.P_copper = zeros(1, NS);
global_output.P_iron = zeros(1, NS);
global_output.id_ref = zeros(1, NS);
global_output.iq_ref = zeros(1, NS);
global_output.control_mode = cell(1, NS);
global_output.efficiency_gen = zeros(1, NS);
global_output.T_winding = zeros(1, NS);

%% ========== MAIN LOOP OVER WIND SPEEDS ========== %%
t_start_global = tic;

last_control_mode = 'CUT-IN';
P_mech_shaft_prev = 0;
mode_persistence = 0;
cut_in_filtered = 0;
cut_in_detected = false;
id_smoothed = 0;
iq_smoothed = 0;

for idx_v0 = 1:NS
    if idx_v0 == 1
        fprintf('Completed 0.0%%, ETA: Calculating...\n');
    else
        elapsed_time = toc(t_start_global);
        avg_time_per_point = elapsed_time / (idx_v0-1);
        remaining_time = (NS - idx_v0) * avg_time_per_point;
        fprintf('Completed %.1f%%, ETA: %s\n', ...
                100*(idx_v0-1)/NS, datestr(now + remaining_time/86400, 'dd-mmm-yyyy HH:MM:SS'));
    end

    V0 = V0_range(idx_v0);
    fprintf('\n\n======= PROCESSING V0 = %.2f m/s (%d/%d) =======\n', V0, idx_v0, NS);

    if V0 < 1
        pitch_offset = -5 * (1 - V0/3);
        beta_deg = beta_deg_base + pitch_offset;
    else
        beta_deg = beta_deg_base;
    end

    TSR = Omega_rotor * R / V0;
    fprintf('  Rotor RPM: %.1f, TSR = %.2f\n', RPM_rotor, TSR);

    converged = false;
    relax_current = 0.2;

    if V0 < 3
        if idx_v0 == 1 || ~converged
            a = 0.05 * ones(size(r));
            ap = 0.02 * ones(size(r));
        end
        relax_current = 0.15;
        current_tol = tol * 25;
        maxIter_local = maxIter * 1.5;
    elseif V0 > 10
        relax_current = 0.25;
        current_tol = tol * 3;
        maxIter_local = maxIter;
    else
        relax_current = 0.2;
        current_tol = tol;
        maxIter_local = maxIter;
    end

    if idx_v0 > 1
        a = global_output.a(idx_v0-1, :);
        ap = global_output.ap(idx_v0-1, :);
    else
        if TSR > 20
            a = 0.01 * ones(size(r));
            ap = 0.01 * ones(size(r));
        else
            a = zeros(size(r)) + 0.3;
            ap = zeros(size(r)) + 0.01;
        end
    end

    f_stall = ones(size(r));
    CL_prev = zeros(size(r));
    CD_prev = zeros(size(r));

    for iter = 1:maxIter_local
        a_old = a;
        ap_old = ap;
        CL_prev = CL_old;
        CD_prev = CD_old;

        for i = 1:length(r)
            if r(i) < min_r_ratio * R, continue; end

            Vax = V0 * (1 - a(i));
            Vtan = Omega_rotor * r(i) * (1 + ap(i));

            if Vax < 0 || Vtan < 0
                a(i) = 1; ap(i) = 0; continue;
            end

            phi_i = atan2(Vax, Vtan);
            alpha_i = rad2deg(phi_i) - beta_deg(i);
            alpha_i = max(-15, min(15, alpha_i));

            g_tip = (B/2) * (R - r(i)) / (r(i) * abs(sin(phi_i)));
            g_root = (B/2) * (r(i) - min_r_ratio*R) / (r(i) * abs(sin(phi_i)));

            if a(i) < 0.5
                factor_tip = 1 + (0.1*B)/(1 - 2*a(i));
                factor_root = 1 + (0.1*B)/(1 - 2*a(i));
            else
                factor_tip = 1; factor_root = 1;
            end

            F_tip = (2/pi) * acos(exp(-g_tip * factor_tip));
            F_root = (2/pi) * acos(exp(-g_root * factor_root));
            F_i = F_tip * F_root;

            W = sqrt(Vax^2 + Vtan^2);
            Re_local = (rho * W * chord(i)) / mu;

            if use_Re_scaling
                CL_max_adj = CL_max * (Re_local / Re_ref)^0.1;
                CL_min_adj = CL_min * (Re_local / Re_ref)^0.1;
            else
                CL_max_adj = CL_max;
                CL_min_adj = CL_min;
            end

            if alpha_i > 15
                viterna_vals = viterna_calc(alpha_i, CL_max_adj);
                CL_i = viterna_vals(1); CD_i = viterna_vals(2);
                f_stall(i) = max(0.1, min(f_stall(i) + 0.02 * (alpha_i - stall_alpha_pos), 1.0));
                CL_i = f_stall(i) * F_CL(15) + (1 - f_stall(i)) * CL_i;
                CD_i = f_stall(i) * F_CD(15) + (1 - f_stall(i)) * CD_i;
            elseif alpha_i < -15
                viterna_vals = viterna_calc(alpha_i, abs(CL_min_adj));
                CL_i = viterna_vals(1); CD_i = viterna_vals(2);
                f_stall(i) = max(0.1, min(f_stall(i) + 0.02 * (abs(alpha_i) - abs(stall_alpha_neg)), 1.0));
                CL_i = f_stall(i) * F_CL(-15) + (1 - f_stall(i)) * CL_i;
                CD_i = f_stall(i) * F_CD(-15) + (1 - f_stall(i)) * CD_i;
            else
                CL_i = F_CL(alpha_i); CD_i = F_CD(alpha_i);
            end

            if iter > 1
                dt = 1/(Omega_rotor/(2*pi));
                tau = 4 * chord(i) / W;
                k = exp(-dt/tau);
                CL_i = k * CL_old(i) + (1 - k) * CL_i;
                CD_i = k * CD_old(i) + (1 - k) * CD_i;
            end
            CL_old(i) = CL_i; CD_old(i) = CD_i;

            sigma = B * chord(i) / (2 * pi * r(i));
            sin_phi = sin(phi_i); cos_phi = cos(phi_i);
            Cn = CL_i * cos_phi + CD_i * sin_phi;
            Ct_force = CL_i * sin_phi - CD_i * cos_phi;

            K = 4 * F_i * sin_phi^2 / (sigma * Cn);
            a_new = 1 / (1 + K);

            CT_val = 4 * F_i * a_new * (1 - a_new);
            CT_crit = 0.96 * F_i;
            if CT_val > CT_crit
                term = CT_val*(50 - 36*F_i) + 12*F_i*(3*F_i - 4);
                if term < 0
                    a_new = 0.5;
                else
                    a_new = (18*F_i - 20 - 3*sqrt(term)) / (36*F_i - 50);
                end
                a_new = max(0.1, min(0.5, a_new));
            end

            if abs(sin_phi) > 1e-3 && abs(Ct_force) > 1e-6
                Kt = 4 * F_i * sin_phi * cos_phi / (sigma * Ct_force);
                ap_new = 1 / (Kt - 1);
            else
                ap_new = 0;
            end

            a(i) = max(0, min(0.5, (1 - relax_current) * a_old(i) + relax_current * a_new));
            ap(i) = max(-0.5, min(0.5, (1 - relax_current) * ap_old(i) + relax_current * ap_new));
        end

        a_abs_diff = max(abs(a - a_old));
        ap_abs_diff = max(abs(ap - ap_old));
        CL_abs_diff = max(abs(CL_old - CL_prev));
        CD_abs_diff = max(abs(CD_old - CD_prev));
        a_rel_diff = max(abs((a - a_old) ./ max(abs(a), 1e-6)));
        ap_rel_diff = max(abs((ap - ap_old) ./ max(abs(ap), 1e-6)));

        convergence_metrics = [a_abs_diff, ap_abs_diff, CL_abs_diff, CD_abs_diff, a_rel_diff, ap_rel_diff];
        convergence_weights = [0.4, 0.4, 0.1, 0.1, 0.05, 0.05];
        weighted_error = sum(convergence_metrics .* convergence_weights);

        if weighted_error < current_tol || ...
           (iter > maxIter_local*0.6 && max(a_abs_diff, ap_abs_diff) < 1e-3)
            converged = true; break;
        end

        if iter == 1
            relax_current = 0.3;
        elseif max(a_abs_diff) > 0.1
            relax_current = 0.1;
        elseif iter > 20 && max(a_abs_diff) < 1e-3
            relax_current = 0.4;
        else
            relax_current = 0.2;
        end
    end

    if ~converged
        warning('BEMT did not converge at V0=%.2f m/s after %d iterations.', V0, iter);
        global_output.convergence_failed(idx_v0) = true;
        global_output.final_errors(idx_v0, :) = convergence_metrics;
    end

    [Q_rotor, T_total, alpha_local, phi_local, F_local, CL_local, CD_local, eff_local, Re_local_vec] = ...
        calculate_performance(r, a, ap, V0, Omega_rotor, beta_deg, B, R, min_r_ratio, ...
        F_CL, F_CD, stall_alpha_pos, stall_alpha_neg, f_stall, dr, chord, polar, mu, rho, use_Re_scaling, Re_ref);

    P_rotor = Q_rotor * Omega_rotor;

    if P_rotor < -1e-3
        fprintf('  Propeller brake state detected at V0=%.2f m/s - Applying brake\n', V0);
        P_rotor = 0; Q_rotor = 0; T_total = 0;
        for i = 1:length(r)
            if r(i) >= min_r_ratio * R
                a(i) = 0.0; ap(i) = 0.0;
            end
        end
        alpha_local(:) = 0; CL_local(:) = 0; CD_local(:) = 0; eff_local(:) = 0;
    end

    P_available = 0.5 * rho * V0^3 * A;

    if P_available > 0
        Cp = P_rotor / P_available;
        Cq = Q_rotor / (0.5 * rho * V0^2 * A * R);
        Ct_global = T_total / (0.5 * rho * V0^2 * A);
    else
        Cp = 0; Cq = 0; Ct_global = 0;
    end

    P_mechanical = P_rotor * efficiency_gearbox;
    Q_generator = P_mechanical / Omega_generator;
    P_generator = P_mechanical;

    power_error = abs(P_mechanical - (Q_generator * Omega_generator));
    if power_error > 1e-6
        Q_generator = P_mechanical / Omega_generator;
    end

    if any(a < -0.1 | a > 1.1), warning('Induction factor out of bounds at V0=%.2f', V0); end
    if any(ap < -0.5 | ap > 0.5), warning('Tangential induction out of bounds at V0=%.2f', V0); end
    if P_rotor < -1e-3, warning('Negative rotor power: %.6f W at V0=%.2f', P_rotor, V0); end
    if P_rotor > P_available * 1.1 && P_available > 1
        warning('Rotor power exceeds available power at V0=%.2f', V0);
    end
    if Cp > 0.593 * 1.1 && P_available > 1
        warning('Cp exceeds Betz limit with tolerance: %.4f at V0=%.2f', Cp, V0);
    end

    global_output.TSR(idx_v0) = TSR;
    global_output.Cp(idx_v0) = Cp;
    global_output.Cq(idx_v0) = Cq;
    global_output.Ct(idx_v0) = Ct_global;
    global_output.P_available(idx_v0) = P_available;
    global_output.P_rotor(idx_v0) = P_rotor;
    global_output.P_mechanical(idx_v0) = P_mechanical;
    global_output.Q_rotor(idx_v0) = Q_rotor;
    global_output.Q_generator(idx_v0) = Q_generator;
    global_output.P_generator(idx_v0) = P_generator;
    global_output.iterations(idx_v0) = iter;
    global_output.converged(idx_v0) = converged;

    global_output.a(idx_v0, :) = a;
    global_output.ap(idx_v0, :) = ap;
    global_output.alpha(idx_v0, :) = alpha_local;
    global_output.phi(idx_v0, :) = phi_local;
    global_output.F(idx_v0, :) = F_local;
    global_output.CL(idx_v0, :) = CL_local;
    global_output.CD(idx_v0, :) = CD_local;
    global_output.efficiency(idx_v0, :) = eff_local;
    global_output.Re_local(idx_v0, :) = Re_local_vec;

    fprintf('  V0 = %.2f m/s: TSR = %.2f, Cp = %.4f\n', V0, TSR, Cp);
    fprintf('  Powers: Available = %.1f W, Rotor = %.1f W, Generator = %.1f W\n', ...
            P_available, P_mechanical, P_generator);
end

%% ========== AUTOMATED CUT-IN DETECTION ========== %%
P_min_viable = 0.02 * P_nominal;
min_TSR_operational = 2.0;   % relaxed (was 3.0) — fixed-RPM turbine works down to TSR~2.3

J_rotor = 0.5 * rho * A * R^3 * 0.1;
J_generator = 0.05;
J_total = J_rotor + J_generator / gear_ratio^2;

Omega_cut_in_min = (2 * pi * V0_range(3)) / (min_TSR_operational * R);

[max_Cp, max_Cp_idx] = max(global_output.Cp);
optimal_TSR = global_output.TSR(max_Cp_idx);
K_opt = (0.5 * rho * A * R^3 * max_Cp) / (optimal_TSR^3);

fprintf('\nAutomated Cut-in Calculation:\n');
fprintf('  Minimum Viable Power: %.1f W\n', P_min_viable);
fprintf('  Minimum Operational TSR: %.1f\n', min_TSR_operational);
fprintf('  System Inertia: %.4f kg·m²\n', J_total);
fprintf('  Optimal Torque Constant K_opt: %.6f (reference only)\n', K_opt);

%% ========== ENHANCED PMSG ELECTRICAL MODEL ========== %%
P_electrical = zeros(1, NS);
I_phase = zeros(1, NS);
V_phase = zeros(1, NS);
T_em = zeros(1, NS);
T_em_ref_arr = zeros(1, NS);   % Torque comandado (antes do smoothing)
T_em_actual  = zeros(1, NS);   % Torque efetivo (após smoothing, consistente com I_phase)
P_copper = zeros(1, NS);
P_iron = zeros(1, NS);
id_ref = zeros(1, NS);
iq_ref = zeros(1, NS);
control_mode = cell(1, NS);
efficiency_gen = zeros(1, NS);
T_winding_vec = zeros(1, NS);

stress_history = zeros(1, NS);
cum_damage = zeros(1, NS);

cut_in_detected = false;
T_mech_prev = 0;
P_mech_prev = 0;

for idx_v0 = 1:NS
    if idx_v0 == 1
        T_winding = T_amb;
    else
        dt = 1.0;
        P_total_losses = P_copper(idx_v0-1) + P_iron(idx_v0-1);
        T_steady_state = T_amb + P_total_losses * R_th;
        tau_th = R_th * C_th;
        T_winding = T_steady_state - (T_steady_state - T_winding_vec(idx_v0-1)) * exp(-dt/tau_th);
    end
    T_winding_vec(idx_v0) = T_winding;

    Rs_actual = Rs * (1 + alpha_cu * (T_winding - T_amb));

    V0 = global_output.V0(idx_v0);
    T_mech = global_output.Q_generator(idx_v0);
    P_mech_shaft = global_output.P_generator(idx_v0);

    omega_e = P * Omega_generator;

    T_mech_available = T_mech;
    P_mech_available = P_mech_shaft;

    TSR_current = global_output.TSR(idx_v0);

    if idx_v0 > 1
        torque_gradient = (T_mech - T_mech_prev) / (V0 - global_output.V0(idx_v0-1));
        positive_torque_gradient = (torque_gradient > 0);
    else
        positive_torque_gradient = true;
    end

    % ========== CUT-IN BASED ON WIND SPEED (REVISED) ==========
    % ===== REVISED: direct cut-in based on wind speed (no filter lag) =====
    cut_in_detected = (V0 >= V_cut_in_to_mppt);
    if idx_v0 > 1 && cut_in_detected && ~strcmp(last_control_mode, 'CUT-IN')
        fprintf('  CUT-IN ACTIVE: V0=%.2f m/s, P_mech=%.1fW\n', V0, P_mech_available);
    end

    if idx_v0 > 1
        dt = 1.0;
        tau_filter = 2.0;
        alpha_filter = dt / (tau_filter + dt);
        P_mech_filtered = alpha_filter * P_mech_shaft + (1 - alpha_filter) * P_mech_shaft_prev;
    else
        P_mech_filtered = P_mech_shaft;
    end
    P_mech_shaft_prev = P_mech_filtered;

    % ========== STATE MACHINE (WIND SPEED BASED) ==========
    if idx_v0 == 1
        current_mode = 'CUT-IN';
        mode_persistence = 0;
    else
        switch last_control_mode
            case 'CUT-IN'
                if V0 >= V_cut_in_to_mppt + V_hysteresis
                    mode_persistence = min(mode_persistence + 1, 3);
                    if mode_persistence >= 2
                        current_mode = 'MPPT';
                    else
                        current_mode = 'CUT-IN';
                    end
                else
                    mode_persistence = max(mode_persistence - 1, 0);
                    current_mode = 'CUT-IN';
                end

            case 'MPPT'
                if V0 >= V_mppt_to_mpta + V_hysteresis
                    mode_persistence = min(mode_persistence + 1, 3);
                    if mode_persistence >= 2
                        current_mode = 'MPTA';
                    else
                        current_mode = 'MPPT';
                    end
                elseif V0 <= V_cut_in_to_mppt - V_hysteresis
                    mode_persistence = max(mode_persistence - 1, 0);
                    current_mode = 'CUT-IN';
                else
                    current_mode = 'MPPT';
                end

            case 'MPTA'
                if V0 <= V_mppt_to_mpta - V_hysteresis
                    mode_persistence = min(mode_persistence + 1, 3);
                    if mode_persistence >= 2
                        current_mode = 'MPPT';
                    else
                        current_mode = 'MPTA';
                    end
                else
                    mode_persistence = 0;
                    current_mode = 'MPTA';
                end
        end
    end

    T_mech_prev = T_mech;
    P_mech_prev = P_mech_shaft;

    last_control_mode = current_mode;
    control_mode{idx_v0} = current_mode;

    % ========== MODE-SPECIFIC CONTROL ==========
    if strcmp(current_mode, 'MPTA')
        % MPTA: limita potência elétrica ao nominal
        P_mech_max = global_output.P_generator(idx_v0);
        P_elec_target = min(P_nominal, P_mech_max * 0.92);
        P_em_target = P_elec_target / 0.97;
        T_em_ref = P_em_target / Omega_generator;
        T_em_ref = min(T_em_ref, T_mech * 0.98);

        iq_ref(idx_v0) = (2/3) * T_em_ref / (P * lambda_m);
        id_ref(idx_v0) = 0;

        [id_ref(idx_v0), iq_ref(idx_v0)] = apply_field_weakening_analytical(...
            id_ref(idx_v0), iq_ref(idx_v0), Rs_actual, Ld, Lq, lambda_m, omega_e, V_ph_max, I_max_continuous);

        I_peak = sqrt(id_ref(idx_v0)^2 + iq_ref(idx_v0)^2);
        if I_peak > I_max_continuous
            [id_ref(idx_v0), iq_ref(idx_v0)] = calculate_mtpa_currents(...
                I_max_continuous, lambda_m, Ld, Lq, P);
            T_em_ref = (3/2) * P * (lambda_m * iq_ref(idx_v0) + (Ld - Lq) * id_ref(idx_v0) * iq_ref(idx_v0));
        end

        P_em_actual = T_em_ref * Omega_generator;
        P_copper(idx_v0) = (3/2) * Rs_actual * (id_ref(idx_v0)^2 + iq_ref(idx_v0)^2);
        P_iron(idx_v0) = calculate_enhanced_iron_losses(omega_e, lambda_m, id_ref(idx_v0), iq_ref(idx_v0), Ld, Lq, K_h, K_e, T_winding);
        P_rot = K_rot * Omega_generator^2;

        P_electrical(idx_v0) = max(0, P_em_actual - P_copper(idx_v0) - P_iron(idx_v0) - P_rot);

        if P_electrical(idx_v0) > P_mech_max * 0.95
            P_electrical(idx_v0) = P_mech_max * 0.95;
        end

        fprintf('  MPTA: Target=%.0fW, Mechanical=%.0fW, Electrical=%.0fW\n', ...
                P_nominal, P_mech_max, P_electrical(idx_v0));

    elseif strcmp(current_mode, 'MPPT')
        % ========== MPPT PARA ROTAÇÃO FIXA ==========
        % Não há lei K_opt·Ω² útil (Ω constante). Extrai o torque disponível.
        if cut_in_detected
            T_em_ref = T_mech * 0.99;   % 5% margem de estabilidade

            iq_ref(idx_v0) = (2/3) * T_em_ref / (P * lambda_m);
            id_ref(idx_v0) = 0;

            [id_ref(idx_v0), iq_ref(idx_v0)] = apply_field_weakening_analytical(...
                id_ref(idx_v0), iq_ref(idx_v0), Rs_actual, Ld, Lq, lambda_m, omega_e, V_ph_max, I_max_continuous);

            fprintf('  MPPT (fixed-RPM): T_mech=%.3f Nm, T_em=%.3f Nm\n', ...
                    T_mech, T_em_ref);
        else
            iq_ref(idx_v0) = 0;
            id_ref(idx_v0) = 0;
            T_em_ref = 0;
        end

    else  % CUT-IN
        iq_ref(idx_v0) = 0;
        id_ref(idx_v0) = 0;
        T_em_ref = 0;
    end

    % ========== CURRENT SMOOTHING ==========
if idx_v0 == 1 || ~strcmp(control_mode{idx_v0}, control_mode{idx_v0-1})
    % At mode transitions: initialize to target (no lag accumulating)
    id_smoothed = id_ref(idx_v0);
    iq_smoothed = iq_ref(idx_v0);
else
    % Within same mode: gentle smoothing only
    smoothing_factor = 0.8;   % was 0.2 — faster convergence
    id_smoothed = (1 - smoothing_factor) * id_smoothed + smoothing_factor * id_ref(idx_v0);
    iq_smoothed = (1 - smoothing_factor) * iq_smoothed + smoothing_factor * iq_ref(idx_v0);
end
    id_actual = id_smoothed;
    iq_actual = iq_smoothed;
    id_ref(idx_v0) = id_actual;
    iq_ref(idx_v0) = iq_actual;

    vd = Rs_actual * id_ref(idx_v0) - omega_e * Lq * iq_ref(idx_v0);
    vq = Rs_actual * iq_ref(idx_v0) + omega_e * Ld * id_ref(idx_v0) + omega_e * lambda_m;

    if ~strcmp(current_mode, 'MPTA')
        P_copper(idx_v0) = (3/2) * Rs_actual * (id_ref(idx_v0)^2 + iq_ref(idx_v0)^2);
        P_iron(idx_v0) = calculate_enhanced_iron_losses(omega_e, lambda_m, id_ref(idx_v0), iq_ref(idx_v0), Ld, Lq, K_h, K_e, T_winding);
        P_rot = K_rot * Omega_generator^2;
        P_em = T_em_ref * Omega_generator;

        if strcmp(current_mode, 'CUT-IN')
            P_electrical(idx_v0) = 0;
        else
            P_electrical(idx_v0) = max(0, P_em - P_copper(idx_v0) - P_iron(idx_v0) - P_rot);
        end
    end

    % ===== RECALCULATE LOSSES WITH POST-SMOOTHING CURRENTS (fixes derating bug) =====
    P_copper(idx_v0) = (3/2) * Rs_actual * (id_ref(idx_v0)^2 + iq_ref(idx_v0)^2);
    P_iron(idx_v0)   = calculate_enhanced_iron_losses(omega_e, lambda_m, ...
                          id_ref(idx_v0), iq_ref(idx_v0), Ld, Lq, K_h, K_e, T_winding);
    P_total_losses   = P_copper(idx_v0) + P_iron(idx_v0);
    T_winding        = T_amb + P_total_losses * R_th;
    T_winding_vec(idx_v0) = T_winding;

    if T_winding > T_derating_start
        derating_factor = 1 - ((T_winding - T_derating_start) / (T_max_winding - T_derating_start));
        derating_factor = max(0.5, derating_factor);

        if derating_factor < 0.98
            P_elec_target_derated = P_nominal * derating_factor;
            if strcmp(current_mode, 'MPTA')
                T_em_ref = P_elec_target_derated / (Omega_generator * eta_gen_est);
                T_em_ref = min(T_em_ref, T_mech);
                iq_ref(idx_v0) = (2/3) * T_em_ref / (P * lambda_m);
                id_ref(idx_v0) = 0;
                P_copper(idx_v0) = (3/2) * Rs_actual * (id_ref(idx_v0)^2 + iq_ref(idx_v0)^2);
                P_em = T_em_ref * Omega_generator;
                P_electrical(idx_v0) = max(0, P_em - P_copper(idx_v0) - P_iron(idx_v0) - P_rot);
                fprintf('  Thermal derating: %.0f%% at %.1f°C\n', derating_factor*100, T_winding);
            end
        end
    end

    % Torque comandado (referência) e torque efetivo (após smoothing)
    T_em_ref_arr(idx_v0) = T_em_ref;
    T_em_actual(idx_v0)  = (3/2) * P * (lambda_m * iq_ref(idx_v0) + (Ld - Lq) * id_ref(idx_v0) * iq_ref(idx_v0));

    I_phase(idx_v0) = sqrt(id_ref(idx_v0)^2 + iq_ref(idx_v0)^2);
    V_phase(idx_v0) = sqrt(vd^2 + vq^2) / sqrt(2);
    T_em(idx_v0) = T_em_actual(idx_v0);   % <-- CHANGED: agora consistente com I_phase

    if P_em > 0
        efficiency_gen(idx_v0) = P_electrical(idx_v0) / P_em;
    else
        efficiency_gen(idx_v0) = 0;
    end

    % ========== FATIGUE ==========
    if idx_v0 == 1
        stress_history(idx_v0) = 0;
        cum_damage(idx_v0) = 0;
    else
        I_current = I_phase(idx_v0);
        I_prev = I_phase(idx_v0-1);

        mechanical_stress_current = base_mechanical_stress * (I_current / I_nominal)^2;
        mechanical_stress_prev = base_mechanical_stress * (I_prev / I_nominal)^2;

        delta_T_current = T_winding_vec(idx_v0) - T_amb;
        delta_T_prev = T_winding_vec(idx_v0-1) - T_amb;

        thermal_stress_current = base_thermal_stress * delta_T_current;
        thermal_stress_prev = base_thermal_stress * delta_T_prev;

        vibration_stress = 0.08e6 * (Omega_generator / Omega_nominal)^1.5;

        total_stress_current = mechanical_stress_current + thermal_stress_current + vibration_stress;
        total_stress_prev = mechanical_stress_prev + thermal_stress_prev + vibration_stress;

        stress_range = abs(total_stress_current - total_stress_prev);
        mean_stress = (total_stress_current + total_stress_prev) / 2;

        if mean_stress < S_N && stress_range > 0
            stress_corrected = stress_range / (1 - mean_stress/S_N);
        else
            stress_corrected = stress_range;
        end

        stress_corrected = max(stress_corrected, stress_min_detectable);
        stress_history(idx_v0) = stress_corrected / 1e6;

        if stress_corrected > damage_threshold * S_N
            N_fail = (S_N / stress_corrected)^m;
            damage_increment = 0.5 / N_fail;
        else
            damage_increment = 1e-8 * (I_current / I_nominal)^2;
        end

        if ~strcmp(control_mode{idx_v0}, control_mode{idx_v0-1})
            transition_damage = 5e-7;
            damage_increment = damage_increment + transition_damage;
        end

        cum_damage(idx_v0) = cum_damage(idx_v0-1) + damage_increment;

        if cum_damage(idx_v0) > 1.0
            cum_damage(idx_v0) = 1.0;
            fprintf('  FATIGUE FAILURE DETECTED at V0=%.2f m/s\n', V0);
        end
    end

    if P_electrical(idx_v0) > P_mech_shaft * 1.01
        warning('Electrical power exceeds mechanical input at V0=%.1f m/s', V0);
        P_electrical(idx_v0) = P_mech_shaft * 0.99;
    end

    if I_phase(idx_v0) > I_max_continuous * 1.1 && I_phase(idx_v0) > 0.1
        warning('Phase current exceeds continuous limit at V0=%.1f m/s', V0);
    end

    if V_phase(idx_v0) > V_ph_max * 1.1 && V_phase(idx_v0) > 0.1
        warning('Phase voltage exceeds maximum at V0=%.1f m/s', V0);
    end
end

%% ========== DATA EXPORT ========== %%
mid_point = round(NS/2);
cfd_data.V0 = global_output.V0(mid_point);
cfd_data.r = r;
cfd_data.a = global_output.a(mid_point, :);
cfd_data.ap = global_output.ap(mid_point, :);
cfd_data.alpha = global_output.alpha(mid_point, :);
cfd_data.phi = global_output.phi(mid_point, :);
cfd_data.CL = global_output.CL(mid_point, :);
cfd_data.CD = global_output.CD(mid_point, :);
cfd_data.F = global_output.F(mid_point, :);
cfd_data.Re_local = global_output.Re_local(mid_point, :);

global_output.P_electrical = P_electrical;
global_output.I_phase = I_phase;
global_output.V_phase = V_phase;
global_output.T_em = T_em;
global_output.T_em_ref    = T_em_ref_arr;
global_output.T_em_actual = T_em_actual;
global_output.P_copper = P_copper;
global_output.P_iron = P_iron;
global_output.id_ref = id_ref;
global_output.iq_ref = iq_ref;
global_output.control_mode = control_mode;
global_output.efficiency_gen = efficiency_gen;
global_output.T_winding = T_winding_vec;
global_output.stress_range = stress_history;
global_output.damage = cum_damage;
global_output.V_cut_in_to_mppt = V_cut_in_to_mppt;
global_output.V_mppt_to_mpta = V_mppt_to_mpta;

%% ========== VALIDATION ========== %%
validation = validate_simulation(global_output, P_nominal, V0_nominal);
if ~validation.passed
    fprintf('\n=== VALIDATION WARNINGS ===\n');
    for i = 1:length(validation.warnings)
        fprintf('⚠ %s\n', validation.warnings{i});
    end
else
    fprintf('\n=== VALIDATION PASSED ===\n');
end

damage_validation = validate_damage_accumulation(global_output, I_nominal);
if ~damage_validation.passed
    fprintf('\n=== DAMAGE VALIDATION WARNINGS ===\n');
    for i = 1:length(damage_validation.messages)
        fprintf('⚠ %s\n', damage_validation.messages{i});
    end
else
    fprintf('\n=== DAMAGE VALIDATION PASSED ===\n');
    fprintf('✓ %s\n', damage_validation.messages{1});
end

cfd_validation = validate_cfd_readiness(global_output);
if ~cfd_validation.passed
    fprintf('\n=== CFD VALIDATION WARNINGS ===\n');
    for i = 1:length(cfd_validation.messages)
        fprintf('⚠ %s\n', cfd_validation.messages{i});
    end
else
    fprintf('\n=== CFD VALIDATION PASSED ===\n');
end

%% ========== EXPORT RESULTS ========== %%
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
resultsFolder = fullfile('Simulation_Results', ['Run_' timestamp]);
if ~isfolder(resultsFolder)
    mkdir(resultsFolder);
end

save(fullfile(resultsFolder, 'Complete_Performance.mat'), 'global_output', '-v7.3');

if isfield(global_output, 'a') && ~isempty(global_output.a) && size(global_output.a, 2) > 0
    mid_point = round(size(global_output.a, 1)/2);
    if isfield(global_output, 'radial_positions') && ~isempty(global_output.radial_positions)
        r_positions = global_output.radial_positions;
    else
        num_radial_points = size(global_output.a, 2);
        r_positions = linspace(min_r_ratio * R, R, num_radial_points);
    end
    cfd_data.V0 = global_output.V0(mid_point);
    cfd_data.r = r_positions;
    cfd_data.a = global_output.a(mid_point, :);
    cfd_data.ap = global_output.ap(mid_point, :);
    cfd_data.alpha = global_output.alpha(mid_point, :);
    cfd_data.phi = global_output.phi(mid_point, :);
    cfd_data.CL = global_output.CL(mid_point, :);
    cfd_data.CD = global_output.CD(mid_point, :);
    cfd_data.F = global_output.F(mid_point, :);
    cfd_data.Re_local = global_output.Re_local(mid_point, :);
    save(fullfile(resultsFolder, 'CFD_Validation_Data.mat'), 'cfd_data');
    fprintf('CFD validation data exported for V0 = %.2f m/s\n', cfd_data.V0);
end

excelFile = fullfile(resultsFolder, 'Complete_Performance_Analysis.xlsx');
headers = {
    'V0 (m/s)', 'TSR', 'Cp', 'Available_Power (W)', ...
    'Rotor_Power (W)', 'Generator_Power_mechanical (W)', 'Rotor_Torque (Nm)', ...
    'Generator_Torque (Nm)', 'Electrical_Power (W)', 'I_phase (A)', ...
    'V_phase (V)', 'T_em_actual (Nm)', 'T_em_ref (Nm)', 'Copper_Loss (W)', ...
    'Iron_Loss (W)', 'Efficiency_generator', 'Control_Mode', ...
    'T_winding (C)', 'Stress_Range (MPa)', 'Cumulative_Damage'
};
NS = length(global_output.V0);
dataCell = cell(NS+1, length(headers));
dataCell(1,:) = headers;
for i = 1:NS
    dataCell{i+1,1} = global_output.V0(i);
    dataCell{i+1,2} = global_output.TSR(i);
    dataCell{i+1,3} = global_output.Cp(i);
    dataCell{i+1,4} = global_output.P_available(i);
    dataCell{i+1,5} = global_output.P_rotor(i);
    dataCell{i+1,6} = global_output.P_mechanical(i);
    dataCell{i+1,7} = global_output.Q_rotor(i);
    dataCell{i+1,8} = global_output.Q_generator(i);
    dataCell{i+1,9} = global_output.P_electrical(i);
    dataCell{i+1,10} = global_output.I_phase(i);
    dataCell{i+1,11} = global_output.V_phase(i);
    dataCell{i+1,12} = global_output.T_em_actual(i);
    dataCell{i+1,13} = global_output.T_em_ref(i);
    dataCell{i+1,14} = global_output.P_copper(i);
    dataCell{i+1,15} = global_output.P_iron(i);
    dataCell{i+1,16} = global_output.efficiency_gen(i);
    dataCell{i+1,17} = global_output.control_mode{i};
    dataCell{i+1,18} = global_output.T_winding(i);
    dataCell{i+1,19} = global_output.stress_range(i);
    dataCell{i+1,20} = global_output.damage(i);
end
writecell(dataCell, excelFile, 'Sheet', 'Results');

fprintf('Results exported to: %s\n', resultsFolder);

%% ========== PLOT SETTINGS ========== %%
plotFolder = fullfile(resultsFolder, 'Figures');
if ~isfolder(plotFolder)
    mkdir(plotFolder);
end

try
    set(groot, 'defaultAxesFontName', 'Times New Roman');
    set(groot, 'defaultTextFontName', 'Times New Roman');
    set(groot, 'defaultUicontrolFontName', 'Times New Roman');
    set(groot, 'defaultUitableFontName', 'Times New Roman');
    set(groot, 'defaultUipanelFontName', 'Times New Roman');
catch
    warning('Times New Roman não disponível. Usando fonte padrão.');
end

set(groot, 'defaultAxesFontSize', 22);
set(groot, 'defaultTextFontSize', 22);
set(groot, 'defaultLineLineWidth', 3.0);
set(groot, 'defaultAxesLineWidth', 2.0);
set(groot, 'defaultAxesGridLineStyle', '-');
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot, 'defaultLegendInterpreter', 'tex');
set(groot, 'defaultTextInterpreter', 'tex');

colors = struct();
colors.blue = [0, 0.4470, 0.7410];
colors.orange = [0.8500, 0.3250, 0.0980];
colors.green = [0.4660, 0.6740, 0.1880];
colors.red = [0.6350, 0.0780, 0.1840];
colors.purple = [0.4940, 0.1840, 0.5560];
colors.yellow = [0.9290, 0.6940, 0.1250];
colors.cyan = [0.3010, 0.7450, 0.9330];
colors.gray = [0.5, 0.5, 0.5];

% Rótulos das linhas de referência (padrão já estabelecido no código)
labelVnom = sprintf('V_{nom} = %.2f m/s', V0_nominal);

fprintf('Gerando todos os gráficos...\n');

%% Figure 1: TSR vs Wind Speed
figure(1);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.V0, global_output.TSR, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'TSR');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('TSR', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'TSR_vs_WindSpeed.svg'));
saveas(gcf, fullfile(plotFolder, 'TSR_vs_WindSpeed.fig'));
close(gcf);

%% Figure 2: TSR vs Cp
TSR_nominal = Omega_rotor * R / V0_nominal;
labelTSRopt = sprintf('TSR_{nom} = %.2f', TSR_nominal);
figure(2);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.TSR, global_output.Cp, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
xline(TSR_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelTSRopt);
xlabel('TSR (\lambda)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('C_p', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([0 15]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'TSR_vs_Cp.svg'));
saveas(gcf, fullfile(plotFolder, 'TSR_vs_Cp.fig'));
close(gcf);

%% Figure 3: Cp vs Wind Speed
figure(3);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.V0, global_output.Cp, '-^', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'C_p');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
yline(0.593, '--k', 'LineWidth', 3.0, 'DisplayName', 'Limite de Betz');
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('C_p', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Cp_vs_WindSpeed.svg'));
saveas(gcf, fullfile(plotFolder, 'Cp_vs_WindSpeed.fig'));
close(gcf);

%% Figure 4: Cq vs Wind Speed
figure(4);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.V0, global_output.Cq, '-d', 'Color', colors.cyan, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.cyan, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'C_q');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('C_q', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Cq_vs_WindSpeed.svg'));
saveas(gcf, fullfile(plotFolder, 'Cq_vs_WindSpeed.fig'));
close(gcf);

%% Figure 5: Available Power
figure(5);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.V0, global_output.P_available, '-v', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Potência Disponível');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Potência Disponível [W]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Available_Power.svg'));
saveas(gcf, fullfile(plotFolder, 'Available_Power.fig'));
close(gcf);

%% Figure 6: Mechanical Power
figure(6);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.V0, global_output.P_rotor, '-o', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Potência Mecânica');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Potência Mecânica [W]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Mechanical_Power.svg'));
saveas(gcf, fullfile(plotFolder, 'Mechanical_Power.fig'));
close(gcf);

%% Figure 7: Mechanical Torque
figure(7);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.V0, global_output.Q_rotor, '-s', 'Color', colors.purple, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Torque Mecânico');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Torque Mecânico [Nm]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Mechanical_Torque.svg'));
saveas(gcf, fullfile(plotFolder, 'Mechanical_Torque.fig'));
close(gcf);

%% Figure 8: Power Comparison
figure(8);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
hold on;
plot(global_output.V0, global_output.P_available, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Disponível');
plot(global_output.V0, global_output.P_rotor, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Mecânica (com perdas)');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Potência [W]', 'FontSize', 24, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Power_Comparison_Aerodynamic.svg'));
saveas(gcf, fullfile(plotFolder, 'Power_Comparison_Aerodynamic.fig'));
close(gcf);

%% Figure 9: Electrical Power vs Wind Speed
figure(9);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
hold on;
plot(global_output.V0, global_output.P_electrical, '-^', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Potência Elétrica');
yline(P_nominal, '--r', 'LineWidth', 3.0, 'DisplayName', sprintf('Potência Nominal (%.0f W)', P_nominal));
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Potência Elétrica [W]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Electrical_Power.svg'));
saveas(gcf, fullfile(plotFolder, 'Electrical_Power.fig'));
close(gcf);

%% Figure 10: Power Comparison (full chain)
figure(10);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
hold on;
plot(global_output.V0, global_output.P_available, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Disponível');
plot(global_output.V0, global_output.P_rotor, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Rotor');
plot(global_output.V0, global_output.P_mechanical, '-^', 'Color', colors.purple, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Gerador (mec.)');
plot(global_output.V0, global_output.P_electrical, '-d', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Elétrica');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Potência [W]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Power_Comparison_Electrical.svg'));
saveas(gcf, fullfile(plotFolder, 'Power_Comparison_Electrical.fig'));
close(gcf);

%% Figure 11: Torque Comparison
figure(11);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
hold on;
plot(global_output.V0, global_output.Q_rotor, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Rotor');
plot(global_output.V0, global_output.Q_generator, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Gerador (após caixa)');
plot(global_output.V0, global_output.T_em, '-^', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Eletromagnético');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Torque [Nm]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Torque_Comparison.svg'));
saveas(gcf, fullfile(plotFolder, 'Torque_Comparison.fig'));
close(gcf);

%% Figure 12: System Efficiency
figure(12);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
hold on;
overall_eff = (global_output.P_electrical ./ global_output.P_available) * 100;
plot(global_output.V0, overall_eff, '-d', 'Color', colors.purple, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Global');
plot(global_output.V0, global_output.Cp*100, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Aerodinâmico (C_p)');
plot(global_output.V0, global_output.efficiency_gen*100, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Gerador');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Eficiência [%]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
ylim([0 100]);
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'System_Efficiency.svg'));
saveas(gcf, fullfile(plotFolder, 'System_Efficiency.fig'));
close(gcf);

%% Figure 13: Electrical Parameters (I vs V)
figure(13);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
yyaxis left;
plot(global_output.V0, global_output.I_phase, '-o', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Corrente');
ylabel('Corrente de Fase [A]', 'FontSize', 24, 'FontWeight', 'bold');
yyaxis right;
plot(global_output.V0, global_output.V_phase, '-s', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Tensão');
ylabel('Tensão de Fase [V]', 'FontSize', 24, 'FontWeight', 'bold');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
ax = gca;
ax.YAxis(1).Color = colors.red;
ax.YAxis(1).Label.Color = colors.red;
ax.YAxis(2).Color = colors.blue;
ax.YAxis(2).Label.Color = colors.blue;
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Electrical_Parameters.svg'));
saveas(gcf, fullfile(plotFolder, 'Electrical_Parameters.fig'));
close(gcf);

%% Figure 14: Control Strategy
figure(14);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
mppt_idx = strcmp(global_output.control_mode, 'MPPT');
mpta_idx = strcmp(global_output.control_mode, 'MPTA');
cutin_idx = strcmp(global_output.control_mode, 'CUT-IN');

plot(global_output.V0(mppt_idx), global_output.T_em(mppt_idx), 'o', ...
    'Color', colors.blue, 'LineWidth', 3.0, 'MarkerSize', 10, ...
    'MarkerFaceColor', colors.blue, 'MarkerEdgeColor', 'k', 'DisplayName', 'MPPT');
hold on;
plot(global_output.V0(mpta_idx), global_output.T_em(mpta_idx), 's', ...
    'Color', colors.red, 'LineWidth', 3.0, 'MarkerSize', 10, ...
    'MarkerFaceColor', colors.red, 'MarkerEdgeColor', 'k', 'DisplayName', 'MPTA');
plot(global_output.V0(cutin_idx), global_output.T_em(cutin_idx), '^', ...
    'Color', colors.gray, 'LineWidth', 3.0, 'MarkerSize', 10, ...
    'MarkerFaceColor', colors.gray, 'MarkerEdgeColor', 'k', 'DisplayName', 'Cut-in');
plot(global_output.V0, global_output.Q_generator, '-', 'Color', colors.green, ...
    'LineWidth', 3.0, 'DisplayName', 'Gerador');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Torque [Nm]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Control_Strategy.svg'));
saveas(gcf, fullfile(plotFolder, 'Control_Strategy.fig'));
close(gcf);

%% Figure 15: Winding Temperature
figure(15);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
plot(global_output.V0, global_output.T_winding, '-o', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Temp. Enrolamento');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
yline(T_derating_start, '--y', 'LineWidth', 3.0, 'DisplayName', 'Início Derating');
yline(T_max_winding, '--r', 'LineWidth', 3.0, 'DisplayName', 'Temp. Máx.');
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Temperatura do Enrolamento [°C]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Generator_Temperature.svg'));
saveas(gcf, fullfile(plotFolder, 'Generator_Temperature.fig'));
close(gcf);

%% Figure 16: Fatigue
figure(16);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');
yyaxis left;
plot(global_output.V0, global_output.stress_range, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Faixa de Tensão');
ylabel('Faixa de Tensão [MPa]', 'FontSize', 24, 'FontWeight', 'bold');
yyaxis right;
plot(global_output.V0, global_output.damage*100, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 10, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Dano');
ylabel('Dano Acumulado [%]', 'FontSize', 24, 'FontWeight', 'bold');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 18);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 18);
box on;
saveas(gcf, fullfile(plotFolder, 'Fatigue_analysis.svg'));
saveas(gcf, fullfile(plotFolder, 'Fatigue_analysis.fig'));
close(gcf);

%% Figure 17: Power Curves Detailed (2x2, sem TSR vs Vento)
figure(17);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
plot(global_output.V0, global_output.P_available, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Disponível');
hold on;
plot(global_output.V0, global_output.P_rotor, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Rotor');
plot(global_output.V0, global_output.P_mechanical, '-^', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Mecânica');
plot(global_output.V0, global_output.P_electrical, '-d', 'Color', colors.purple, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Elétrica');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Potência [W]', 'FontSize', 18, 'FontWeight', 'bold');
title('Curvas de Potência', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
plot(global_output.V0, global_output.Cp, '-o', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('C_p', 'FontSize', 18, 'FontWeight', 'bold');
title('C_p vs Vento', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
plot(global_output.TSR, global_output.Cp, '-^', 'Color', colors.cyan, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.cyan, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(TSR_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelTSRopt);
xlabel('TSR', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('C_p', 'FontSize', 18, 'FontWeight', 'bold');
title('C_p vs TSR', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([0 15]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Power_Curves_Electrical_Detailed.svg'));
saveas(gcf, fullfile(plotFolder, 'Power_Curves_Electrical_Detailed.fig'));
close(gcf);

%% Figure 18: Efficiency Curves Detailed (2x2)
figure(18);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
plot(global_output.V0, global_output.efficiency_gen * 100, '-o', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Eficiência Gerador [%]', 'FontSize', 18, 'FontWeight', 'bold');
title('Eficiência do Gerador', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
ylim([0 100]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
system_efficiency = global_output.P_electrical ./ global_output.P_available * 100;
system_efficiency(isinf(system_efficiency) | isnan(system_efficiency)) = 0;
plot(global_output.V0, system_efficiency, '-s', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Eficiência Global [%]', 'FontSize', 18, 'FontWeight', 'bold');
title('Eficiência Global', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
ylim([0 50]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
drivetrain_efficiency = global_output.P_mechanical ./ global_output.P_rotor * 100;
drivetrain_efficiency(isinf(drivetrain_efficiency) | isnan(drivetrain_efficiency)) = 0;
plot(global_output.V0, drivetrain_efficiency, '-^', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Eficiência Trem de Força [%]', 'FontSize', 18, 'FontWeight', 'bold');
title('Trem de Força', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
ylim([95 100]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
plot(global_output.P_mechanical, global_output.efficiency_gen * 100, 'o', ...
    'MarkerSize', 8, 'Color', colors.purple, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k');
xlabel('Potência Mecânica [W]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Eficiência Gerador [%]', 'FontSize', 18, 'FontWeight', 'bold');
title('Eficiência vs Potência', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Efficiency_Curves_Detailed.svg'));
saveas(gcf, fullfile(plotFolder, 'Efficiency_Curves_Detailed.fig'));
close(gcf);

%% Figure 19: Electrical Parameters Detailed (2x2)
figure(19);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
plot(global_output.V0, global_output.I_phase, '-o', 'Color', colors.orange, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.orange, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Corrente de Fase');
hold on;
yline(I_max_continuous, '--r', 'LineWidth', 3.0, 'DisplayName', 'I_{max,cont}');
yline(I_max_peak, '--k', 'LineWidth', 3.0, 'DisplayName', 'I_{max,pico}');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Corrente de Fase [A]', 'FontSize', 18, 'FontWeight', 'bold');
title('Corrente de Fase', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
plot(global_output.V0, global_output.V_phase, '-s', 'Color', colors.cyan, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.cyan, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Tensão de Fase');
hold on;
yline(V_ph_max, '--r', 'LineWidth', 3.0, 'DisplayName', 'V_{max}');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Tensão de Fase [V]', 'FontSize', 18, 'FontWeight', 'bold');
title('Tensão de Fase', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
plot(global_output.V0, global_output.T_em, '-^', 'Color', colors.purple, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Torque Eletromagnético [Nm]', 'FontSize', 18, 'FontWeight', 'bold');
title('Torque Eletromagnético', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
plot(global_output.V0, global_output.id_ref, '-d', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'i_d');
hold on;
plot(global_output.V0, global_output.iq_ref, '-v', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'i_q');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Correntes de Referência [A]', 'FontSize', 18, 'FontWeight', 'bold');
title('Correntes d-q', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Electrical_Parameters_Detailed.svg'));
saveas(gcf, fullfile(plotFolder, 'Electrical_Parameters_Detailed.fig'));
close(gcf);

%% Figure 20: Control Strategy Analysis (2x2)
figure(20);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

colors_mode = zeros(length(global_output.control_mode), 3);
for i = 1:length(global_output.control_mode)
    switch global_output.control_mode{i}
        case 'MPPT'
            colors_mode(i,:) = colors.green;
        case 'MPTA'
            colors_mode(i,:) = colors.orange;
        case 'CUT-IN'
            colors_mode(i,:) = colors.gray;
        otherwise
            colors_mode(i,:) = [0.5, 0.5, 0.5];
    end
end

subplot(2,2,1);
scatter(global_output.V0, global_output.P_mechanical, 50, colors_mode, 'filled');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Potência Mecânica [W]', 'FontSize', 18, 'FontWeight', 'bold');
title('Estratégia de Controle', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
scatter(global_output.P_mechanical, global_output.efficiency_gen * 100, 50, colors_mode, 'filled');
xlabel('Potência Mecânica [W]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Eficiência Gerador [%]', 'FontSize', 18, 'FontWeight', 'bold');
title('Controle vs Eficiência', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
mode_numeric = zeros(size(global_output.control_mode));
for i = 1:length(global_output.control_mode)
    switch global_output.control_mode{i}
        case 'CUT-IN', mode_numeric(i) = 0;
        case 'MPPT', mode_numeric(i) = 1;
        case 'MPTA', mode_numeric(i) = 2;
        otherwise, mode_numeric(i) = 0;
    end
end
plot(global_output.V0, mode_numeric, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Modo de Controle', 'FontSize', 18, 'FontWeight', 'bold');
yticks([0 1 2]);
yticklabels({'CUT-IN', 'MPPT', 'MPTA'});
title('Transições de Modo', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
histogram(categorical(global_output.control_mode), 'FaceColor', colors.blue);
xlabel('Modo de Controle', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequência', 'FontSize', 18, 'FontWeight', 'bold');
title('Distribuição dos Modos', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Control_Strategy_Analysis.svg'));
saveas(gcf, fullfile(plotFolder, 'Control_Strategy_Analysis.fig'));
close(gcf);

%% Figure 21: Thermal Analysis (2x2)
figure(21);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
plot(global_output.V0, global_output.T_winding, '-o', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Temp. Enrolamento');
hold on;
yline(T_max_winding, '--r', 'LineWidth', 3.0, 'DisplayName', 'Temp. Máx.');
yline(T_derating_start, '--y', 'LineWidth', 3.0, 'DisplayName', 'Início Derating');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Temperatura [°C]', 'FontSize', 18, 'FontWeight', 'bold');
title('Temperatura do Enrolamento', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
plot(global_output.V0, global_output.P_copper, '-s', 'Color', colors.orange, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.orange, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Perdas Enrolamento');
hold on;
plot(global_output.V0, global_output.P_iron, '-^', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Perdas Núcleo');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Perdas [W]', 'FontSize', 18, 'FontWeight', 'bold');
title('Perdas no Gerador', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
plot(global_output.I_phase, global_output.T_winding, 'o', 'MarkerSize', 8, ...
    'Color', colors.purple, 'MarkerFaceColor', colors.purple, 'MarkerEdgeColor', 'k');
xlabel('Corrente [A]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Temperatura [°C]', 'FontSize', 18, 'FontWeight', 'bold');
title('Temperatura vs Corrente', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
total_losses = global_output.P_copper + global_output.P_iron;
plot(global_output.V0, total_losses, '-d', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Perdas Totais [W]', 'FontSize', 18, 'FontWeight', 'bold');
title('Perdas Totais', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Thermal_Analysis_Detailed.svg'));
saveas(gcf, fullfile(plotFolder, 'Thermal_Analysis_Detailed.fig'));
close(gcf);

%% Figure 22: Fatigue Analysis (2x2)
figure(22);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
semilogy(global_output.V0, max(global_output.damage, 1e-12), '-o', 'Color', colors.purple, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Dano');
hold on;
yline(damage_threshold, '--r', 'LineWidth', 3.0, 'DisplayName', 'Limite');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Dano Acumulado', 'FontSize', 18, 'FontWeight', 'bold');
title('Dano de Fadiga', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
plot(global_output.V0, global_output.stress_range, '-s', 'Color', colors.orange, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.orange, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Faixa de Tensão [MPa]', 'FontSize', 18, 'FontWeight', 'bold');
title('Faixa de Tensão', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
plot(global_output.I_phase, global_output.stress_range, 'o', 'MarkerSize', 8, ...
    'Color', colors.cyan, 'MarkerFaceColor', colors.cyan, 'MarkerEdgeColor', 'k');
xlabel('Corrente [A]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Tensão [MPa]', 'FontSize', 18, 'FontWeight', 'bold');
title('Tensão vs Corrente', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
plot(global_output.T_winding, global_output.stress_range, 'o', 'MarkerSize', 8, ...
    'Color', colors.red, 'MarkerFaceColor', colors.red, 'MarkerEdgeColor', 'k');
xlabel('Temperatura [°C]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Tensão [MPa]', 'FontSize', 18, 'FontWeight', 'bold');
title('Tensão vs Temperatura', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Fatigue_Analysis_Detailed.svg'));
saveas(gcf, fullfile(plotFolder, 'Fatigue_Analysis_Detailed.fig'));
close(gcf);

%% Figure 23: Torque Curves Detailed (2x2)
figure(23);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
plot(global_output.V0, global_output.Q_rotor, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Rotor');
hold on;
plot(global_output.V0, global_output.Q_generator, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Gerador');
plot(global_output.V0, global_output.T_em, '-^', 'Color', colors.green, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.green, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'EM');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Torque [Nm]', 'FontSize', 18, 'FontWeight', 'bold');
title('Comparação de Torques', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
plot(RPM_generator * ones(size(global_output.V0)), global_output.T_em, '-d', ...
    'Color', colors.purple, 'LineWidth', 3.0, 'MarkerSize', 8, ...
    'MarkerFaceColor', colors.purple, 'MarkerEdgeColor', 'k');
xlabel('RPM Gerador', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Torque EM [Nm]', 'FontSize', 18, 'FontWeight', 'bold');
title('Torque vs Velocidade (RPM fixo)', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
plot(global_output.V0, global_output.Ct, '-v', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('C_t', 'FontSize', 18, 'FontWeight', 'bold');
title('Coeficiente de Empuxo', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
plot(global_output.V0, global_output.Cq, '-^', 'Color', colors.cyan, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.cyan, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Velocidade do Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('C_q', 'FontSize', 18, 'FontWeight', 'bold');
title('Coeficiente de Torque', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Torque_Curves_Detailed.svg'));
saveas(gcf, fullfile(plotFolder, 'Torque_Curves_Detailed.fig'));
close(gcf);

%% Figure 24: Current-Voltage Characteristics
figure(24);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
plot(global_output.V0, global_output.I_phase, '-o', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Corrente [A]', 'FontSize', 18, 'FontWeight', 'bold');
title('Corrente vs Vento', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
plot(global_output.V0, global_output.V_phase, '-s', 'Color', colors.red, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.red, ...
    'MarkerEdgeColor', 'k');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Tensão [V]', 'FontSize', 18, 'FontWeight', 'bold');
title('Tensão vs Vento', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
plot(global_output.I_phase, global_output.V_phase, 'o', 'MarkerSize', 8, ...
    'Color', colors.green, 'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k');
xlabel('Corrente [A]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Tensão [V]', 'FontSize', 18, 'FontWeight', 'bold');
title('V vs I', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
plot(global_output.P_mechanical, global_output.I_phase, 'o', 'MarkerSize', 8, ...
    'Color', colors.purple, 'MarkerFaceColor', colors.purple, 'MarkerEdgeColor', 'k');
xlabel('Potência Mecânica [W]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Corrente [A]', 'FontSize', 18, 'FontWeight', 'bold');
title('Corrente vs Potência', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Current_Voltage_Characteristics.svg'));
saveas(gcf, fullfile(plotFolder, 'Current_Voltage_Characteristics.fig'));
close(gcf);

%% Figure 25: Loss Breakdown
figure(25);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

V0_v = global_output.V0(:)';
P_copper_v = global_output.P_copper(:)';
P_iron_v = global_output.P_iron(:)';
P_mechanical_v = global_output.P_mechanical(:)';
P_electrical_v = global_output.P_electrical(:)';

total_losses = P_copper_v + P_iron_v;

subplot(2,2,1);
Y_losses = [P_copper_v; P_iron_v]';
h_area = area(V0_v, Y_losses, 'EdgeColor', 'none');
h_area(1).FaceColor = colors.orange;
h_area(2).FaceColor = colors.blue;
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
legend('Enrolamento', 'Núcleo', labelVnom, 'Location', 'best', 'FontSize', 14);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Perdas [W]', 'FontSize', 18, 'FontWeight', 'bold');
title('Decomposição das Perdas', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
P_mech_safe = max(P_mechanical_v, 1e-6);
useful_percent = min((P_electrical_v ./ P_mech_safe) * 100, 100);
losses_percent = min((total_losses ./ P_mech_safe) * 100, 100);
Y_power_dist = [useful_percent; losses_percent]';
h_area2 = area(V0_v, Y_power_dist, 'EdgeColor', 'none');
h_area2(1).FaceColor = colors.green;
h_area2(2).FaceColor = colors.red;
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
legend('Útil', 'Perdas', labelVnom, 'Location', 'best', 'FontSize', 14);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Distribuição [%]', 'FontSize', 18, 'FontWeight', 'bold');
title('Distribuição de Potência', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
xlim([1 13]);
ylim([0 100]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
mean_copper = mean(P_copper_v);
mean_iron = mean(P_iron_v);
if (mean_copper + mean_iron) > 0
    h_pie = pie([mean_copper, mean_iron], {'Enrolamento', 'Núcleo'});
    if numel(h_pie) >= 4
        set(h_pie(1), 'FaceColor', colors.orange);
        set(h_pie(3), 'FaceColor', colors.blue);
        set(h_pie(2), 'FontSize', 16, 'FontWeight', 'bold');
        set(h_pie(4), 'FontSize', 16, 'FontWeight', 'bold');
    end
else
    text(0.5, 0.5, 'Sem perdas', 'HorizontalAlignment', 'center', ...
        'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
end
title('Perdas Médias', 'FontSize', 20, 'FontWeight', 'bold');

subplot(2,2,4);
total_safe = max(total_losses, 1e-6);
copper_pct = (P_copper_v ./ total_safe) * 100;
iron_pct = (P_iron_v ./ total_safe) * 100;
plot(V0_v, copper_pct, '-o', 'Color', colors.orange, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.orange, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Enrolamento');
hold on;
plot(V0_v, iron_pct, '-s', 'Color', colors.blue, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Núcleo');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('%', 'FontSize', 18, 'FontWeight', 'bold');
title('Distribuição vs Vento', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
ylim([0 100]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Loss_Breakdown_Analysis.svg'));
saveas(gcf, fullfile(plotFolder, 'Loss_Breakdown_Analysis.fig'));
close(gcf);

%% Figure 26: Convergence Analysis
figure(26);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
plot(global_output.V0, global_output.iterations, '-o', 'Color', colors.purple, ...
    'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', colors.purple, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Iterações');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Iterações', 'FontSize', 18, 'FontWeight', 'bold');
title('Convergência BEMT', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
converged_colors = zeros(length(global_output.converged), 3);
converged_colors(global_output.converged, :) = repmat(colors.green, sum(global_output.converged), 1);
converged_colors(~global_output.converged, :) = repmat(colors.red, sum(~global_output.converged), 1);
scatter(global_output.V0, double(global_output.converged), 60, converged_colors, 'filled');
hold on;
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Status', 'FontSize', 18, 'FontWeight', 'bold');
yticks([0 1]);
yticklabels({'Não', 'Sim'});
title('Status de Convergência', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
if isfield(global_output, 'final_errors') && ~isempty(global_output.final_errors)
    semilogy(global_output.V0, max(global_output.final_errors(:,1), 1e-12), '-o', ...
        'Color', colors.blue, 'LineWidth', 3.0, 'MarkerSize', 8, ...
        'MarkerFaceColor', colors.blue, 'MarkerEdgeColor', 'k', 'DisplayName', 'erro a');
    hold on;
    semilogy(global_output.V0, max(global_output.final_errors(:,2), 1e-12), '-s', ...
        'Color', colors.red, 'LineWidth', 3.0, 'MarkerSize', 8, ...
        'MarkerFaceColor', colors.red, 'MarkerEdgeColor', 'k', 'DisplayName', 'erro a''');
    xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
    xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Erro Final', 'FontSize', 18, 'FontWeight', 'bold');
    title('Erros de Convergência', 'FontSize', 20, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 14);
    grid on;
    xlim([1 13]);
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);
end

subplot(2,2,4);
convergence_rate = mean(global_output.converged) * 100;
avg_iterations = mean(global_output.iterations(global_output.converged));
max_iterations = max(global_output.iterations);
stats_text = sprintf(['Taxa Convergência: %.1f%%\n' ...
                     'Média Iterações: %.1f\n' ...
                     'Máx Iterações: %d\n' ...
                     'Não Convergidos: %d'], ...
                     convergence_rate, avg_iterations, max_iterations, ...
                     sum(~global_output.converged));
text(0.1, 0.7, stats_text, 'FontSize', 16, 'FontWeight', 'bold', 'Units', 'normalized');
axis off;
title('Estatísticas', 'FontSize', 20, 'FontWeight', 'bold');

saveas(gcf, fullfile(plotFolder, 'Convergence_Analysis_Detailed.svg'));
saveas(gcf, fullfile(plotFolder, 'Convergence_Analysis_Detailed.fig'));
close(gcf);

%% Figure 27: Operating Analysis
figure(27);
set(gcf, 'Position', [100 100 1800 1050], 'Color', 'white', 'Renderer', 'opengl');

subplot(2,2,1);
low_wind = global_output.V0 <= 4;
medium_wind = (global_output.V0 > 4) & (global_output.V0 <= 11);
high_wind = global_output.V0 > 11;
scatter(global_output.V0(low_wind), global_output.P_electrical(low_wind), 60, ...
    colors.gray, 'filled', 'DisplayName', 'Cut-in (≤4 m/s)');
hold on;
scatter(global_output.V0(medium_wind), global_output.P_electrical(medium_wind), 60, ...
    colors.green, 'filled', 'DisplayName', 'MPPT (4-11 m/s)');
scatter(global_output.V0(high_wind), global_output.P_electrical(high_wind), 60, ...
    colors.red, 'filled', 'DisplayName', 'MPTA (>11 m/s)');
xline(V0_nominal, '--k', 'LineWidth', 3.0, 'DisplayName', labelVnom);
xlabel('Vento [m/s]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Potência Elétrica [W]', 'FontSize', 18, 'FontWeight', 'bold');
title('Pontos de Operação', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
grid on;
xlim([1 13]);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,2);
histogram(global_output.Cp, 20, 'FaceColor', colors.cyan, 'FaceAlpha', 0.7);
xlabel('C_p', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequência', 'FontSize', 18, 'FontWeight', 'bold');
title('Distribuição de C_p', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,3);
scatter(global_output.P_electrical, global_output.efficiency_gen * 100, 50, ...
    global_output.V0, 'filled');
xlabel('Potência Elétrica [W]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Eficiência Gerador [%]', 'FontSize', 18, 'FontWeight', 'bold');
title('Eficiência vs Potência', 'FontSize', 20, 'FontWeight', 'bold');
colorbar;
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

subplot(2,2,4);
mppt_currents = global_output.I_phase(strcmp(global_output.control_mode, 'MPPT'));
mpta_currents = global_output.I_phase(strcmp(global_output.control_mode, 'MPTA'));

if ~isempty(mppt_currents) && ~isempty(mpta_currents)
    if isrow(mppt_currents), mppt_currents = mppt_currents'; end
    if isrow(mpta_currents), mpta_currents = mpta_currents'; end
    all_data = {mppt_currents, mpta_currents};
    group_labels = {'MPPT', 'MPTA'};
    groups = [];
    data_combined = [];
    for i = 1:length(all_data)
        if ~isempty(all_data{i})
            groups = [groups; i * ones(length(all_data{i}), 1)];
            data_combined = [data_combined; all_data{i}];
        end
    end
    if ~isempty(data_combined)
        boxplot(data_combined, groups, 'Labels', group_labels, 'Colors', [colors.blue; colors.red]);
    else
        text(0.5, 0.5, 'Sem dados', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
else
    if isempty(mppt_currents) && isempty(mpta_currents)
        text(0.5, 0.5, 'Sem dados', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    elseif isempty(mppt_currents)
        boxplot(mpta_currents, 'Labels', {'MPTA'}, 'Colors', colors.red);
    else
        boxplot(mppt_currents, 'Labels', {'MPPT'}, 'Colors', colors.blue);
    end
end
ylabel('Corrente [A]', 'FontSize', 18, 'FontWeight', 'bold');
title('Corrente por Modo', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'FontSize', 16);

saveas(gcf, fullfile(plotFolder, 'Operating_Points_Analysis.svg'));
saveas(gcf, fullfile(plotFolder, 'Operating_Points_Analysis.fig'));
close(gcf);

fprintf('\n\n======= SIMULAÇÃO CONCLUÍDA COM SUCESSO =======\n');
fprintf('Gráficos salvos em: %s\n', plotFolder);
fprintf('Resultados salvos em: %s\n', resultsFolder);
fprintf('Tempo total: %.2f s\n', toc(t_start_global));

%% ========== HELPER FUNCTIONS ========== %%

function validation_results = validate_simulation(global_output, P_nominal, V0_nominal)
    validation_results.passed = true;
    validation_results.warnings = {};
    if any(global_output.stress_range > 50)
        validation_results.warnings{end+1} = 'Stress range exceeds 50 MPa';
        validation_results.passed = false;
    end
    if global_output.damage(end) > 0.1
        validation_results.warnings{end+1} = 'Damage accumulation too rapid';
        validation_results.passed = false;
    end
    active_indices = global_output.V0 >= 4;
    active_stress = global_output.stress_range(active_indices);
    if length(active_stress) > 3 && std(active_stress) < 0.1
        validation_results.warnings{end+1} = 'Stress range insufficient variability';
    end
    if global_output.damage(end) == 0 && max(global_output.I_phase) > 2.0
        validation_results.warnings{end+1} = 'No damage accumulation';
        validation_results.passed = false;
    end
    mpta_points = sum(strcmp(global_output.control_mode, 'MPTA'));
    if mpta_points == 0 && max(global_output.P_electrical) > P_nominal * 0.7
        validation_results.warnings{end+1} = 'No MPTA operation detected';
    end
    power_violations = global_output.P_electrical > global_output.P_generator * 1.01;
    if any(power_violations)
        validation_results.warnings{end+1} = 'Electrical power exceeds mechanical input';
        validation_results.passed = false;
    end
    if any(global_output.Cp > 0.60)
        validation_results.warnings{end+1} = 'Cp exceeds realistic limits';
    end
    check_smooth_transitions(global_output);
end

function damage_validation = validate_damage_accumulation(global_output, I_nominal)
    damage_validation.passed = false;
    damage_validation.messages = {};
    operational_indices = global_output.V0 >= 4 & global_output.V0 <= 12;
    operational_damage = global_output.damage(operational_indices);
    if max(operational_damage) == 0
        damage_validation.messages{end+1} = 'No damage during normal operation';
        return;
    end
    damage_rate = max(diff(global_output.damage(operational_indices)));
    if damage_rate > 0.01
        damage_validation.messages{end+1} = sprintf('Damage rate too high: %.4f', damage_rate);
        return;
    elseif damage_rate < 1e-8
        damage_validation.messages{end+1} = 'Damage rate too low';
    end
    high_current_indices = global_output.I_phase > I_nominal * 0.5;
    high_current_stress = global_output.stress_range(high_current_indices);
    if mean(high_current_stress) < 0.5
        damage_validation.messages{end+1} = 'Stress too low at high current';
        return;
    end
    if isempty(damage_validation.messages)
        damage_validation.passed = true;
        damage_validation.messages{1} = sprintf('Validated: max=%.6f, final=%.6f', ...
            max(global_output.damage), global_output.damage(end));
    end
end

function check_smooth_transitions(global_output)
    torque_gradients = diff(global_output.T_em);
    abrupt_changes = find(abs(torque_gradients) > 2.0);
    if ~isempty(abrupt_changes)
        fprintf('⚠ %d abrupt torque changes\n', length(abrupt_changes));
    end
    modes = global_output.control_mode;
    fprintf('Control mode analysis:\n');
    fprintf('  MPPT: %d\n', sum(strcmp(modes, 'MPPT')));
    fprintf('  MPTA: %d\n', sum(strcmp(modes, 'MPTA')));
    fprintf('  CUT-IN: %d\n', sum(strcmp(modes, 'CUT-IN')));
end

function cfd_validation = validate_cfd_readiness(global_output)
    cfd_validation.passed = true;
    cfd_validation.messages = {};
    cp_smoothness = std(diff(global_output.Cp(global_output.V0 > 3)));
    if cp_smoothness > 0.08
        cfd_validation.messages{end+1} = sprintf('Cp variability (std=%.3f)', cp_smoothness);
        cfd_validation.passed = false;
    end
    a_max = max(global_output.a(:));
    a_min = min(global_output.a(:));
    if a_max > 0.6 || a_min < -0.1
        cfd_validation.messages{end+1} = sprintf('Induction factors out: [%.3f, %.3f]', a_min, a_max);
        cfd_validation.passed = false;
    end
    for i = 1:size(global_output.a, 1)
        a_profile = global_output.a(i, :);
        if any(diff(a_profile(2:end)) > 0.3)
            cfd_validation.messages{end+1} = sprintf('Non-physical a at V0=%.1f', global_output.V0(i));
            cfd_validation.passed = false;
            break;
        end
    end
end

function [id_fw, iq_fw] = apply_field_weakening_analytical(id, iq, Rs, Ld, Lq, lambda_m, omega_e, V_max, I_max)
    id_fw = id;
    iq_fw = iq;

    if omega_e < 1e-3
        return;
    end

    % ===== NEW: early exit when voltage is well within limits =====
    vd_test = Rs * id_fw - omega_e * Lq * iq_fw;
    vq_test = Rs * iq_fw + omega_e * (Ld * id_fw + lambda_m);
    V_test  = sqrt(vd_test^2 + vq_test^2);
    if V_test <= V_max * 0.98
        return;   % No field weakening needed
    end

    V_lim = V_max / omega_e;
    a = Lq^2 - Ld^2;
    b = 2 * Ld * lambda_m;
    c = lambda_m^2 - V_lim^2;
    if abs(a) > 1e-6
        discriminant = b^2 - 4*a*c;
        if discriminant >= 0
            id_candidate = (-b - sqrt(discriminant)) / (2*a);
            if abs(id_candidate) <= I_max
                iq_candidate = sqrt(I_max^2 - id_candidate^2);
                V_test = sqrt((Lq*iq_candidate)^2 + (Ld*id_candidate + lambda_m)^2);
                if V_test <= V_lim * 1.01
                    id_fw = id_candidate;
                    iq_fw = iq_candidate;
                    return;
                end
            end
        end
    end
    max_iter = 20;
    for iter = 1:max_iter
        term = (V_max/omega_e)^2 - (Lq*iq_fw)^2;
        if term > 0
            id_fw = (sqrt(term) - lambda_m) / Ld;
        else
            id_fw = -lambda_m/Ld;
        end
        I_mag = sqrt(id_fw^2 + iq_fw^2);
        if I_mag > I_max
            scaling = I_max / I_mag;
            id_fw = id_fw * scaling;
            iq_fw = iq_fw * scaling;
        end
        Vd = Rs * id_fw - omega_e * Lq * iq_fw;
        Vq = Rs * iq_fw + omega_e * (Ld * id_fw + lambda_m);
        V_mag = sqrt(Vd^2 + Vq^2);
        if V_mag <= V_max * 1.01
            break;
        end
        if iter > 10
            iq_fw = iq_fw * 0.95;
        end
    end
end

function P_iron = calculate_enhanced_iron_losses(omega_e, lambda_m, id, iq, Ld, Lq, K_h, K_e, T_winding)
    psi_d = Ld * id + lambda_m;
    psi_q = Lq * iq;
    psi_mag = sqrt(psi_d^2 + psi_q^2);
    T_factor = 1 + 0.005 * (T_winding - 25);
    P_hyst = K_h * omega_e * psi_mag^2 * T_factor;
    P_eddy = K_e * omega_e^2 * psi_mag^2 * T_factor;
    P_excess = 0.015 * omega_e^1.5 * psi_mag^1.5 * T_factor;
    P_iron = P_hyst + P_eddy + P_excess;
end

function [id_opt, iq_opt] = calculate_mtpa_currents(I_max, lambda_m, Ld, Lq, P)
    if abs(Ld - Lq) < 1e-6
        id_opt = 0;
        iq_opt = I_max;
    else
        fun = @(iq) -calculate_torque_from_iq(iq, I_max, lambda_m, Ld, Lq, P);
        iq_opt = fminbnd(fun, 0, I_max, optimset('Display','off'));
        id_opt = -sqrt(I_max^2 - iq_opt^2);
    end
end

function T = calculate_torque_from_iq(iq, I_max, lambda_m, Ld, Lq, P)
    if iq > I_max
        T = 0;
    else
        id = -sqrt(I_max^2 - iq^2);
        T = (3/2) * P * (lambda_m * iq + (Ld - Lq) * id * iq);
    end
end

function [Q_rotor, T_total, alpha, phi, F, CL, CD, eff, Re_local_vec] = ...
    calculate_performance(r, a, ap, V0, omega, beta_deg, B, R, min_r_ratio, ...
                          F_CL, F_CD, stall_alpha_pos, stall_alpha_neg, f_stall, dr, chord, polar, mu, rho, use_Re_scaling, Re_ref)
    n = length(r);
    Q_rotor = 0;
    T_total = 0;
    alpha = zeros(n,1);
    phi = zeros(n,1);
    F = zeros(n,1);
    CL = zeros(n,1);
    CD = zeros(n,1);
    eff = zeros(n,1);
    Re_local_vec = zeros(n,1);
    CL_max = max(polar.CL);
    CL_min = min(polar.CL);
    for i = 1:n
        if r(i) < min_r_ratio * R
            continue;
        end
        Vax = V0 * (1 - a(i));
        Vtan = omega * r(i) * (1 + ap(i));
        phi_i = atan2(Vax, Vtan);
        alpha_i = rad2deg(phi_i) - beta_deg(i);
        alpha_i = max(-20, min(20, alpha_i));
        g_tip = (B/2) * (R - r(i)) / (r(i) * abs(sin(phi_i)));
        g_root = (B/2) * (r(i) - min_r_ratio*R) / (r(i) * abs(sin(phi_i)));
        if a(i) < 0.5
            factor_tip = 1 + (0.1*B)/(1 - 2*a(i));
            factor_root = 1 + (0.1*B)/(1 - 2*a(i));
        else
            factor_tip = 1;
            factor_root = 1;
        end
        F_tip = (2/pi) * acos(exp(-g_tip * factor_tip));
        F_root = (2/pi) * acos(exp(-g_root * factor_root));
        F_i = F_tip * F_root;
        W = sqrt(Vax^2 + Vtan^2);
        Re_local = (rho * W * chord(i)) / mu;
        Re_local_vec(i) = Re_local;
        if use_Re_scaling
            CL_max_adj = CL_max * (Re_local / Re_ref)^0.1;
            CL_min_adj = CL_min * (Re_local / Re_ref)^0.1;
        else
            CL_max_adj = CL_max;
            CL_min_adj = CL_min;
        end
        if alpha_i > 15
            viterna_vals = [
                (CL_max_adj/2) * sind(2*alpha_i),
                CL_max_adj * (sind(alpha_i))^2
            ];
            CL_i = viterna_vals(1);
            CD_i = viterna_vals(2);
            CL_i = f_stall(i) * F_CL(15) + (1 - f_stall(i)) * CL_i;
            CD_i = f_stall(i) * F_CD(15) + (1 - f_stall(i)) * CD_i;
        elseif alpha_i < -15
            viterna_vals = [
                (abs(CL_min_adj)/2 * sind(2*alpha_i)),
                abs(CL_min_adj) * (sind(alpha_i))^2
            ];
            CL_i = viterna_vals(1);
            CD_i = viterna_vals(2);
            CL_i = f_stall(i) * F_CL(-15) + (1 - f_stall(i)) * CL_i;
            CD_i = f_stall(i) * F_CD(-15) + (1 - f_stall(i)) * CD_i;
        else
            CL_i = F_CL(alpha_i);
            CD_i = F_CD(alpha_i);
        end
        alpha(i) = alpha_i;
        phi(i) = rad2deg(phi_i);
        F(i) = F_i;
        CL(i) = CL_i;
        CD(i) = CD_i;
        if CD_i > 0
            eff(i) = CL_i / CD_i;
        else
            eff(i) = 0;
        end
        Fn = 0.5 * rho * W^2 * chord(i) * (CL_i * cos(phi_i) + CD_i * sin(phi_i));
        Ft = 0.5 * rho * W^2 * chord(i) * (CL_i * sin(phi_i) - CD_i * cos(phi_i));
        dT = Fn * dr(i) * B;
        dQ = Ft * r(i) * dr(i) * B;
        T_total = T_total + dT;
        Q_rotor = Q_rotor + dQ;
    end
end