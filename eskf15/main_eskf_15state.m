clear; clc; close all;

%% ============================================================
% SYNTHETIC FIXED-WING-LIKE DATA MODE - 18 STATE ESKF
% ============================================================
%
% State:
%   nominal:
%     p_n, v_n, q_nb, b_g, b_a, b_baro, wind_ne
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%      = 18x1
%
% Pipeline:
%   - fixed-wing-like synthetic truth
%   - NED frame
%   - IMU + GPS pos + GPS vel + baro + TAS
%   - wind_ref truth and TAS measurement
%   - GPS schedule / dropout support
%   - baro bias estimated only when recent GPS height exists
%   - GPS outage: nominal b_baro freeze, no manual P manipulation
%   - GPS outage: process noise scaling, covariance grows through Qd
%   - innovation gates
% ============================================================

%% ============================================================
% 1) Init filter
% ============================================================

[state, ~, params] = init_filter();

%% ---------------- Script paths ----------------
scriptDir = string(fileparts(mfilename("fullpath")));
if strlength(scriptDir) == 0
    scriptDir = string(pwd);
end

%% ---------------- USER FLAGS ----------------
USE_BARO = true;
USE_GPS_VEL = true;
USE_AIRSPEED = true;

USE_TRUTH_ATT_INIT_FOR_DEBUG = true;
USE_GPS_VEL_INIT = true;

USE_OBSERVABILITY_ANALYSIS = true;

SAVE_RUN_OUTPUT = true;
save_file = "synthetic_fixedwing_like_run.mat";

SAVE_FIGURES = true;             % true: figure export aktif
SAVE_ONLY_FINAL_FIGS = false;    % true: sadece listelenen final figure'lar kaydedilir
SAVE_FIG_FORMATS = ["png"];      % sentetik main icin varsayilan PNG

FINAL_FIGURES_ALWAYS = [ ...
    "3D Trajectory"
    "Position NED"
    "Velocity NED"
    "Attitude"
    "Gyro Bias"
    "Accel Bias"
    "Covariance Diagonal Groups"
    "Measurement Residuals"
    "Observation Consistency Errors"
    "Observability Analysis"
    "Synthetic GPS Schedule and Baro Bias"
    "Synthetic Innovation Gate Flags"
    "Synthetic Estimated Wind"
    "Synthetic Airspeed Residual"
];

%% ---------------- GPS MEASUREMENT SCHEDULE ----------------
% "always"  : GPS position + GPS velocity tüm koşu boyunca açık
% "never"   : GPS position + GPS velocity tüm koşu boyunca kapalı
% "windows" : sadece verilen zaman pencerelerinde açık
GPS_MEAS_MODE = "always";

% Dropout denemesi:
% GPS_MEAS_MODE = "windows";
% GPS_MEAS_WINDOWS = [
%     0    120
%     150  375
% ];
GPS_MEAS_WINDOWS = [
    0 inf
];

%% ---------------- BARO BIAS OBSERVABILITY ----------------
ESTIMATE_BARO_BIAS_WITH_GPS = true;
GPS_HEIGHT_ANCHOR_TIMEOUT_S = 0.50;

%% ---------------- GPS OUTAGE PROCESS-NOISE HANDLING ----------------
% P'ye elle müdahale yok.
% GPS yokken covariance Qd üzerinden büyüsün diye process noise artırılır.
GPS_OUTAGE_SIGMA_A_SCALE     = 3.0;
GPS_OUTAGE_SIGMA_G_SCALE     = 1.5;
GPS_OUTAGE_SIGMA_BA_RW_SCALE = 4.0;

%% ============================================================
% 1.5) Synthetic fixed-wing-like 18-state tuning
% ============================================================
%
% Bu main dosyası, patched simulate_truth_and_sensors.m içindeki
% generate_fixedwing_like_truth() fonksiyonuyla beraber kullanılacak şekilde
% ayarlandı.
%
% Hedef:
%   - Gerçek logun birebir kopyası değil.
%   - Fixed-wing için mantıklı TAS / roll / pitch / wind ölçeği.
%   - Wind + TAS observability için yeterli dönüş ve hız değişimi.
%   - Sensör noise tarafında filtreyi gereksiz yere dövmeyen ama gerçekçi R.

params.sensor_profile = 'Synthetic Fixed-Wing-Like: IMU + GPS + baro + TAS';

%% Frame / timing
params.g_n = [0; 0; 9.81];       % NED: D positive down
params.sim.dt = 0.004;           % 250 Hz IMU-like truth grid
params.sim.T_end = 375.0;

%% Required truth defaults
params.truth.bg0 = [0; 0; 0];
params.truth.ba0 = [0; 0; 0];

%% Fixed-wing-like truth scenario
params.truth.scenario_mode = "fixedwing_like";

% Initial position: [N; E; D] m. D negative means above local origin.
% Real logun birebir başlangıcı değil; sadece makul uçuş irtifası.
params.truth.fixedwing.p0_ned = [0.0; 0.0; -45.0];

% Airspeed profile bounds. Patched truth generator Va'yı bu aralıkta clamp eder.
params.truth.fixedwing.Va_trim = 14.5;     % [m/s]
params.truth.fixedwing.Va_min  = 10.5;     % [m/s]
params.truth.fixedwing.Va_max  = 18.5;     % [m/s]

% Heading / maneuver profile. Generator içinde loiter/racetrack/S-turn benzeri
% schedule üretir; bu değerler sınır ve başlangıç yönünü belirler.
params.truth.fixedwing.yaw0_deg = -80.0;
params.truth.fixedwing.max_yaw_rate_deg_s = 32.0;

% Vertical motion / bank limits.
params.truth.fixedwing.max_gamma_deg = 8.0;
params.truth.fixedwing.max_roll_deg  = 45.0;

%% IMU noise, fixed-wing real-like
params.sigma_g = deg2rad(0.20);
params.sigma_a = 0.080;

params.sigma_bg_rw = deg2rad(0.010);
params.sigma_ba_rw = 0.010;

%% Synthetic sensor rates
params.synthetic_gps_rate_hz = 5;
params.synthetic_baro_rate_hz = 20;
params.synthetic_airspeed_rate_hz = 20;

%% GPS noise, fixed-wing real-like
% Eski [3 3 5] / [0.50 0.50 0.80] çok kaba kalıyordu.
% Hala noisy, ama filtreyi sensör gürültüsüyle gereksiz yere boğmuyor.
params.sigma_gps_pos = [2.0; 2.0; 3.0];
params.R_gps_pos = diag(params.sigma_gps_pos.^2);

params.sigma_gps_vel = [0.35; 0.35; 0.50];
params.R_gps_vel = diag(params.sigma_gps_vel.^2);

%% Baro noise / bias
params.sigma_baro = 1.5;
params.R_baro = params.sigma_baro^2;

% True synthetic baro bias.
params.truth.b_baro0 = 6.0;
params.truth.b_baro_rw = 0.002;

% Filter-side baro bias uncertainty / random walk.
params.sigma_baro_bias0 = 10.0;
params.sigma_baro_bias_rw = 0.003;

params.max_baro_update_rate_hz = 20;
params.estimate_baro_bias = false;

%% Wind + TAS
params.synthetic.enable_wind_tas = USE_AIRSPEED;
params.synthetic.wind_seed = 42;

% Gerçek fixed-wing log ölçeğine yakın ama aynı olmayan nominal wind.
% Patched truth generator bu base etrafında yavaş değişen wind üretir.
params.synthetic.wind_base_ne = [0.65; -0.85];

% TAS measurement noise. Eski 2.5 m/s yüzünden 14 m/s truth'tan 4 m/s
% gibi uçak-olmayan measurement çıkabiliyordu. Buna gerek yok, medeniyet çöktü zaten.
params.sigma_tas = 1.0;
params.R_tas = params.sigma_tas^2;

params.max_airspeed_update_rate_hz = 20;

% Wind state prior / random walk
params.sigma_wind0 = [4.0; 4.0];
params.sigma_wind_rw = 0.020;

params.use_joseph_form = true;

%% Innovation gates
params.gps_pos_gate_chi2 = 25.0;
params.gps_vel_gate_chi2 = 25.0;
params.baro_gate_chi2    = 9.0;
params.tas_gate_chi2     = 9.0;

%% Observability options
params.obs.enable = USE_OBSERVABILITY_ANALYSIS;
params.obs.sample_step = 100;
params.obs.window_sec = 5.0;
params.obs.rank_tol = 1e-8;

%% ============================================================
% 2) Generate synthetic sim
% ============================================================

fprintf('[main synthetic] Fixed-wing-like sentetik simulasyon modu.\n');

sim = simulate_truth_and_sensors(params);

t = sim.t;
N = numel(t);

fprintf("\n=== MAIN SYNTHETIC FIXED-WING-LIKE 18-STATE ESKF ===\n");
fprintf("N = %d samples\n", N);
fprintf("t range = %.3f to %.3f s\n", t(1), t(end));

fprintf('\n--- GPS MEASUREMENT WINDOWS ---\n');
disp(GPS_MEAS_WINDOWS);

tmp_use_gps = false(1,N);

for kk = 1:N
    tmp_use_gps(kk) = is_gps_schedule_active(sim.t(kk), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);
end

sim.gps_meas_scheduled = tmp_use_gps;

fprintf('Expected GPS active time [s]   : %.3f\n', sum(tmp_use_gps) * median(diff(sim.t), "omitnan"));
fprintf('Expected GPS active samples    : %d / %d\n', sum(tmp_use_gps), N);

%% ============================================================
% 3) Initial nominal state
% ============================================================

idx_gps0 = find(sim.gps_pos_available & all(isfinite(sim.gps_pos),1), 1, 'first');

if isempty(idx_gps0)
    warning('[main synthetic] GPS ölçümü yok, p0 = truth p0 kullanılıyor.');
    state.p_n = sim.p_true(:,1);
else
    state.p_n = sim.gps_pos(:,idx_gps0);
end

idx_gps_vel0 = find(sim.gps_vel_available & all(isfinite(sim.gps_vel),1), 1, 'first');

if USE_GPS_VEL_INIT && ~isempty(idx_gps_vel0)
    state.v_n = sim.gps_vel(:,idx_gps_vel0);
    fprintf('[main synthetic] Initial velocity GPS velocity üzerinden başlatıldı.\n');
else
    state.v_n = [0;0;0];
    fprintf('[main synthetic] Initial velocity v0 = [0;0;0] seçildi.\n');
end

if USE_TRUTH_ATT_INIT_FOR_DEBUG
    state.q_nb = sim.q_true(:,1);
    state.q_nb = state.q_nb / norm(state.q_nb);
    fprintf('[main synthetic] DEBUG: attitude q0 truth üzerinden başlatıldı.\n');
else
    state.q_nb = [1;0;0;0];
    fprintf('[main synthetic] REAL-LIKE: attitude q0 identity seçildi.\n');
end

state.b_g = [0;0;0];
state.b_a = [0;0;0];

state.b_baro = 0;
state.wind_ne = [0;0];

fprintf('[main synthetic] Initial p = [%.3f %.3f %.3f]^T m\n', state.p_n);
fprintf('[main synthetic] Initial v = [%.3f %.3f %.3f]^T m/s\n', state.v_n);
fprintf('[main synthetic] Initial b_baro = %.6f m\n', state.b_baro);
fprintf('[main synthetic] Initial wind NE = [%.6f %.6f]^T m/s\n', state.wind_ne);
fprintf('[main synthetic] GPS_MEAS_MODE = %s\n', GPS_MEAS_MODE);

%% ============================================================
% 4) Initial covariance, 18-state
% ============================================================

idx_p      = 1:3;
idx_v      = 4:6;
idx_th     = 7:9;
idx_bg     = 10:12;
idx_ba     = 13:15;
idx_bbaro  = 16;
idx_wind   = 17:18;

sigma_p0  = [3.0; 3.0; 5.0];
sigma_v0  = [1.0; 1.0; 1.5];

if USE_TRUTH_ATT_INIT_FOR_DEBUG
    sigma_th0 = deg2rad([1.0; 1.0; 5.0]);
else
    sigma_th0 = deg2rad([5.0; 5.0; 45.0]);
end

sigma_bg0 = deg2rad([0.5; 0.5; 0.8]);
sigma_ba0 = [0.30; 0.30; 0.50];

P = zeros(18,18);

P(idx_p, idx_p)         = diag(sigma_p0.^2);
P(idx_v, idx_v)         = diag(sigma_v0.^2);
P(idx_th, idx_th)       = diag(sigma_th0.^2);
P(idx_bg, idx_bg)       = diag(sigma_bg0.^2);
P(idx_ba, idx_ba)       = diag(sigma_ba0.^2);
P(idx_bbaro, idx_bbaro) = params.sigma_baro_bias0^2;
P(idx_wind, idx_wind)   = diag(params.sigma_wind0.^2);

%% ============================================================
% 5) Logs
% ============================================================

log_p = nan(3, N);
log_v = nan(3, N);
log_q = nan(4, N);
log_bg = nan(3, N);
log_ba = nan(3, N);
log_bbaro = nan(1, N);
log_wind_ne = nan(2, N);
log_Pdiag = nan(18, N);

log_res_gps_pos = nan(3, N);
log_res_gps_vel = nan(3, N);
log_res_baro = nan(1, N);
log_res_airspeed = nan(1,N);

log_use_gps_meas = false(1,N);
log_estimate_baro_bias = false(1,N);

log_gate_gps_pos = false(1,N);
log_gate_gps_vel = false(1,N);
log_gate_baro = false(1,N);
log_gate_airspeed = false(1,N);

used_gps_pos_updates = 0;
used_gps_vel_updates = 0;
used_baro_updates = 0;
used_airspeed_updates = 0;

accepted_gps_pos_updates = 0;
accepted_gps_vel_updates = 0;
accepted_baro_updates = 0;
accepted_airspeed_updates = 0;

rejected_gps_pos_updates = 0;
rejected_gps_vel_updates = 0;
rejected_baro_updates = 0;
rejected_airspeed_updates = 0;

log_p(:,1) = state.p_n;
log_v(:,1) = state.v_n;
log_q(:,1) = state.q_nb;
log_bg(:,1) = state.b_g;
log_ba(:,1) = state.b_a;
log_bbaro(1) = state.b_baro;
log_wind_ne(:,1) = state.wind_ne;
log_Pdiag(:,1) = diag(P);

%% ============================================================
% 6) Rate limiters and dynamic state
% ============================================================

last_baro_update_t = -inf;

if isfield(params, 'max_baro_update_rate_hz') && params.max_baro_update_rate_hz > 0
    baro_min_dt = 1 / params.max_baro_update_rate_hz;
else
    baro_min_dt = 0;
end

last_airspeed_update_t = -inf;

if isfield(params, 'max_airspeed_update_rate_hz') && params.max_airspeed_update_rate_hz > 0
    airspeed_min_dt = 1 / params.max_airspeed_update_rate_hz;
else
    airspeed_min_dt = 0;
end

prev_use_gps_meas = is_gps_schedule_active(t(1), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);

last_gps_pos_update_t = -inf;
last_gps_vel_update_t = -inf;

baro_bias_hold_value = state.b_baro;

log_use_gps_meas(1) = prev_use_gps_meas;
log_estimate_baro_bias(1) = false;

%% ============================================================
% 7) Main ESKF loop
% ============================================================

fprintf("\n[main synthetic] ESKF replay başlıyor...\n");

for k = 2:N

    dt_k = sim.t(k) - sim.t(k-1);

    if dt_k <= 0 || dt_k > 0.2
        warning("[main synthetic] Bad dt at k=%d: %.6f s, propagation skipped.", k, dt_k);
        continue;
    end

    %% ---------------- GPS measurement schedule ----------------
    use_gps_meas_now = is_gps_schedule_active(sim.t(k), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);

    if prev_use_gps_meas && ~use_gps_meas_now
        baro_bias_hold_value = state.b_baro;
        fprintf('[main synthetic] t=%.2f s: GPS OFF, freezing b_baro = %.6f m\n', ...
            sim.t(k), baro_bias_hold_value);
    end

    if ~prev_use_gps_meas && use_gps_meas_now
        fprintf('[main synthetic] t=%.2f s: GPS ON, re-enabling b_baro estimation\n', ...
            sim.t(k));
    end

    if ~use_gps_meas_now
        state.b_baro = baro_bias_hold_value;
    end

    %% ---------------- Per-step params ----------------
    params_k = params;

    if ~use_gps_meas_now
        params_k.sigma_baro_bias_rw = 0.0;

        params_k.sigma_a     = params.sigma_a     * GPS_OUTAGE_SIGMA_A_SCALE;
        params_k.sigma_g     = params.sigma_g     * GPS_OUTAGE_SIGMA_G_SCALE;
        params_k.sigma_ba_rw = params.sigma_ba_rw * GPS_OUTAGE_SIGMA_BA_RW_SCALE;
    end

    params_k.estimate_baro_bias = false;

    %% ---------------- IMU propagation ----------------
    imu.gyro_m = sim.imu_gyro(:,k);
    imu.accel_m = sim.imu_accel(:,k);

    state = propagate_nominal(state, imu, params_k, dt_k);

    [F, G, Qd] = compute_F_G_Qd(state, imu, params_k, dt_k);
    P = propagate_covariance(P, F, G, Qd, dt_k);

    %% ---------------- GPS position update ----------------
    if use_gps_meas_now && sim.gps_pos_available(k)
        z_gps_pos = sim.gps_pos(:,k);

        if all(isfinite(z_gps_pos))
            [state, P, residual, ~, ~, accepted] = update_gnss_pos(state, P, z_gps_pos, params_k);

            log_res_gps_pos(:,k) = residual;
            used_gps_pos_updates = used_gps_pos_updates + 1;

            if accepted
                accepted_gps_pos_updates = accepted_gps_pos_updates + 1;
                last_gps_pos_update_t = sim.t(k);
            else
                rejected_gps_pos_updates = rejected_gps_pos_updates + 1;
                log_gate_gps_pos(k) = true;
            end
        end
    end

    %% ---------------- GPS velocity update ----------------
    if USE_GPS_VEL && use_gps_meas_now && sim.gps_vel_available(k)
        z_gps_vel = sim.gps_vel(:,k);

        if all(isfinite(z_gps_vel))
            [state, P, residual, ~, ~, accepted] = update_gnss_vel(state, P, z_gps_vel, params_k);

            log_res_gps_vel(:,k) = residual;
            used_gps_vel_updates = used_gps_vel_updates + 1;

            if accepted
                accepted_gps_vel_updates = accepted_gps_vel_updates + 1;
                last_gps_vel_update_t = sim.t(k);
            else
                rejected_gps_vel_updates = rejected_gps_vel_updates + 1;
                log_gate_gps_vel(k) = true;
            end
        end
    end

    %% ---------------- Is baro bias observable now? ----------------
    gps_height_recent = use_gps_meas_now && ...
        ((sim.t(k) - last_gps_pos_update_t) <= GPS_HEIGHT_ANCHOR_TIMEOUT_S);

    estimate_baro_bias_now = ESTIMATE_BARO_BIAS_WITH_GPS && gps_height_recent;

    params_k.estimate_baro_bias = estimate_baro_bias_now;

    if ~estimate_baro_bias_now
        params_k.sigma_baro_bias_rw = 0.0;
    end

    %% ---------------- Barometer update ----------------
    if USE_BARO && sim.baro_available(k)

        if sim.t(k) - last_baro_update_t >= baro_min_dt

            z_baro = sim.baro(k);

            if isfinite(z_baro)
                [state, P, residual, ~, ~, accepted] = update_baro(state, P, z_baro, params_k);

                log_res_baro(k) = residual;
                used_baro_updates = used_baro_updates + 1;
                last_baro_update_t = sim.t(k);

                if accepted
                    accepted_baro_updates = accepted_baro_updates + 1;
                else
                    rejected_baro_updates = rejected_baro_updates + 1;
                    log_gate_baro(k) = true;
                end
            end
        end
    end

    %% ---------------- Airspeed update, 18-state wind model ----------------
    if USE_AIRSPEED && isfield(sim, "airspeed_available") && sim.airspeed_available(k)

        if sim.t(k) - last_airspeed_update_t >= airspeed_min_dt

            z_tas = sim.airspeed(k);

            if isfinite(z_tas) && z_tas > 3
                [state, P, residual, ~, ~, accepted] = update_airspeed(state, P, z_tas, params_k);

                if isfinite(residual)
                    log_res_airspeed(k) = residual;
                    used_airspeed_updates = used_airspeed_updates + 1;
                    last_airspeed_update_t = sim.t(k);

                    if accepted
                        accepted_airspeed_updates = accepted_airspeed_updates + 1;
                    else
                        rejected_airspeed_updates = rejected_airspeed_updates + 1;
                        log_gate_airspeed(k) = true;
                    end
                end
            end
        end
    end

    %% ---------------- Enforce nominal b_baro freeze during GPS outage ----------------
    if ~use_gps_meas_now
        state.b_baro = baro_bias_hold_value;
    end

    %% ---------------- Log ----------------
    log_p(:,k) = state.p_n;
    log_v(:,k) = state.v_n;
    log_q(:,k) = state.q_nb;
    log_bg(:,k) = state.b_g;
    log_ba(:,k) = state.b_a;
    log_bbaro(k) = state.b_baro;
    log_wind_ne(:,k) = state.wind_ne;
    log_Pdiag(:,k) = diag(P);

    log_use_gps_meas(k) = use_gps_meas_now;
    log_estimate_baro_bias(k) = estimate_baro_bias_now;

    prev_use_gps_meas = use_gps_meas_now;
end

fprintf("[main synthetic] ESKF replay bitti.\n");

%% ============================================================
% 8) Plot results
% ============================================================

plot_results(sim, t, log_p, log_v, log_q, log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params);

figure('Name','Synthetic GPS Schedule and Baro Bias');

yyaxis left;
plot(t, log_bbaro, 'LineWidth', 1.3);
ylabel('b_{baro} [m]');
grid on;

yyaxis right;
stairs(t, double(log_use_gps_meas), 'LineWidth', 1.1);
hold on;
stairs(t, double(log_estimate_baro_bias), '--', 'LineWidth', 1.1);
ylabel('Flag');
ylim([-0.1 1.1]);

xlabel('Time [s]');
legend('b_{baro}', 'GPS pos+vel scheduled', 'baro bias estimated', 'Location','best');
title('Synthetic GPS Schedule and Barometer Bias Observability');

figure('Name','Synthetic Innovation Gate Flags');
plot(t, double(log_gate_gps_pos), 'LineWidth', 1.0);
hold on;
plot(t, double(log_gate_gps_vel), 'LineWidth', 1.0);
plot(t, double(log_gate_baro), 'LineWidth', 1.0);
plot(t, double(log_gate_airspeed), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Gate rejected flag');
legend('GPS pos','GPS vel','Baro','TAS','Location','best');
title('Synthetic Innovation Gate Rejections');

figure('Name','Synthetic Estimated Wind');
plot(t, log_wind_ne(1,:), 'LineWidth', 1.3);
hold on;
plot(t, log_wind_ne(2,:), 'LineWidth', 1.3);

if isfield(sim, "wind_ref_available") && any(sim.wind_ref_available)
    plot(t(sim.wind_ref_available), sim.wind_ref(1,sim.wind_ref_available), '--', 'LineWidth', 1.1);
    plot(t(sim.wind_ref_available), sim.wind_ref(2,sim.wind_ref_available), '--', 'LineWidth', 1.1);
    legend('ESKF wind_N','ESKF wind_E','True wind_N','True wind_E','Location','best');
else
    legend('ESKF wind_N','ESKF wind_E','Location','best');
end

grid on;
xlabel('Time [s]');
ylabel('Wind [m/s]');
title('Synthetic Wind Estimation');

figure('Name','Synthetic Airspeed Residual');
idx_air = isfinite(log_res_airspeed);
plot(t(idx_air), log_res_airspeed(idx_air), '.');
grid on;
xlabel('Time [s]');
ylabel('TAS residual [m/s]');
title('Synthetic TAS Innovation');

%% ============================================================
% 9) Observability analysis
% ============================================================

if isfield(params, 'obs') && isfield(params.obs, 'enable') && params.obs.enable
    obs = analyze_observability(sim, params);

    disp('--- Observability summary ---');
    fprintf('State dimension        : %d\n', obs.n_err);
    fprintf('Mean local rank        : %.2f\n', mean(obs.local_rank, 'omitnan'));
    fprintf('Mean windowed rank     : %.2f\n', mean(obs.window_rank, 'omitnan'));
    fprintf('Min local sigma_min    : %.3e\n', min(obs.local_sigma_min));
    fprintf('Min windowed sigma_min : %.3e\n', min(obs.window_sigma_min));
end

%% ============================================================
% 10) Error metrics
% ============================================================

err_p = log_p - sim.p_true;
err_v = log_v - sim.v_true;

pos_err_norm = vecnorm(err_p, 2, 1);
vel_err_norm = vecnorm(err_v, 2, 1);

idx_gps_pos_all = sim.gps_pos_available & all(isfinite(sim.gps_pos),1);
idx_gps_vel_all = sim.gps_vel_available & all(isfinite(sim.gps_vel),1);

idx_gps_pos_used = idx_gps_pos_all & log_use_gps_meas;
idx_gps_pos_holdout = idx_gps_pos_all & ~log_use_gps_meas;

idx_gps_vel_used = idx_gps_vel_all & log_use_gps_meas;
idx_gps_vel_holdout = idx_gps_vel_all & ~log_use_gps_meas;

fprintf('\n--- SYNTHETIC FIXED-WING-LIKE RUN FINISHED ---\n');

fprintf('Mean pos error norm [m]  : %.4f\n', mean(pos_err_norm, 'omitnan'));
fprintf('Final pos error norm [m] : %.4f\n', norm(err_p(:,end)));
fprintf('Max pos error norm [m]   : %.4f\n', max(pos_err_norm, [], 'omitnan'));

fprintf('Mean vel error norm [m/s]  : %.4f\n', mean(vel_err_norm, 'omitnan'));
fprintf('Final vel error norm [m/s] : %.4f\n', norm(err_v(:,end)));
fprintf('Max vel error norm [m/s]   : %.4f\n', max(vel_err_norm, [], 'omitnan'));

if any(idx_gps_pos_used)
    gps_used_err = log_p(:,idx_gps_pos_used) - sim.gps_pos(:,idx_gps_pos_used);

    fprintf('GPS position USED RMSE norm [m]    : %.4f\n', ...
        sqrt(mean(sum(gps_used_err.^2,1), 'omitnan')));
end

if any(idx_gps_pos_holdout)
    gps_hold_err = log_p(:,idx_gps_pos_holdout) - sim.gps_pos(:,idx_gps_pos_holdout);

    fprintf('GPS position HOLDOUT RMSE norm [m] : %.4f\n', ...
        sqrt(mean(sum(gps_hold_err.^2,1), 'omitnan')));
end

if any(idx_gps_vel_used)
    gpsv_used_err = log_v(:,idx_gps_vel_used) - sim.gps_vel(:,idx_gps_vel_used);

    fprintf('GPS velocity USED RMSE norm [m/s]    : %.4f\n', ...
        sqrt(mean(sum(gpsv_used_err.^2,1), 'omitnan')));
end

if any(idx_gps_vel_holdout)
    gpsv_hold_err = log_v(:,idx_gps_vel_holdout) - sim.gps_vel(:,idx_gps_vel_holdout);

    fprintf('GPS velocity HOLDOUT RMSE norm [m/s] : %.4f\n', ...
        sqrt(mean(sum(gpsv_hold_err.^2,1), 'omitnan')));
end

%% Baro residual check
idx_baro = sim.baro_available & isfinite(sim.baro);

if any(idx_baro)
    baro_raw_err = log_p(3,idx_baro) - sim.baro(idx_baro);
    baro_model_err = log_p(3,idx_baro) + log_bbaro(idx_baro) - sim.baro(idx_baro);

    fprintf('Baro raw RMSE [m]       : %.4f\n', sqrt(mean(baro_raw_err.^2, 'omitnan')));
    fprintf('Baro raw mean [m]       : %.4f\n', mean(baro_raw_err, 'omitnan'));
    fprintf('Baro modeled RMSE [m]   : %.4f\n', sqrt(mean(baro_model_err.^2, 'omitnan')));
    fprintf('Baro modeled mean [m]   : %.4f\n', mean(baro_model_err, 'omitnan'));
    fprintf('Baro modeled std [m]    : %.4f\n', std(baro_model_err, 0, 'omitnan'));
end

%% Wind reference comparison
wind_ref_rmse = nan(2,1);
wind_ref_rmse_norm = nan;
wind_ref_mean_error = nan(2,1);

if isfield(sim, "wind_ref_available") && any(sim.wind_ref_available)
    idx_wind_ref = sim.wind_ref_available & all(isfinite(sim.wind_ref),1);

    if any(idx_wind_ref)
        wind_err = log_wind_ne(:,idx_wind_ref) - sim.wind_ref(:,idx_wind_ref);

        wind_ref_rmse = sqrt(mean(wind_err.^2, 2, 'omitnan'));
        wind_ref_rmse_norm = sqrt(mean(sum(wind_err.^2,1), 'omitnan'));
        wind_ref_mean_error = mean(wind_err, 2, 'omitnan');

        fprintf('Wind RMSE N/E [m/s]       : [%.4f %.4f]\n', wind_ref_rmse);
        fprintf('Wind RMSE norm [m/s]      : %.4f\n', wind_ref_rmse_norm);
        fprintf('Wind mean error N/E [m/s] : [%.4f %.4f]\n', wind_ref_mean_error);
    end
end

idx_air = isfinite(log_res_airspeed);

if any(idx_air)
    fprintf('Airspeed residual RMSE [m/s] : %.4f\n', ...
        sqrt(mean(log_res_airspeed(idx_air).^2, 'omitnan')));
    fprintf('Airspeed residual mean [m/s] : %.4f\n', ...
        mean(log_res_airspeed(idx_air), 'omitnan'));
    fprintf('Airspeed residual std [m/s]  : %.4f\n', ...
        std(log_res_airspeed(idx_air), 0, 'omitnan'));
end

fprintf('\n--- UPDATE COUNTS ---\n');

fprintf('GPS pos updates available : %d\n', sum(sim.gps_pos_available));
fprintf('GPS pos updates used      : %d\n', used_gps_pos_updates);
fprintf('GPS pos updates accepted  : %d\n', accepted_gps_pos_updates);
fprintf('GPS pos updates rejected  : %d\n', rejected_gps_pos_updates);

fprintf('GPS vel updates available : %d\n', sum(sim.gps_vel_available));
fprintf('GPS vel updates used      : %d\n', used_gps_vel_updates);
fprintf('GPS vel updates accepted  : %d\n', accepted_gps_vel_updates);
fprintf('GPS vel updates rejected  : %d\n', rejected_gps_vel_updates);

fprintf('BARO updates available    : %d\n', sum(sim.baro_available & isfinite(sim.baro)));
fprintf('BARO updates used         : %d\n', used_baro_updates);
fprintf('BARO updates accepted     : %d\n', accepted_baro_updates);
fprintf('BARO updates rejected     : %d\n', rejected_baro_updates);

fprintf('Airspeed updates available : %d\n', sum(sim.airspeed_available));
fprintf('Airspeed updates used      : %d\n', used_airspeed_updates);
fprintf('Airspeed updates accepted  : %d\n', accepted_airspeed_updates);
fprintf('Airspeed updates rejected  : %d\n', rejected_airspeed_updates);

fprintf('GPS pos skipped by schedule : %d\n', sum(idx_gps_pos_holdout));
fprintf('GPS vel skipped by schedule : %d\n', sum(idx_gps_vel_holdout));
fprintf('GPS measurement scheduled samples : %d\n', sum(log_use_gps_meas));
fprintf('Baro bias estimated samples       : %d\n', sum(log_estimate_baro_bias));

fprintf('\n--- FINAL STATES ---\n');

fprintf('Final estimated p [m]       : [%.4f %.4f %.4f]^T\n', log_p(:,end));
fprintf('Final truth p [m]           : [%.4f %.4f %.4f]^T\n', sim.p_true(:,end));
fprintf('Final estimated v [m/s]     : [%.4f %.4f %.4f]^T\n', log_v(:,end));
fprintf('Final truth v [m/s]         : [%.4f %.4f %.4f]^T\n', sim.v_true(:,end));
fprintf('Final gyro bias [rad/s]     : [%.6f %.6f %.6f]^T\n', log_bg(:,end));
fprintf('Final accel bias [m/s^2]    : [%.6f %.6f %.6f]^T\n', log_ba(:,end));
fprintf('Final baro offset [m]       : %.6f\n', log_bbaro(end));
fprintf('Final true baro offset [m]  : %.6f\n', sim.b_baro_true(end));
fprintf('Final wind NE [m/s]         : [%.6f %.6f]^T\n', log_wind_ne(:,end));

if isfield(sim, 'wind_ref')
    fprintf('Final true wind NE [m/s]    : [%.6f %.6f]^T\n', sim.wind_ref(:,end));
end

fprintf('\n--- CONFIG USED ---\n');

fprintf('Sensor profile : %s\n', params.sensor_profile);
fprintf('USE_BARO       : %d\n', USE_BARO);
fprintf('USE_GPS_VEL    : %d\n', USE_GPS_VEL);
fprintf('USE_AIRSPEED   : %d\n', USE_AIRSPEED);
fprintf('USE_TRUTH_ATT_INIT_DEBUG : %d\n', USE_TRUTH_ATT_INIT_FOR_DEBUG);
fprintf('GPS_MEAS_MODE  : %s\n', GPS_MEAS_MODE);
fprintf('GPS_HEIGHT_ANCHOR_TIMEOUT_S : %.3f\n', GPS_HEIGHT_ANCHOR_TIMEOUT_S);
fprintf('GPS_OUTAGE_SIGMA_A_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_A_SCALE);
fprintf('GPS_OUTAGE_SIGMA_G_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_G_SCALE);
fprintf('GPS_OUTAGE_SIGMA_BA_RW_SCALE : %.3f\n', GPS_OUTAGE_SIGMA_BA_RW_SCALE);
fprintf('sigma_gps_pos : [%.3f %.3f %.3f]\n', params.sigma_gps_pos);
fprintf('sigma_gps_vel : [%.3f %.3f %.3f]\n', params.sigma_gps_vel);
fprintf('sigma_baro    : %.3f\n', params.sigma_baro);
fprintf('sigma_tas     : %.3f\n', params.sigma_tas);

%% ============================================================
% 11) Save run output
% ============================================================

if SAVE_RUN_OUTPUT
    fprintf("\n[main synthetic] Saving run output to %s ...\n", save_file);

    save(save_file, ...
        "sim", "t", ...
        "log_p", "log_v", "log_q", ...
        "log_bg", "log_ba", "log_bbaro", "log_wind_ne", "log_Pdiag", ...
        "log_res_gps_pos", "log_res_gps_vel", "log_res_baro", "log_res_airspeed", ...
        "log_use_gps_meas", "log_estimate_baro_bias", ...
        "log_gate_gps_pos", "log_gate_gps_vel", "log_gate_baro", "log_gate_airspeed", ...
        "used_gps_pos_updates", "used_gps_vel_updates", "used_baro_updates", "used_airspeed_updates", ...
        "accepted_gps_pos_updates", "accepted_gps_vel_updates", "accepted_baro_updates", "accepted_airspeed_updates", ...
        "rejected_gps_pos_updates", "rejected_gps_vel_updates", "rejected_baro_updates", "rejected_airspeed_updates", ...
        "wind_ref_rmse", "wind_ref_rmse_norm", "wind_ref_mean_error", ...
        "USE_BARO", "USE_GPS_VEL", "USE_AIRSPEED", "USE_TRUTH_ATT_INIT_FOR_DEBUG", ...
        "GPS_MEAS_MODE", "GPS_MEAS_WINDOWS", ...
        "GPS_HEIGHT_ANCHOR_TIMEOUT_S", ...
        "ESTIMATE_BARO_BIAS_WITH_GPS", ...
        "GPS_OUTAGE_SIGMA_A_SCALE", "GPS_OUTAGE_SIGMA_G_SCALE", ...
        "GPS_OUTAGE_SIGMA_BA_RW_SCALE", ...
        "SAVE_FIGURES", "SAVE_ONLY_FINAL_FIGS", "SAVE_FIG_FORMATS", ...
        "params", "-v7.3");

    fprintf("[main synthetic] Saved: %s\n", save_file);
end

%% ============================================================
% 12) Save figures, optional
% ============================================================

if SAVE_FIGURES
    figRoot = fullfile(scriptDir, "figures_saved");
    runTag = "synthetic_fixedwing_eskf_" + string(GPS_MEAS_MODE) + "_" + ...
        string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    figDir = fullfile(figRoot, runTag);

    finalFigureNames = FINAL_FIGURES_ALWAYS;

    if SAVE_ONLY_FINAL_FIGS
        close_non_final_figures(finalFigureNames);
    end

    save_final_open_figures( ...
        figDir, ...
        char("synthetic_fixedwing_eskf_" + string(GPS_MEAS_MODE)), ...
        SAVE_FIG_FORMATS);

    fprintf("[main synthetic] Figures saved to:\n%s\n", figDir);
else
    fprintf("[main synthetic] SAVE_FIGURES = false, figure export skipped.\n");
end

%% ============================================================
% Local helper: GPS schedule
% ============================================================

function active = is_gps_schedule_active(t_now, mode, windows)
    switch mode
        case "always"
            active = true;

        case "never"
            active = false;

        case "windows"
            active = false;

            if isempty(windows)
                return;
            end

            for iw = 1:size(windows,1)
                if t_now >= windows(iw,1) && t_now <= windows(iw,2)
                    active = true;
                    return;
                end
            end

        otherwise
            error('Unknown GPS_MEAS_MODE: %s', mode);
    end
end

function close_non_final_figures(keepNames)
    %CLOSE_NON_FINAL_FIGURES
    % keepNames icinde olmayan figure'lari kapatir.
    keepNames = string(keepNames);
    figs = findobj(groot, "Type", "figure");

    for i = 1:numel(figs)
        fig = figs(i);

        if isempty(fig) || ~isvalid(fig)
            continue;
        end

        try
            figName = string(fig.Name);
        catch
            figName = "";
        end

        keep = false;

        for k = 1:numel(keepNames)
            if contains(lower(figName), lower(keepNames(k)))
                keep = true;
                break;
            end
        end

        if ~keep
            fprintf("[main synthetic] Closing non-final figure: %s\n", figName);
            close(fig);
        end
    end
end

function save_final_open_figures(outDir, prefix, formats)
    %SAVE_FINAL_OPEN_FIGURES
    % Acik kalan figure'lari istenen formatlarda kaydeder.

    if nargin < 1 || strlength(string(outDir)) == 0
        outDir = fullfile(pwd, "figures_saved");
    end

    if nargin < 2 || strlength(string(prefix)) == 0
        prefix = "fig";
    end

    if nargin < 3 || isempty(formats)
        formats = ["png"];
    end

    formats = string(formats);

    if ~exist(outDir, "dir")
        mkdir(outDir);
    end

    figs = findobj(groot, "Type", "figure");

    if isempty(figs)
        warning("Kaydedilecek acik figure bulunamadi.");
        return;
    end

    figNums = nan(numel(figs), 1);

    for i = 1:numel(figs)
        if isvalid(figs(i)) && isprop(figs(i), "Number")
            figNums(i) = figs(i).Number;
        else
            figNums(i) = inf;
        end
    end

    [~, idx] = sort(figNums, "ascend");
    figs = figs(idx);

    fprintf("\n=== Saving figures to: %s ===\n", outDir);
    savedCount = 0;

    for i = 1:numel(figs)
        fig = figs(i);

        if isempty(fig) || ~isvalid(fig)
            continue;
        end

        try
            figure(fig);
            drawnow;
        catch
            continue;
        end

        try
            set(fig, "Color", "w");
            set(fig, "InvertHardcopy", "off");
            set(fig, "PaperPositionMode", "auto");
        catch
        end

        try
            if isprop(fig, "Name") && ~isempty(fig.Name)
                figName = string(fig.Name);
            else
                figName = "Figure_" + string(i);
            end
        catch
            figName = "Figure_" + string(i);
        end

        cleanName = regexprep(figName, "[^\w\d\-]+", "_");
        cleanName = regexprep(cleanName, "_+", "_");
        cleanName = strip(cleanName, "_");

        if strlength(cleanName) == 0
            cleanName = "Figure_" + string(i);
        end

        baseName = sprintf("%s_%02d_%s", prefix, i, cleanName);
        basePath = fullfile(outDir, baseName);

        okAny = false;

        if any(formats == "png")
            try
                print(fig, char(basePath + ".png"), "-dpng", "-r300");
                okAny = true;
            catch ME
                warning("PNG kaydi basarisiz: %s\n%s", baseName, ME.message);
            end
        end

        if any(formats == "pdf")
            try
                set(fig, "Renderer", "opengl");
                print(fig, char(basePath + ".pdf"), "-dpdf", "-bestfit");
                okAny = true;
            catch ME
                warning("PDF kaydi basarisiz: %s\n%s", baseName, ME.message);
            end
        end

        if any(formats == "fig")
            try
                savefig(fig, char(basePath + ".fig"));
                okAny = true;
            catch ME
                warning("FIG kaydi basarisiz: %s\n%s", baseName, ME.message);
            end
        end

        if okAny
            savedCount = savedCount + 1;
            fprintf("Saved: %s\n", baseName);
        end
    end

    fprintf("=== Done. Saved %d figures. ===\n\n", savedCount);
end
