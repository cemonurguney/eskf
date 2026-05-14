clear; clc; close all;

%% ============================================================
% Fixed-wing files
% ============================================================

fixedwing_file = "0002_20.48_fixedwing_eskf.mat";
combined_csv_file = "0002_20.48.csv";

%% ---------------- USER FLAGS ----------------
USE_BARO = true;
USE_AIRSPEED = true;

USE_ATT_INIT_FOR_DEBUG = true;
USE_COMBINED_PX4_REFERENCE = true;

SAVE_RUN_OUTPUT = true;
save_file = "fixedwing_run_gps_dropout_observability.mat";

%% ---------------- GPS MEASUREMENT SCHEDULE ----------------
% "always"  : GPS position + GPS velocity tüm uçuş boyunca açık
% "never"   : GPS position + GPS velocity tüm uçuş boyunca kapalı
% "windows" : sadece verilen zaman pencerelerinde açık
GPS_MEAS_MODE = "windows";

% Resetlenmiş sim zamanı üzerinden saniye cinsinden.
% Örnek: 120-150 s arası GPS kesintisi istiyorsan:
% GPS_MEAS_WINDOWS = [0 120; 150 375];
%
% Örnek: İlk 30 s GPS var, 180-220 s tekrar var:
% GPS_MEAS_WINDOWS = [0 30; 180 220];

GPS_MEAS_WINDOWS = [
    0    120
    150  375
];

%% ---------------- BARO BIAS OBSERVABILITY ----------------
% Baro bias sadece GPS position height yakın zamanda geldiyse estimate edilir.
% GPS yokken b_baro freeze kalır.
ESTIMATE_BARO_BIAS_WITH_GPS = true;

% GPS 5 Hz civarı. 0.5 s, yaklaşık 2-3 GPS position sample aralığı.
GPS_HEIGHT_ANCHOR_TIMEOUT_S = 0.50;

% GPS gidince bias state'i freeze edilir.
BARO_BIAS_FREEZE_STD = 0.05;      % [m]

% GPS geri gelince bias tekrar oynayabilsin diye covariance açılır.
BARO_BIAS_REACQUIRE_STD = 1.0;    % [m]

%% ---------------- GPS OUTAGE COVARIANCE HANDLING ----------------
GPS_OUTAGE_POS_INFLATE_NE = [6.0; 12.0];      % [m]
GPS_OUTAGE_VEL_INFLATE_NE = [1.0; 1.2];       % [m/s]

GPS_OUTAGE_SIGMA_A_SCALE     = 3.0;
GPS_OUTAGE_SIGMA_G_SCALE     = 1.5;
GPS_OUTAGE_SIGMA_BA_RW_SCALE = 4.0;

%% ============================================================
% 1) Init filter
% ============================================================

[state, ~, params] = init_filter();

%% ============================================================
% 1.5) Fixed-wing 18-state tuning
% ============================================================

params.sensor_profile = 'PX4 Fixed-Wing ULog: IMU + GPS + baro + TAS';

% NED frame
params.g_n = [0; 0; 9.81];

%% IMU noise
params.sigma_g = deg2rad(0.20);
params.sigma_a = 0.080;

params.sigma_bg_rw = deg2rad(0.010);
params.sigma_ba_rw = 0.010;

%% GPS noise
params.sigma_gps_pos = [3.0; 3.0; 5.0];
params.R_gps_pos = diag(params.sigma_gps_pos.^2);

params.sigma_gps_vel = [0.50; 0.50; 0.80];
params.R_gps_vel = diag(params.sigma_gps_vel.^2);

%% Barometer offset tuning
params.sigma_baro = 3.0;
params.R_baro = params.sigma_baro^2;

params.sigma_baro_bias0 = 10.0;
params.sigma_baro_bias_rw = 0.01;

params.max_baro_update_rate_hz = 20;

% Loop içinde dinamik güncellenecek.
params.estimate_baro_bias = false;

%% Wind / airspeed tuning
params.sigma_tas = 2.5;
params.R_tas = params.sigma_tas^2;

params.sigma_wind0 = [8.0; 8.0];
params.sigma_wind_rw = 0.02;

params.max_airspeed_update_rate_hz = 30;

params.use_joseph_form = true;

%% ============================================================
% 2) Build sim from fixed-wing MAT
% ============================================================

t_start = 600;
t_end   = 975;

sim = build_sim_from_fixedwing_mat(fixedwing_file, t_start, t_end);

t = sim.t;
N = numel(t);

fprintf("\n=== MAIN FIXED-WING 18-STATE ESKF, GPS DROPOUT OBSERVABILITY TEST ===\n");
fprintf("N = %d samples\n", N);
fprintf("t range = %.3f to %.3f s\n", t(1), t(end));

fprintf('\n--- GPS MEASUREMENT WINDOWS ---\n');
disp(GPS_MEAS_WINDOWS);

tmp_use_gps = false(1,N);
for kk = 1:N
    tmp_use_gps(kk) = is_gps_schedule_active(sim.t(kk), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);
end

fprintf('Expected GPS active time [s]   : %.3f\n', sum(tmp_use_gps) * median(diff(sim.t), "omitnan"));
fprintf('Expected GPS active samples    : %d / %d\n', sum(tmp_use_gps), N);

%% ============================================================
% 2.5) Barometer bias initial value
% ============================================================

% Bu deneyde başta kalibrasyon yok.
% b_baro = 0 başlıyor, GPS measurement pencerelerinde öğreniliyor,
% GPS gidince freeze ediliyor.
b_baro0 = 0;

fprintf('[fixedwing main] Dynamic GPS dropout test: b_baro starts at zero.\n');
fprintf('[fixedwing main] b_baro will be estimated only after recent GPS position height updates.\n');

%% ============================================================
% 3) Initial nominal state
% ============================================================

idx_gps0 = find(sim.gps_pos_available & all(isfinite(sim.gps_pos),1), 1, 'first');

if isempty(idx_gps0)
    warning('[fixedwing main] GPS ölçümü yok, p0 = [0;0;0] kullanılıyor.');
    state.p_n = [0;0;0];
else
    state.p_n = sim.gps_pos(:,idx_gps0);
end

idx_gps_vel0 = find(sim.gps_vel_available & all(isfinite(sim.gps_vel),1), 1, 'first');

if ~isempty(idx_gps_vel0)
    state.v_n = sim.gps_vel(:,idx_gps_vel0);
    fprintf('[fixedwing main] Initial velocity GPS velocity üzerinden başlatıldı.\n');
else
    state.v_n = [0;0;0];
    fprintf('[fixedwing main] Initial velocity v0 = [0;0;0] seçildi.\n');
end

if USE_ATT_INIT_FOR_DEBUG && ~isempty(sim.q_ref)
    state.q_nb = sim.q_ref(:,1);
    state.q_nb = state.q_nb / norm(state.q_nb);
    fprintf('[fixedwing main] DEBUG: attitude q0 log attitude üzerinden başlatıldı.\n');
else
    state.q_nb = [1;0;0;0];
    fprintf('[fixedwing main] REAL-LIKE: attitude q0 identity seçildi.\n');
end

state.b_g = [0;0;0];
state.b_a = [0;0;0];

% 16th nominal state
state.b_baro = b_baro0;

% 17-18th nominal states
state.wind_ne = [0;0];

fprintf('[fixedwing main] Initial p = [%.3f %.3f %.3f]^T m\n', state.p_n);
fprintf('[fixedwing main] Initial v = [%.3f %.3f %.3f]^T m/s\n', state.v_n);
fprintf('[fixedwing main] Initial bg = [%.6f %.6f %.6f]^T rad/s\n', state.b_g);
fprintf('[fixedwing main] Initial ba = [%.6f %.6f %.6f]^T m/s^2\n', state.b_a);
fprintf('[fixedwing main] Initial b_baro = %.6f m\n', state.b_baro);
fprintf('[fixedwing main] Initial wind NE = [%.6f %.6f]^T m/s\n', state.wind_ne);
fprintf('[fixedwing main] USE_BARO = %d\n', USE_BARO);
fprintf('[fixedwing main] USE_AIRSPEED = %d\n', USE_AIRSPEED);
fprintf('[fixedwing main] GPS_MEAS_MODE = %s\n', GPS_MEAS_MODE);
fprintf('[fixedwing main] ESTIMATE_BARO_BIAS_WITH_GPS = %d\n', ESTIMATE_BARO_BIAS_WITH_GPS);

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

if USE_ATT_INIT_FOR_DEBUG
    sigma_th0 = deg2rad([3.0; 3.0; 8.0]);
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
% 5) Allocate logs
% ============================================================

log_p = nan(3,N);
log_v = nan(3,N);
log_q = nan(4,N);
log_bg = nan(3,N);
log_ba = nan(3,N);
log_bbaro = nan(1,N);
log_wind_ne = nan(2,N);
log_Pdiag = nan(18,N);

log_res_gps_pos = nan(3,N);
log_res_gps_vel = nan(3,N);
log_res_baro = nan(1,N);
log_res_airspeed = nan(1,N);

log_use_gps_meas = false(1,N);
log_estimate_baro_bias = false(1,N);

used_gps_pos_updates = 0;
used_gps_vel_updates = 0;
used_baro_updates = 0;
used_airspeed_updates = 0;

log_p(:,1) = state.p_n;
log_v(:,1) = state.v_n;
log_q(:,1) = state.q_nb;
log_bg(:,1) = state.b_g;
log_ba(:,1) = state.b_a;
log_bbaro(1) = state.b_baro;
log_wind_ne(:,1) = state.wind_ne;
log_Pdiag(:,1) = diag(P);

%% ============================================================
% 6) Update rate limiters
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

%% ============================================================
% 6.5) Dynamic GPS / baro-bias observability state
% ============================================================

prev_use_gps_meas = is_gps_schedule_active(t(1), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);

last_gps_pos_update_t = -inf;
last_gps_vel_update_t = -inf;

baro_bias_hold_value = state.b_baro;

log_use_gps_meas(1) = prev_use_gps_meas;
log_estimate_baro_bias(1) = false;

%% ============================================================
% 7) Main ESKF replay loop
% ============================================================

fprintf("\n[fixedwing main] ESKF replay başlıyor...\n");

for k = 2:N

    dt_k = sim.t(k) - sim.t(k-1);

    if dt_k <= 0 || dt_k > 0.1
        warning("Bad dt at k=%d: %.6f s, propagation skipped.", k, dt_k);
        continue;
    end

    %% ---------------- GPS measurement schedule ----------------
    use_gps_meas_now = is_gps_schedule_active(sim.t(k), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);

    % GPS ON -> OFF:
    % b_baro artık gözlenebilir değil, o andaki değeri freeze ediyoruz.
    % Aynı anda yatay P şişiriyoruz, çünkü GPS yokken N/E drift beklenir.
    if prev_use_gps_meas && ~use_gps_meas_now
        baro_bias_hold_value = state.b_baro;

        P(16,:) = 0;
        P(:,16) = 0;
        P(16,16) = BARO_BIAS_FREEZE_STD^2;

        % GPS outage covariance inflation
        P(1:2,1:2) = P(1:2,1:2) + diag(GPS_OUTAGE_POS_INFLATE_NE.^2);
        P(4:5,4:5) = P(4:5,4:5) + diag(GPS_OUTAGE_VEL_INFLATE_NE.^2);
        P = 0.5 * (P + P.');

        fprintf('[fixedwing main] t=%.2f s: GPS OFF, freezing b_baro = %.6f m\n', ...
            sim.t(k), baro_bias_hold_value);
    end

    % GPS OFF -> ON:
    % Bias tekrar gözlenebilir olacak. Covariance'ı biraz açıyoruz.
    if ~prev_use_gps_meas && use_gps_meas_now
        P(16,:) = 0;
        P(:,16) = 0;
        P(16,16) = BARO_BIAS_REACQUIRE_STD^2;

        fprintf('[fixedwing main] t=%.2f s: GPS ON, re-enabling b_baro estimation\n', ...
            sim.t(k));
    end

    % GPS yokken bias gerçekten sabit kalsın.
    if ~use_gps_meas_now
        state.b_baro = baro_bias_hold_value;
    end

    %% ---------------- Per-step params ----------------
    params_k = params;

    % GPS yokken process noise artırılır.
    % Ama baro bias random walk kapatılır, çünkü b_baro freeze.
    if ~use_gps_meas_now
        params_k.sigma_baro_bias_rw = 0.0;

        params_k.sigma_a     = params.sigma_a     * GPS_OUTAGE_SIGMA_A_SCALE;
        params_k.sigma_g     = params.sigma_g     * GPS_OUTAGE_SIGMA_G_SCALE;
        params_k.sigma_ba_rw = params.sigma_ba_rw * GPS_OUTAGE_SIGMA_BA_RW_SCALE;
    end

    params_k.estimate_baro_bias = false;

    %% ---------------- IMU propagation ----------------
    imu.gyro_m  = sim.imu_gyro(:,k);
    imu.accel_m = sim.imu_accel(:,k);

    state = propagate_nominal(state, imu, params_k, dt_k);

    [F, G, Qd] = compute_F_G_Qd(state, imu, params_k, dt_k);
    P = propagate_covariance(P, F, G, Qd, dt_k);

    %% ---------------- GPS position update ----------------
    if use_gps_meas_now && sim.gps_pos_available(k)
        z_gps_pos = sim.gps_pos(:,k);

        if all(isfinite(z_gps_pos))
            [state, P, residual, ~, ~] = update_gnss_pos(state, P, z_gps_pos, params_k);
            log_res_gps_pos(:,k) = residual;
            used_gps_pos_updates = used_gps_pos_updates + 1;
            last_gps_pos_update_t = sim.t(k);
        end
    end

    %% ---------------- GPS velocity update ----------------
    if use_gps_meas_now && sim.gps_vel_available(k)
        z_gps_vel = sim.gps_vel(:,k);

        if all(isfinite(z_gps_vel))
            [state, P, residual, ~, ~] = update_gnss_vel(state, P, z_gps_vel, params_k);
            log_res_gps_vel(:,k) = residual;
            used_gps_vel_updates = used_gps_vel_updates + 1;
            last_gps_vel_update_t = sim.t(k);
        end
    end

    %% ---------------- Is baro bias observable now? ----------------
    % Hz farkı burada önemli:
    % GPS 5 Hz, baro ~20 Hz kullanılıyor.
    % Baro bias'ı sadece son GPS position height update'i yakın zamandaysa açıyoruz.
    gps_height_recent = use_gps_meas_now && ...
        ((sim.t(k) - last_gps_pos_update_t) <= GPS_HEIGHT_ANCHOR_TIMEOUT_S);

    estimate_baro_bias_now = ESTIMATE_BARO_BIAS_WITH_GPS && gps_height_recent;

    params_k.estimate_baro_bias = estimate_baro_bias_now;

    if ~estimate_baro_bias_now
        params_k.sigma_baro_bias_rw = 0.0;
    end

    %% ---------------- Baro update ----------------
    if USE_BARO && sim.baro_available(k)

        if sim.t(k) - last_baro_update_t >= baro_min_dt

            z_baro = sim.baro(k);

            if isfinite(z_baro)
                [state, P, residual, ~, ~] = update_baro(state, P, z_baro, params_k);

                log_res_baro(k) = residual;
                used_baro_updates = used_baro_updates + 1;
                last_baro_update_t = sim.t(k);
            end
        end
    end

    %% ---------------- Airspeed update, 18-state wind model ----------------
    if USE_AIRSPEED && isfield(sim, "airspeed_available") && sim.airspeed_available(k)

        if sim.t(k) - last_airspeed_update_t >= airspeed_min_dt

            z_tas = sim.airspeed(k);

            if isfinite(z_tas) && z_tas > 3
                [state, P, residual, ~, ~] = update_airspeed(state, P, z_tas, params_k);

                if isfinite(residual)
                    log_res_airspeed(k) = residual;
                    used_airspeed_updates = used_airspeed_updates + 1;
                    last_airspeed_update_t = sim.t(k);
                end
            end
        end
    end

    %% ---------------- Enforce b_baro freeze during GPS outage ----------------
    if ~use_gps_meas_now
        state.b_baro = baro_bias_hold_value;

        P(16,:) = 0;
        P(:,16) = 0;
        P(16,16) = BARO_BIAS_FREEZE_STD^2;
    end

    %% ---------------- Logging ----------------
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

fprintf("[fixedwing main] ESKF replay bitti.\n");

%% ============================================================
% 8) Plot observation/reference results
% ============================================================

plot_results_fixedwing_obs(sim, t, log_p, log_v, log_q, ...
    log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params);

%% ============================================================
% 9) Metrics
% ============================================================

idx_gps_pos_all = sim.gps_pos_available & all(isfinite(sim.gps_pos),1);
idx_gps_vel_all = sim.gps_vel_available & all(isfinite(sim.gps_vel),1);
idx_baro = sim.baro_available & isfinite(sim.baro);
idx_air = isfinite(log_res_airspeed);

idx_gps_pos_used = idx_gps_pos_all & log_use_gps_meas;
idx_gps_pos_holdout = idx_gps_pos_all & ~log_use_gps_meas;

idx_gps_vel_used = idx_gps_vel_all & log_use_gps_meas;
idx_gps_vel_holdout = idx_gps_vel_all & ~log_use_gps_meas;

gps_pos_rmse = nan(3,1);
gps_pos_rmse_norm = nan;
gps_vel_rmse = nan(3,1);
gps_vel_rmse_norm = nan;

gps_pos_rmse_used = nan(3,1);
gps_pos_rmse_used_norm = nan;
gps_pos_rmse_holdout = nan(3,1);
gps_pos_rmse_holdout_norm = nan;

gps_vel_rmse_used = nan(3,1);
gps_vel_rmse_used_norm = nan;
gps_vel_rmse_holdout = nan(3,1);
gps_vel_rmse_holdout_norm = nan;

baro_rmse = nan;
baro_mean_error = nan;
baro_std_error = nan;

baro_raw_rmse = nan;
baro_raw_mean_error = nan;
baro_raw_std_error = nan;

baro_model_rmse = nan;
baro_model_mean_error = nan;
baro_model_std_error = nan;

airspeed_rmse = nan;
airspeed_mean_error = nan;
airspeed_std_error = nan;

fprintf('\n--- FIXED-WING RUN FINISHED ---\n');

if any(idx_gps_pos_all)
    gps_pos_err = log_p(:,idx_gps_pos_all) - sim.gps_pos(:,idx_gps_pos_all);
    gps_pos_rmse = sqrt(mean(gps_pos_err.^2, 2, 'omitnan'));
    gps_pos_rmse_norm = sqrt(mean(sum(gps_pos_err.^2,1), 'omitnan'));

    fprintf('GPS position ALL RMSE N/E/D [m]      : [%.3f %.3f %.3f]\n', gps_pos_rmse);
    fprintf('GPS position ALL RMSE norm [m]       : %.3f\n', gps_pos_rmse_norm);
end

if any(idx_gps_pos_used)
    err = log_p(:,idx_gps_pos_used) - sim.gps_pos(:,idx_gps_pos_used);
    gps_pos_rmse_used = sqrt(mean(err.^2, 2, 'omitnan'));
    gps_pos_rmse_used_norm = sqrt(mean(sum(err.^2,1), 'omitnan'));

    fprintf('GPS position USED RMSE N/E/D [m]     : [%.3f %.3f %.3f]\n', gps_pos_rmse_used);
    fprintf('GPS position USED RMSE norm [m]      : %.3f\n', gps_pos_rmse_used_norm);
end

if any(idx_gps_pos_holdout)
    err = log_p(:,idx_gps_pos_holdout) - sim.gps_pos(:,idx_gps_pos_holdout);
    gps_pos_rmse_holdout = sqrt(mean(err.^2, 2, 'omitnan'));
    gps_pos_rmse_holdout_norm = sqrt(mean(sum(err.^2,1), 'omitnan'));

    fprintf('GPS position HOLDOUT RMSE N/E/D [m]  : [%.3f %.3f %.3f]\n', gps_pos_rmse_holdout);
    fprintf('GPS position HOLDOUT RMSE norm [m]   : %.3f\n', gps_pos_rmse_holdout_norm);
end

if any(idx_gps_vel_all)
    gps_vel_err = log_v(:,idx_gps_vel_all) - sim.gps_vel(:,idx_gps_vel_all);
    gps_vel_rmse = sqrt(mean(gps_vel_err.^2, 2, 'omitnan'));
    gps_vel_rmse_norm = sqrt(mean(sum(gps_vel_err.^2,1), 'omitnan'));

    fprintf('GPS velocity ALL RMSE N/E/D [m/s]      : [%.3f %.3f %.3f]\n', gps_vel_rmse);
    fprintf('GPS velocity ALL RMSE norm [m/s]       : %.3f\n', gps_vel_rmse_norm);
end

if any(idx_gps_vel_used)
    err = log_v(:,idx_gps_vel_used) - sim.gps_vel(:,idx_gps_vel_used);
    gps_vel_rmse_used = sqrt(mean(err.^2, 2, 'omitnan'));
    gps_vel_rmse_used_norm = sqrt(mean(sum(err.^2,1), 'omitnan'));

    fprintf('GPS velocity USED RMSE N/E/D [m/s]     : [%.3f %.3f %.3f]\n', gps_vel_rmse_used);
    fprintf('GPS velocity USED RMSE norm [m/s]      : %.3f\n', gps_vel_rmse_used_norm);
end

if any(idx_gps_vel_holdout)
    err = log_v(:,idx_gps_vel_holdout) - sim.gps_vel(:,idx_gps_vel_holdout);
    gps_vel_rmse_holdout = sqrt(mean(err.^2, 2, 'omitnan'));
    gps_vel_rmse_holdout_norm = sqrt(mean(sum(err.^2,1), 'omitnan'));

    fprintf('GPS velocity HOLDOUT RMSE N/E/D [m/s]  : [%.3f %.3f %.3f]\n', gps_vel_rmse_holdout);
    fprintf('GPS velocity HOLDOUT RMSE norm [m/s]   : %.3f\n', gps_vel_rmse_holdout_norm);
end

if any(idx_baro)
    baro_raw_err = log_p(3,idx_baro) - sim.baro(idx_baro);
    baro_model_err = log_p(3,idx_baro) + log_bbaro(idx_baro) - sim.baro(idx_baro);

    baro_raw_rmse = sqrt(mean(baro_raw_err.^2, 'omitnan'));
    baro_raw_mean_error = mean(baro_raw_err, 'omitnan');
    baro_raw_std_error = std(baro_raw_err, 0, 'omitnan');

    baro_model_rmse = sqrt(mean(baro_model_err.^2, 'omitnan'));
    baro_model_mean_error = mean(baro_model_err, 'omitnan');
    baro_model_std_error = std(baro_model_err, 0, 'omitnan');

    baro_rmse = baro_model_rmse;
    baro_mean_error = baro_model_mean_error;
    baro_std_error = baro_model_std_error;

    fprintf('Baro raw Down RMSE [m]              : %.3f\n', baro_raw_rmse);
    fprintf('Baro raw Down mean error [m]        : %.3f\n', baro_raw_mean_error);
    fprintf('Baro raw Down std error [m]         : %.3f\n', baro_raw_std_error);

    fprintf('Baro modeled Down RMSE [m]          : %.3f\n', baro_model_rmse);
    fprintf('Baro modeled Down mean error [m]    : %.3f\n', baro_model_mean_error);
    fprintf('Baro modeled Down std error [m]     : %.3f\n', baro_model_std_error);
end

if any(idx_air)
    airspeed_rmse = sqrt(mean(log_res_airspeed(idx_air).^2, 'omitnan'));
    airspeed_mean_error = mean(log_res_airspeed(idx_air), 'omitnan');
    airspeed_std_error = std(log_res_airspeed(idx_air), 0, 'omitnan');

    fprintf('Airspeed residual RMSE [m/s]        : %.3f\n', airspeed_rmse);
    fprintf('Airspeed residual mean [m/s]        : %.3f\n', airspeed_mean_error);
    fprintf('Airspeed residual std [m/s]         : %.3f\n', airspeed_std_error);
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

        fprintf('Wind vs PX4 ref RMSE N/E [m/s]     : [%.3f %.3f]\n', wind_ref_rmse);
        fprintf('Wind vs PX4 ref RMSE norm [m/s]    : %.3f\n', wind_ref_rmse_norm);
        fprintf('Wind vs PX4 ref mean N/E [m/s]     : [%.3f %.3f]\n', wind_ref_mean_error);
    end
end

fprintf('\n--- UPDATE COUNTS ---\n');
fprintf('GPS pos updates available : %d\n', sum(sim.gps_pos_available));
fprintf('GPS pos updates used      : %d\n', used_gps_pos_updates);
fprintf('GPS vel updates available : %d\n', sum(sim.gps_vel_available));
fprintf('GPS vel updates used      : %d\n', used_gps_vel_updates);
fprintf('BARO updates available    : %d\n', sum(sim.baro_available & isfinite(sim.baro)));
fprintf('BARO updates used         : %d\n', used_baro_updates);

if isfield(sim, "airspeed_available")
    fprintf('Airspeed updates available : %d\n', sum(sim.airspeed_available));
end
fprintf('Airspeed updates used      : %d\n', used_airspeed_updates);

fprintf('GPS measurement scheduled samples       : %d\n', sum(log_use_gps_meas));
fprintf('Baro bias estimated samples             : %d\n', sum(log_estimate_baro_bias));

fprintf('\n--- FINAL STATES ---\n');
fprintf('Final estimated p NED [m]     : [%.4f %.4f %.4f]^T\n', log_p(:,end));
fprintf('Final estimated v NED [m/s]   : [%.4f %.4f %.4f]^T\n', log_v(:,end));
fprintf('Final gyro bias [rad/s]       : [%.6f %.6f %.6f]^T\n', log_bg(:,end));
fprintf('Final accel bias [m/s^2]      : [%.6f %.6f %.6f]^T\n', log_ba(:,end));
fprintf('Final baro offset [m]         : %.6f\n', log_bbaro(end));
fprintf('Final wind NE [m/s]           : [%.6f %.6f]^T\n', log_wind_ne(:,end));

fprintf('\n--- CONFIG USED ---\n');
fprintf('Sensor profile : %s\n', params.sensor_profile);
fprintf('USE_BARO       : %d\n', USE_BARO);
fprintf('USE_AIRSPEED   : %d\n', USE_AIRSPEED);
fprintf('ATT DEBUG INIT : %d\n', USE_ATT_INIT_FOR_DEBUG);
fprintf('GPS_MEAS_MODE  : %s\n', GPS_MEAS_MODE);
fprintf('GPS_HEIGHT_ANCHOR_TIMEOUT_S : %.3f\n', GPS_HEIGHT_ANCHOR_TIMEOUT_S);
fprintf('EST BARO BIAS WITH GPS      : %d\n', ESTIMATE_BARO_BIAS_WITH_GPS);
fprintf('GPS_OUTAGE_POS_INFLATE_NE   : [%.3f %.3f]\n', GPS_OUTAGE_POS_INFLATE_NE);
fprintf('GPS_OUTAGE_VEL_INFLATE_NE   : [%.3f %.3f]\n', GPS_OUTAGE_VEL_INFLATE_NE);
fprintf('GPS_OUTAGE_SIGMA_A_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_A_SCALE);
fprintf('GPS_OUTAGE_SIGMA_G_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_G_SCALE);
fprintf('GPS_OUTAGE_SIGMA_BA_RW_SCALE : %.3f\n', GPS_OUTAGE_SIGMA_BA_RW_SCALE);
fprintf('sigma_gps_pos  : [%.3f %.3f %.3f]\n', params.sigma_gps_pos);
fprintf('sigma_gps_vel  : [%.3f %.3f %.3f]\n', params.sigma_gps_vel);
fprintf('sigma_baro     : %.3f\n', params.sigma_baro);
fprintf('sigma_tas      : %.3f\n', params.sigma_tas);
fprintf('sigma_wind0    : [%.3f %.3f]\n', params.sigma_wind0);
fprintf('sigma_wind_rw  : %.6f\n', params.sigma_wind_rw);

%% ============================================================
% 10) PX4 reference comparison from combined CSV
% ============================================================

if USE_COMBINED_PX4_REFERENCE
    px4_cmp = compare_with_px4_reference_from_csv_autoalign( ...
        combined_csv_file, t, log_p, log_v, sim);
else
    px4_cmp = struct();
    px4_cmp.available = false;
end

%% ============================================================
% 11) Extra plots
% ============================================================

figure('Name','Estimated Barometer Offset');
plot(t, log_bbaro, 'LineWidth', 1.3);
grid on;
xlabel('Time [s]');
ylabel('b_{baro} [m]');
title('Estimated Barometer Offset State');

figure('Name','GPS Measurement Schedule and Baro Bias Observability');
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
title('GPS Dropout Schedule and Barometer Bias Observability');

figure('Name','Estimated Wind NE');
plot(t, log_wind_ne(1,:), 'LineWidth', 1.3);
hold on;
plot(t, log_wind_ne(2,:), 'LineWidth', 1.3);

if isfield(sim, "wind_ref_available") && any(sim.wind_ref_available)
    plot(t(sim.wind_ref_available), sim.wind_ref(1,sim.wind_ref_available), '.');
    plot(t(sim.wind_ref_available), sim.wind_ref(2,sim.wind_ref_available), '.');
    legend('ESKF wind_N','ESKF wind_E','PX4 wind_N ref','PX4 wind_E ref','Location','best');
else
    legend('ESKF wind_N','ESKF wind_E','Location','best');
end

grid on;
xlabel('Time [s]');
ylabel('Wind [m/s]');
title('Estimated Horizontal Wind');

figure('Name','Airspeed Residual');
plot(t(idx_air), log_res_airspeed(idx_air), '.');
grid on;
xlabel('Time [s]');
ylabel('TAS residual [m/s]');
title('Airspeed Innovation: z_{TAS} - ||v - w||');

figure('Name','GPS Position Holdout Error');
if any(idx_gps_pos_holdout)
    err_hold = log_p(:,idx_gps_pos_holdout) - sim.gps_pos(:,idx_gps_pos_holdout);
    th = t(idx_gps_pos_holdout);

    plot(th, err_hold(1,:), '.'); hold on;
    plot(th, err_hold(2,:), '.');
    plot(th, err_hold(3,:), '.');
    grid on;
    xlabel('Time [s]');
    ylabel('Position error wrt held-out GPS [m]');
    legend('N','E','D');
    title('GPS Position Holdout Error During GPS Dropout');
end

%% ============================================================
% 12) Save run output
% ============================================================

if SAVE_RUN_OUTPUT
    fprintf("\n[fixedwing main] Saving run output to %s ...\n", save_file);

    save(save_file, ...
        "fixedwing_file", "combined_csv_file", ...
        "sim", "t", ...
        "log_p", "log_v", "log_q", ...
        "log_bg", "log_ba", "log_bbaro", "log_wind_ne", "log_Pdiag", ...
        "log_res_gps_pos", "log_res_gps_vel", "log_res_baro", "log_res_airspeed", ...
        "log_use_gps_meas", "log_estimate_baro_bias", ...
        "gps_pos_rmse", "gps_pos_rmse_norm", ...
        "gps_vel_rmse", "gps_vel_rmse_norm", ...
        "gps_pos_rmse_used", "gps_pos_rmse_used_norm", ...
        "gps_pos_rmse_holdout", "gps_pos_rmse_holdout_norm", ...
        "gps_vel_rmse_used", "gps_vel_rmse_used_norm", ...
        "gps_vel_rmse_holdout", "gps_vel_rmse_holdout_norm", ...
        "baro_rmse", "baro_mean_error", "baro_std_error", ...
        "baro_raw_rmse", "baro_raw_mean_error", "baro_raw_std_error", ...
        "baro_model_rmse", "baro_model_mean_error", "baro_model_std_error", ...
        "airspeed_rmse", "airspeed_mean_error", "airspeed_std_error", ...
        "wind_ref_rmse", "wind_ref_rmse_norm", "wind_ref_mean_error", ...
        "used_gps_pos_updates", "used_gps_vel_updates", ...
        "used_baro_updates", "used_airspeed_updates", ...
        "USE_BARO", "USE_AIRSPEED", "USE_ATT_INIT_FOR_DEBUG", ...
        "GPS_MEAS_MODE", "GPS_MEAS_WINDOWS", ...
        "GPS_HEIGHT_ANCHOR_TIMEOUT_S", ...
        "ESTIMATE_BARO_BIAS_WITH_GPS", ...
        "BARO_BIAS_FREEZE_STD", "BARO_BIAS_REACQUIRE_STD", ...
        "GPS_OUTAGE_POS_INFLATE_NE", "GPS_OUTAGE_VEL_INFLATE_NE", ...
        "GPS_OUTAGE_SIGMA_A_SCALE", "GPS_OUTAGE_SIGMA_G_SCALE", ...
        "GPS_OUTAGE_SIGMA_BA_RW_SCALE", ...
        "USE_COMBINED_PX4_REFERENCE", ...
        "params", "px4_cmp", "-v7.3");

    fprintf("[fixedwing main] Saved: %s\n", save_file);
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