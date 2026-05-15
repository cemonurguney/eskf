clear; clc; close all;

%% ============================================================
% RAW MRS / QUAD DATA MODE - 18 STATE ESKF
% ============================================================
% Inputs:
% - raw IMU
% - raw RTK/GPS position
% - pseudo GPS velocity
% - raw barometer
%
% Reference / truth:
% - sim.p_true, sim.v_true, sim.q_true only for plotting/error metrics.
%
% Important:
% Truth is NOT used as measurement update.
%
% State:
%   nominal:
%     p_n, v_n, q_nb, b_g, b_a, b_baro, wind_ne
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%      = 18x1
%
% For quad/raw data:
%   - wind state exists for common 18-state compatibility.
%   - no TAS update is used, so wind is not physically interpreted here.
% ============================================================

real_data_file = "rectangle_raw_clean.mat";

%% ---------------- USER FLAGS ----------------
USE_BARO = true;
USE_PSEUDO_GPS_VEL = true;

% Debug attitude initialization:
% true  -> q0 truth/reference üzerinden başlar. Measurement değildir.
% false -> q0 identity; daha gerçekçi ama yaw bilinmezse bozulabilir.
USE_TRUTH_ATT_INIT_FOR_DEBUG = true;

% Initial velocity:
% true  -> ilk pseudo GPS velocity varsa v0 olarak kullanır.
% false -> v0 = [0;0;0]
USE_PSEUDO_GPS_VEL_INIT = true;

SAVE_RUN_OUTPUT = true;
save_file = "raw_quad_run_gps_dropout_observability.mat";

%% ---------------- GPS MEASUREMENT SCHEDULE ----------------
% "always"  : GPS position + GPS velocity tüm koşu boyunca açık
% "never"   : GPS position + GPS velocity tüm koşu boyunca kapalı
% "windows" : sadece verilen zaman pencerelerinde açık
GPS_MEAS_MODE = "always";

% Dropout denemesi yapmak istersen:
% GPS_MEAS_MODE = "windows";
% GPS_MEAS_WINDOWS = [
%     0    60
%     90   inf
% ];
GPS_MEAS_WINDOWS = [
    0 inf
];

%% ---------------- BARO BIAS OBSERVABILITY ----------------
% Baro bias sadece GPS position height yakın zamanda geldiyse estimate edilir.
% GPS yokken b_baro nominal state olarak son değerde freeze kalır.
ESTIMATE_BARO_BIAS_WITH_GPS = true;

% GPS 5 Hz ise 0.5 s yaklaşık 2-3 GPS sample toleransı.
GPS_HEIGHT_ANCHOR_TIMEOUT_S = 0.50;

%% ---------------- GPS OUTAGE PROCESS-NOISE HANDLING ----------------
% P'ye elle müdahale yok.
% GPS yokken covariance Qd üzerinden büyüsün diye process noise artırılır.
GPS_OUTAGE_SIGMA_A_SCALE     = 3.0;
GPS_OUTAGE_SIGMA_G_SCALE     = 1.5;
GPS_OUTAGE_SIGMA_BA_RW_SCALE = 4.0;

fprintf('[main raw] Raw gerçek veri modu: %s yükleniyor...\n', real_data_file);

%% ============================================================
% 1) Init filter
% ============================================================

[state, ~, params] = init_filter();

%% ============================================================
% 1.5) RAW / QUAD SENSOR-INFORMED TUNING
% ============================================================

params.sensor_profile = 'MRS MAS Rectangle: Pixhawk4 IMU + MS5611 baro + Emlid Reach RTK';

%% ---------------- Frame convention ----------------
% ENU convention:
% p = [E; N; U]
% gravity points downward, so z/up gravity is negative.
params.g_n = [0; 0; -9.81];

%% ---------------- IMU process noise ----------------
params.sigma_g = deg2rad(0.20);
params.sigma_a = 0.060;

params.sigma_bg_rw = deg2rad(0.010);
params.sigma_ba_rw = 0.008;

%% ---------------- Raw RTK/GPS position measurement ----------------
params.sigma_gps_pos = [0.40; 0.40; 0.60];     % [E; N; U] [m]
params.R_gps_pos = diag(params.sigma_gps_pos.^2);

%% ---------------- Pseudo GPS velocity measurement ----------------
params.sigma_gps_vel = [0.45; 0.45; 0.90];     % [vE; vN; vU] [m/s]
params.R_gps_vel = diag(params.sigma_gps_vel.^2);

%% ---------------- Barometer measurement ----------------
params.sigma_baro = 3.0;                       % [m]
params.R_baro = params.sigma_baro^2;

params.sigma_baro_bias0 = 10.0;                % [m]
params.sigma_baro_bias_rw = 0.01;              % [m/sqrt(s)]
params.max_baro_update_rate_hz = 20;

params.estimate_baro_bias = false;

%% ---------------- Wind state, common 18-state compatibility ----------------
% Quad/raw hattında TAS yok. Wind state yalnızca ortak 18-state yapısını
% bozmayalım diye duruyor. Fiziksel wind sonucu gibi yorumlanmaz.
params.sigma_wind0 = [8.0; 8.0];
params.sigma_wind_rw = 0.0;

params.use_joseph_form = true;

%% ---------------- Innovation gates ----------------
params.gps_pos_gate_chi2 = 25.0;
params.gps_vel_gate_chi2 = 25.0;
params.baro_gate_chi2    = 9.0;

%% ============================================================
% 2) Build sim from raw MAT
% ============================================================

sim = build_sim_from_raw_mat(real_data_file);

t = sim.t;
N = numel(t);

fprintf('\n=== MAIN RAW/QUAD 18-STATE ESKF ===\n');
fprintf('N = %d samples\n', N);
fprintf('t range = %.3f to %.3f s\n', t(1), t(end));

fprintf('\n--- GPS MEASUREMENT WINDOWS ---\n');
disp(GPS_MEAS_WINDOWS);

tmp_use_gps = false(1,N);
for kk = 1:N
    tmp_use_gps(kk) = is_gps_schedule_active(sim.t(kk), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);
end

fprintf('Expected GPS active time [s]   : %.3f\n', sum(tmp_use_gps) * median(diff(sim.t), "omitnan"));
fprintf('Expected GPS active samples    : %d / %d\n', sum(tmp_use_gps), N);

%% ============================================================
% 3) Start index
% ============================================================

k_start = 1;
fprintf('[main raw] Başlangıç indeksi: k_start = %d\n', k_start);

%% ============================================================
% 4) Trim sim if needed
% ============================================================

if k_start > 1
    t_offset = sim.t(k_start);

    sim.t = sim.t(k_start:end) - t_offset;
    sim.imu_gyro = sim.imu_gyro(:,k_start:end);
    sim.imu_accel = sim.imu_accel(:,k_start:end);

    sim.p_true = sim.p_true(:,k_start:end);
    sim.v_true = sim.v_true(:,k_start:end);
    sim.q_true = sim.q_true(:,k_start:end);

    sim.gps_pos = sim.gps_pos(:,k_start:end);
    sim.gps_pos_available = sim.gps_pos_available(k_start:end);

    sim.gps_vel = sim.gps_vel(:,k_start:end);
    sim.gps_vel_available = sim.gps_vel_available(k_start:end);

    sim.baro = sim.baro(k_start:end);
    sim.baro_available = sim.baro_available(k_start:end);

    if isfield(sim,'range')
        sim.range = sim.range(k_start:end);
        sim.range_available = sim.range_available(k_start:end);
    end
end

t = sim.t;
N = numel(t);

%% ============================================================
% 5) Initial nominal state
% ============================================================

%% ---------------- Position init ----------------
idx_gps0 = find(sim.gps_pos_available & all(isfinite(sim.gps_pos),1), 1, 'first');

if isempty(idx_gps0)
    warning('[main raw] GPS ölçümü yok, p0 = [0;0;0] kullanılıyor.');
    state.p_n = [0;0;0];
else
    state.p_n = sim.gps_pos(:,idx_gps0);
end

%% ---------------- Velocity init ----------------
idx_gps_vel0 = find(sim.gps_vel_available & all(isfinite(sim.gps_vel),1), 1, 'first');

if USE_PSEUDO_GPS_VEL_INIT && ~isempty(idx_gps_vel0)
    state.v_n = sim.gps_vel(:,idx_gps_vel0);
    fprintf('[main raw] Initial velocity pseudo GPS velocity üzerinden başlatıldı.\n');
else
    state.v_n = [0;0;0];
    fprintf('[main raw] Initial velocity v0 = [0;0;0] seçildi.\n');
end

%% ---------------- Attitude init ----------------
if USE_TRUTH_ATT_INIT_FOR_DEBUG
    state.q_nb = sim.q_true(:,1);
    state.q_nb = state.q_nb / norm(state.q_nb);
    fprintf('[main raw] DEBUG: attitude q0 truth/reference üzerinden başlatıldı.\n');
else
    state.q_nb = [1;0;0;0];
    fprintf('[main raw] REAL-LIKE: attitude q0 identity seçildi. Yaw sapması beklenebilir.\n');
end

%% ---------------- Bias init ----------------
state.b_g = [0;0;0];
state.b_a = [0;0;0];

% Baro offset: başta kalibrasyon yok. GPS height geldikçe öğrenilir.
state.b_baro = 0;

% Wind state: quad/raw hattında fiziksel yorum yok.
state.wind_ne = [0;0];

fprintf('[main raw] Bias init zero kullanıldı.\n');
fprintf('[main raw] Initial p = [%.3f %.3f %.3f]^T m\n', state.p_n);
fprintf('[main raw] Initial v = [%.3f %.3f %.3f]^T m/s\n', state.v_n);
fprintf('[main raw] Initial bg = [%.6f %.6f %.6f]^T rad/s\n', state.b_g);
fprintf('[main raw] Initial ba = [%.6f %.6f %.6f]^T m/s^2\n', state.b_a);
fprintf('[main raw] Initial b_baro = %.6f m\n', state.b_baro);
fprintf('[main raw] USE_BARO = %d\n', USE_BARO);
fprintf('[main raw] USE_PSEUDO_GPS_VEL = %d\n', USE_PSEUDO_GPS_VEL);
fprintf('[main raw] USE_TRUTH_ATT_INIT_FOR_DEBUG = %d\n', USE_TRUTH_ATT_INIT_FOR_DEBUG);
fprintf('[main raw] GPS_MEAS_MODE = %s\n', GPS_MEAS_MODE);

%% ============================================================
% 6) Initial covariance, 18-state
% ============================================================

idx_p      = 1:3;
idx_v      = 4:6;
idx_th     = 7:9;
idx_bg     = 10:12;
idx_ba     = 13:15;
idx_bbaro  = 16;
idx_wind   = 17:18;

sigma_p0  = [0.50; 0.50; 0.80];
sigma_v0  = [0.35; 0.35; 0.60];

if USE_TRUTH_ATT_INIT_FOR_DEBUG
    sigma_th0 = deg2rad([0.5; 0.5; 2.0]);
else
    sigma_th0 = deg2rad([3.0; 3.0; 30.0]);
end

sigma_bg0 = deg2rad([0.50; 0.50; 0.80]);
sigma_ba0 = [0.20; 0.20; 0.30];

P = zeros(18,18);
P(idx_p, idx_p)         = diag(sigma_p0.^2);
P(idx_v, idx_v)         = diag(sigma_v0.^2);
P(idx_th, idx_th)       = diag(sigma_th0.^2);
P(idx_bg, idx_bg)       = diag(sigma_bg0.^2);
P(idx_ba, idx_ba)       = diag(sigma_ba0.^2);
P(idx_bbaro, idx_bbaro) = params.sigma_baro_bias0^2;
P(idx_wind, idx_wind)   = diag(params.sigma_wind0.^2);

%% ============================================================
% 7) Logs
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

log_use_gps_meas = false(1,N);
log_estimate_baro_bias = false(1,N);

log_gate_gps_pos = false(1,N);
log_gate_gps_vel = false(1,N);
log_gate_baro = false(1,N);

used_gps_pos_updates = 0;
used_gps_vel_updates = 0;
used_baro_updates = 0;

accepted_gps_pos_updates = 0;
accepted_gps_vel_updates = 0;
accepted_baro_updates = 0;

rejected_gps_pos_updates = 0;
rejected_gps_vel_updates = 0;
rejected_baro_updates = 0;

log_p(:,1) = state.p_n;
log_v(:,1) = state.v_n;
log_q(:,1) = state.q_nb;
log_bg(:,1) = state.b_g;
log_ba(:,1) = state.b_a;
log_bbaro(1) = state.b_baro;
log_wind_ne(:,1) = state.wind_ne;
log_Pdiag(:,1) = diag(P);

%% ============================================================
% 8) Rate limiters and dynamic state
% ============================================================

last_baro_update_t = -inf;
if isfield(params, 'max_baro_update_rate_hz') && params.max_baro_update_rate_hz > 0
    baro_min_dt = 1 / params.max_baro_update_rate_hz;
else
    baro_min_dt = 0;
end

prev_use_gps_meas = is_gps_schedule_active(t(1), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);

last_gps_pos_update_t = -inf;
last_gps_vel_update_t = -inf;

baro_bias_hold_value = state.b_baro;

log_use_gps_meas(1) = prev_use_gps_meas;
log_estimate_baro_bias(1) = false;

%% ============================================================
% 9) Main ESKF loop
% ============================================================

fprintf("\n[main raw] ESKF replay başlıyor...\n");

for k = 2:N

    dt_k = sim.t(k) - sim.t(k-1);

    if dt_k <= 0 || dt_k > 0.2
        warning("[main raw] Bad dt at k=%d: %.6f s, propagation skipped.", k, dt_k);
        continue;
    end

    %% ---------------- GPS measurement schedule ----------------
    use_gps_meas_now = is_gps_schedule_active(sim.t(k), GPS_MEAS_MODE, GPS_MEAS_WINDOWS);

    if prev_use_gps_meas && ~use_gps_meas_now
        baro_bias_hold_value = state.b_baro;
        fprintf('[main raw] t=%.2f s: GPS OFF, freezing b_baro = %.6f m\n', ...
            sim.t(k), baro_bias_hold_value);
    end

    if ~prev_use_gps_meas && use_gps_meas_now
        fprintf('[main raw] t=%.2f s: GPS ON, re-enabling b_baro estimation\n', ...
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

    %% ---------------- GPS / RTK position update ----------------
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

    %% ---------------- Pseudo GPS velocity update ----------------
    if USE_PSEUDO_GPS_VEL && use_gps_meas_now && sim.gps_vel_available(k)
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

fprintf("[main raw] ESKF replay bitti.\n");

%% ============================================================
% 10) Plot results
% ============================================================

plot_results_raw(sim, t, log_p, log_v, log_q, log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params);

figure('Name','Raw GPS Schedule and Baro Bias');
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
title('Raw/Quad GPS Schedule and Barometer Bias Observability');

figure('Name','Raw Innovation Gate Flags');
plot(t, double(log_gate_gps_pos), 'LineWidth', 1.0); hold on;
plot(t, double(log_gate_gps_vel), 'LineWidth', 1.0);
plot(t, double(log_gate_baro), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Gate rejected flag');
legend('GPS pos','GPS vel','Baro','Location','best');
title('Raw/Quad Innovation Gate Rejections');

%% ============================================================
% 11) Error metrics
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

fprintf('\n--- RAW RUN FINISHED ---\n');
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

fprintf('GPS pos skipped by schedule : %d\n', sum(idx_gps_pos_holdout));
fprintf('GPS vel skipped by schedule : %d\n', sum(idx_gps_vel_holdout));
fprintf('GPS measurement scheduled samples : %d\n', sum(log_use_gps_meas));
fprintf('Baro bias estimated samples       : %d\n', sum(log_estimate_baro_bias));

fprintf('\n--- RESIDUAL SUMMARY ---\n');
gps_pos_res_norm = vecnorm(log_res_gps_pos, 2, 1);
gps_vel_res_norm = vecnorm(log_res_gps_vel, 2, 1);

fprintf('Mean GPS pos residual norm [m]     : %.4f\n', mean(gps_pos_res_norm, 'omitnan'));
fprintf('Max GPS pos residual norm [m]      : %.4f\n', max(gps_pos_res_norm, [], 'omitnan'));
fprintf('Mean GPS vel residual norm [m/s]   : %.4f\n', mean(gps_vel_res_norm, 'omitnan'));
fprintf('Max GPS vel residual norm [m/s]    : %.4f\n', max(gps_vel_res_norm, [], 'omitnan'));
fprintf('Mean baro residual abs [m]         : %.4f\n', mean(abs(log_res_baro), 'omitnan'));
fprintf('Max baro residual abs [m]          : %.4f\n', max(abs(log_res_baro), [], 'omitnan'));

fprintf('\n--- FINAL STATES ---\n');
fprintf('Final estimated p [m]       : [%.4f %.4f %.4f]^T\n', log_p(:,end));
fprintf('Final truth p [m]           : [%.4f %.4f %.4f]^T\n', sim.p_true(:,end));
fprintf('Final estimated v [m/s]     : [%.4f %.4f %.4f]^T\n', log_v(:,end));
fprintf('Final truth v [m/s]         : [%.4f %.4f %.4f]^T\n', sim.v_true(:,end));
fprintf('Final gyro bias [rad/s]     : [%.6f %.6f %.6f]^T\n', log_bg(:,end));
fprintf('Final accel bias [m/s^2]    : [%.6f %.6f %.6f]^T\n', log_ba(:,end));
fprintf('Final baro offset [m]       : %.6f\n', log_bbaro(end));
fprintf('Final wind NE [m/s]         : [%.6f %.6f]^T\n', log_wind_ne(:,end));

fprintf('\n--- CONFIG USED ---\n');
fprintf('Sensor profile : %s\n', params.sensor_profile);
fprintf('USE_BARO       : %d\n', USE_BARO);
fprintf('USE_PSEUDO_GPS_VEL : %d\n', USE_PSEUDO_GPS_VEL);
fprintf('USE_TRUTH_ATT_INIT_DEBUG : %d\n', USE_TRUTH_ATT_INIT_FOR_DEBUG);
fprintf('GPS_MEAS_MODE  : %s\n', GPS_MEAS_MODE);
fprintf('GPS_HEIGHT_ANCHOR_TIMEOUT_S : %.3f\n', GPS_HEIGHT_ANCHOR_TIMEOUT_S);
fprintf('GPS_OUTAGE_SIGMA_A_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_A_SCALE);
fprintf('GPS_OUTAGE_SIGMA_G_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_G_SCALE);
fprintf('GPS_OUTAGE_SIGMA_BA_RW_SCALE : %.3f\n', GPS_OUTAGE_SIGMA_BA_RW_SCALE);
fprintf('sigma_gps_pos : [%.3f %.3f %.3f]\n', params.sigma_gps_pos);
fprintf('sigma_gps_vel : [%.3f %.3f %.3f]\n', params.sigma_gps_vel);
fprintf('sigma_baro    : %.3f\n', params.sigma_baro);

%% ============================================================
% 12) Save run output
% ============================================================

if SAVE_RUN_OUTPUT
    fprintf("\n[main raw] Saving run output to %s ...\n", save_file);

    save(save_file, ...
        "real_data_file", "sim", "t", ...
        "log_p", "log_v", "log_q", ...
        "log_bg", "log_ba", "log_bbaro", "log_wind_ne", "log_Pdiag", ...
        "log_res_gps_pos", "log_res_gps_vel", "log_res_baro", ...
        "log_use_gps_meas", "log_estimate_baro_bias", ...
        "log_gate_gps_pos", "log_gate_gps_vel", "log_gate_baro", ...
        "used_gps_pos_updates", "used_gps_vel_updates", "used_baro_updates", ...
        "accepted_gps_pos_updates", "accepted_gps_vel_updates", "accepted_baro_updates", ...
        "rejected_gps_pos_updates", "rejected_gps_vel_updates", "rejected_baro_updates", ...
        "USE_BARO", "USE_PSEUDO_GPS_VEL", "USE_TRUTH_ATT_INIT_FOR_DEBUG", ...
        "GPS_MEAS_MODE", "GPS_MEAS_WINDOWS", ...
        "GPS_HEIGHT_ANCHOR_TIMEOUT_S", ...
        "ESTIMATE_BARO_BIAS_WITH_GPS", ...
        "GPS_OUTAGE_SIGMA_A_SCALE", "GPS_OUTAGE_SIGMA_G_SCALE", ...
        "GPS_OUTAGE_SIGMA_BA_RW_SCALE", ...
        "params", "-v7.3");

    fprintf("[main raw] Saved: %s\n", save_file);
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