clear; clc; close all;

%% ============================================================
% FIXED-WING REAL DATA MODE - 18 STATE ESKF
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
%   - PX4 fixed-wing ULog-derived MAT
%   - IMU + GPS position + GPS velocity + baro + TAS
%   - GPS schedule / dropout support
%   - baro bias estimated only when recent GPS height exists
%   - GPS outage: nominal b_baro freeze, no manual P manipulation
%   - GPS outage: process-noise scaling, covariance grows through Qd
%   - innovation gates
%   - observability analysis on estimated nominal trajectory
% ============================================================

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

USE_OBSERVABILITY_ANALYSIS = true;

SAVE_RUN_OUTPUT = true;
save_file = "fixedwing_run_gps_dropout_observability.mat";

%% ---------------- GPS MEASUREMENT SCHEDULE ----------------
% "always"  : GPS position + GPS velocity tüm uçuş boyunca açık
% "never"   : GPS position + GPS velocity tüm uçuş boyunca kapalı
% "windows" : sadece verilen zaman pencerelerinde açık

GPS_MEAS_MODE = "always";

% Resetlenmiş sim zamanı üzerinden saniye cinsinden.
% Örnek: 120-150 s arası GPS kesintisi.
GPS_MEAS_WINDOWS = [
    0    120
    150  375
];

%% ---------------- BARO BIAS OBSERVABILITY ----------------
% Baro bias sadece GPS position height yakın zamanda geldiyse estimate edilir.
% GPS yokken b_baro nominal state olarak son değerde freeze kalır.

ESTIMATE_BARO_BIAS_WITH_GPS = true;

% GPS 5 Hz civarı. 0.5 s, yaklaşık 2-3 GPS position sample aralığı.
GPS_HEIGHT_ANCHOR_TIMEOUT_S = 0.50;

%% ---------------- GPS OUTAGE PROCESS-NOISE HANDLING ----------------
% P'ye elle müdahale yok.
% GPS yokken covariance Qd üzerinden kendi büyüsün diye process noise artırılır.

GPS_OUTAGE_SIGMA_A_SCALE     = 3.0;
GPS_OUTAGE_SIGMA_G_SCALE     = 1.5;
GPS_OUTAGE_SIGMA_BA_RW_SCALE = 4.0;

%% ============================================================
% 1) Init filter
% ============================================================

[state, ~, params] = init_filter();

%% ============================================================
% 1.5) Fixed-wing 18-state tuning, dataset-informed
% ============================================================

params.sensor_profile = 'PX4 Fixed-Wing ULog: IMU + GPS + baro + TAS';

%% Frame
params.g_n = [0; 0; 9.81];

%% IMU process noise
% Dataset: fixed-wing ULog, approx 250 Hz IMU.
% Bunlar filtre propagation tarafını belirler, ölçüm residuallarını birebir
% ezmek için kullanılmaz.
params.sigma_g = deg2rad(0.20);
params.sigma_a = 0.080;

params.sigma_bg_rw = deg2rad(0.010);
params.sigma_ba_rw = 0.010;

%% GPS measurement noise, dataset-informed
% Önceki residuallara göre generic [3 3 5] biraz gevşek kalıyordu.
% Bu değerler GPS'i hâlâ "mutlak referans" yapmıyor ama fazla da salmıyor.
params.sigma_gps_pos = [2.0; 2.0; 3.0];      % [N; E; D] [m]
params.R_gps_pos = diag(params.sigma_gps_pos.^2);

% GPS velocity residualları ~0.2-0.3 m/s bandında.
% Biraz pay bırakıyoruz.
params.sigma_gps_vel = [0.35; 0.35; 0.50];   % [N; E; D] [m/s]
params.R_gps_vel = diag(params.sigma_gps_vel.^2);

%% Barometer measurement noise / bias
% Raw baro offset büyük, ama b_baro state'i ile modeled residual düşük.
% sigma_baro çok küçük olursa D eksenini baroya yapıştırır.
params.sigma_baro = 1.5;                     % [m]
params.R_baro = params.sigma_baro^2;

params.sigma_baro_bias0 = 10.0;              % [m]
params.sigma_baro_bias_rw = 0.005;           % [m/sqrt(s)]

params.max_baro_update_rate_hz = 20;
params.estimate_baro_bias = false;

%% TAS / wind measurement tuning
% TAS residualında sistematik mean olabilir, o yüzden residual RMSE kadar
% düşük sigma vermiyoruz. Fazla güvenmek wind'i TAS modeline zorlar.
params.sigma_tas = 1.8;                      % [m/s]
params.R_tas = params.sigma_tas^2;

params.sigma_wind0 = [5.0; 5.0];             % [m/s]
params.sigma_wind_rw = 0.015;                % [m/sqrt(s)]

params.max_airspeed_update_rate_hz = 30;

params.use_joseph_form = true;

%% Innovation gates
params.gps_pos_gate_chi2 = 25.0;
params.gps_vel_gate_chi2 = 25.0;
params.baro_gate_chi2    = 9.0;
params.tas_gate_chi2     = 9.0;

%% Observability options
params.obs.enable = USE_OBSERVABILITY_ANALYSIS;

params.obs.sample_step = 1000;
params.obs.window_sec  = 3.0;
params.obs.rank_tol    = 1e-8;

params.obs.show_local = false;
params.obs.show_gps_schedule = true;
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

% Observability analysis de aynı GPS schedule'ı görsün.
sim.gps_meas_scheduled = tmp_use_gps;

fprintf('Expected GPS active time [s]   : %.3f\n', ...
    sum(tmp_use_gps) * median(diff(sim.t), "omitnan"));
fprintf('Expected GPS active samples    : %d / %d\n', sum(tmp_use_gps), N);

%% ============================================================
% 2.5) Barometer bias initial value
% ============================================================

% Başta kalibrasyon yok.
% b_baro = 0 başlar, GPS measurement pencerelerinde öğrenilir,
% GPS gidince nominal state son değerde freeze edilir.

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

if USE_ATT_INIT_FOR_DEBUG && isfield(sim, "q_ref") && ~isempty(sim.q_ref)
    state.q_nb = sim.q_ref(:,1);
    state.q_nb = state.q_nb / norm(state.q_nb);
    fprintf('[fixedwing main] DEBUG: attitude q0 log attitude üzerinden başlatıldı.\n');
else
    state.q_nb = [1;0;0;0];
    fprintf('[fixedwing main] REAL-LIKE: attitude q0 identity seçildi.\n');
end

state.b_g = [0;0;0];
state.b_a = [0;0;0];

state.b_baro = b_baro0;
state.b_baro = b_baro0;

%% ---------------- Initial wind estimate from GPS velocity + TAS ----------------
% Pitot/TAS tek başına wind vector vermez.
% Ama GPS ground velocity + TAS + attitude ile ilk wind_NE tahmini yapılabilir.
%
% Model:
%   v_ground^n = v_air^n + wind^n
%   v_air^n ≈ TAS * body_x^n
%   wind_NE ≈ v_ground_NE - TAS * body_x_NE
%
% Not:
%   PX4 wind_ref burada KULLANILMIYOR.
%   Çünkü PX4 wind_ref daha sonra karşılaştırma referansı olarak kullanılacak.

WIND_INIT_FROM_TAS_GPS = true;
WIND_INIT_WINDOW_S     = 10.0;   % ilk 10 saniye
WIND_INIT_TAS_MAX_DT_S = 0.10;   % TAS nearest sample toleransı
WIND_INIT_MAX_NORM     = 20.0;   % saçma adayları ele

state.wind_ne = [0;0];
wind_init_success = false;
wind_candidates = [];

if WIND_INIT_FROM_TAS_GPS && USE_AIRSPEED && ...
        isfield(sim, "airspeed_available") && ...
        isfield(sim, "airspeed") && ...
        isfield(sim, "gps_vel_available") && ...
        isfield(sim, "gps_vel")

    idx_gps_init = find( ...
        sim.t <= WIND_INIT_WINDOW_S & ...
        sim.gps_vel_available & ...
        all(isfinite(sim.gps_vel), 1));

    for ii = 1:numel(idx_gps_init)
        kk = idx_gps_init(ii);

        % GPS schedule kapalıysa initialization için de kullanma.
        % Bizim dropout testinde başlangıçta GPS açık olduğu için sorun yok.
        if ~is_gps_schedule_active(sim.t(kk), GPS_MEAS_MODE, GPS_MEAS_WINDOWS)
            continue;
        end

        z_vg = sim.gps_vel(:,kk);

        [z_tas, tas_ok] = nearest_airspeed_for_wind_init_local( ...
            sim, kk, WIND_INIT_TAS_MAX_DT_S);

        [q_nb_k, q_ok] = attitude_for_wind_init_local(sim, kk, state);

        if ~tas_ok || ~q_ok
            continue;
        end

        if ~isfinite(z_tas) || z_tas <= 3
            continue;
        end

        q_nb_k = q_nb_k(:) / norm(q_nb_k);
        R_nb_k = quat_to_rotmat(q_nb_k);

        % body x-axis in navigation frame
        body_x_n = R_nb_k(:,1);

        if ~all(isfinite(body_x_n)) || norm(body_x_n(1:2)) < 0.2
            continue;
        end

        w_ne_k = z_vg(1:2) - z_tas * body_x_n(1:2);

        if all(isfinite(w_ne_k)) && norm(w_ne_k) < WIND_INIT_MAX_NORM
            wind_candidates = [wind_candidates, w_ne_k]; %#ok<AGROW>
        end
    end

    if ~isempty(wind_candidates)
        state.wind_ne = median(wind_candidates, 2, "omitnan");
        wind_init_success = true;

        fprintf('[fixedwing main] Initial wind NE from GPS vel + TAS = [%.3f %.3f]^T m/s\n', ...
            state.wind_ne(1), state.wind_ne(2));

        fprintf('[fixedwing main] Wind init candidate count = %d\n', ...
            size(wind_candidates,2));

        % Initial wind artık sıfırdan sallanmıyor, o yüzden covariance'ı biraz daraltabiliriz.
        % Çok küçültmüyoruz; çünkü TAS yön modeli sideslip/AoA ihmal ediyor.
        params.sigma_wind0 = [3.0; 3.0];

    else
        fprintf('[fixedwing main] Initial wind NE = [0 0]^T, no valid GPS+TAS init samples.\n');
    end
else
    fprintf('[fixedwing main] Initial wind NE = [0 0]^T, TAS/GPS wind init disabled or unavailable.\n');
end

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

    if prev_use_gps_meas && ~use_gps_meas_now
        baro_bias_hold_value = state.b_baro;

        fprintf('[fixedwing main] t=%.2f s: GPS OFF, freezing b_baro = %.6f m\n', ...
            sim.t(k), baro_bias_hold_value);
    end

    if ~prev_use_gps_meas && use_gps_meas_now
        fprintf('[fixedwing main] t=%.2f s: GPS ON, re-enabling b_baro estimation\n', ...
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
    imu.gyro_m  = sim.imu_gyro(:,k);
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
    if use_gps_meas_now && sim.gps_vel_available(k)
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

    %% ---------------- Baro update ----------------
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
% 7.5) Observability analysis on real fixed-wing estimated trajectory
% ============================================================

obs = struct();
obs.available = false;

if isfield(params, 'obs') && isfield(params.obs, 'enable') && params.obs.enable

    % Real data'da truth yok.
    % Bu yüzden observability analizi estimated nominal trajectory etrafında yapılır.
    obs_sim = sim;

    obs_sim.p_true = log_p;
    obs_sim.v_true = log_v;
    obs_sim.q_true = log_q;

    obs_sim.bg_true = log_bg;
    obs_sim.ba_true = log_ba;
    obs_sim.b_baro_true = log_bbaro;

    % TAS ölçümünün wind Jacobian'ı için estimated wind kullanılır.
    obs_sim.wind_ref = log_wind_ne;
    obs_sim.wind_ref_available = all(isfinite(log_wind_ne), 1);

    % GPS schedule ve baro bias observability flagleri.
    obs_sim.gps_meas_scheduled = log_use_gps_meas;
    obs_sim.estimate_baro_bias_scheduled = log_estimate_baro_bias;

    obs = analyze_observability(obs_sim, params);
    obs.available = true;

    valid_obs = obs.t_eval <= (sim.t(end) - params.obs.window_sec);

    gps_on_obs = interp1(sim.t, double(log_use_gps_meas), ...
        obs.t_eval, 'nearest', 'extrap') > 0.5;

    gps_off_obs = ~gps_on_obs;

    valid_on  = valid_obs & gps_on_obs;
    valid_off = valid_obs & gps_off_obs;

    fprintf('\n--- FIXED-WING OBSERVABILITY SUMMARY ---\n');
    fprintf('State dimension              : %d\n', obs.n_err);
    fprintf('Mean local rank              : %.2f\n', mean(obs.local_rank, 'omitnan'));
    fprintf('Mean windowed rank           : %.2f\n', mean(obs.window_rank, 'omitnan'));
    fprintf('Min local sigma_min          : %.3e\n', min(obs.local_sigma_min));
    fprintf('Min windowed sigma_min       : %.3e\n', min(obs.window_sigma_min));

    if any(valid_obs)
        fprintf('Mean windowed rank, valid    : %.2f\n', ...
            mean(obs.window_rank(valid_obs), 'omitnan'));
        fprintf('Min windowed sigma min, valid: %.3e\n', ...
            min(obs.window_sigma_min(valid_obs)));
    end

    fprintf('\n--- OBSERVABILITY BY GPS AVAILABILITY ---\n');

    if any(valid_on)
        fprintf('Windowed rank mean, GPS ON        : %.2f\n', ...
            mean(obs.window_rank(valid_on), 'omitnan'));
        fprintf('Windowed rank min,  GPS ON        : %.2f\n', ...
            min(obs.window_rank(valid_on)));
        fprintf('Windowed sigma min, GPS ON        : %.3e\n', ...
            min(obs.window_sigma_min(valid_on)));
    else
        fprintf('No valid GPS ON observability samples.\n');
    end

    if any(valid_off)
        fprintf('Windowed rank mean, GPS OFF       : %.2f\n', ...
            mean(obs.window_rank(valid_off), 'omitnan'));
        fprintf('Windowed rank min,  GPS OFF       : %.2f\n', ...
            min(obs.window_rank(valid_off)));
        fprintf('Windowed sigma min, GPS OFF       : %.3e\n', ...
            min(obs.window_sigma_min(valid_off)));
    else
        fprintf('No valid GPS OFF observability samples.\n');
    end

else
    fprintf('\n--- Observability analysis skipped ---\n');
    fprintf('Set USE_OBSERVABILITY_ANALYSIS = true to enable it.\n');
end

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

if isfield(sim, "airspeed_available")
    fprintf('Airspeed updates available : %d\n', sum(sim.airspeed_available));
end

fprintf('Airspeed updates used      : %d\n', used_airspeed_updates);
fprintf('Airspeed updates accepted  : %d\n', accepted_airspeed_updates);
fprintf('Airspeed updates rejected  : %d\n', rejected_airspeed_updates);

fprintf('GPS pos skipped by schedule : %d\n', sum(idx_gps_pos_holdout));
fprintf('GPS vel skipped by schedule : %d\n', sum(idx_gps_vel_holdout));
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
fprintf('GPS_OUTAGE_SIGMA_A_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_A_SCALE);
fprintf('GPS_OUTAGE_SIGMA_G_SCALE     : %.3f\n', GPS_OUTAGE_SIGMA_G_SCALE);
fprintf('GPS_OUTAGE_SIGMA_BA_RW_SCALE : %.3f\n', GPS_OUTAGE_SIGMA_BA_RW_SCALE);
fprintf('gps_pos_gate_chi2 : %.3f\n', params.gps_pos_gate_chi2);
fprintf('gps_vel_gate_chi2 : %.3f\n', params.gps_vel_gate_chi2);
fprintf('baro_gate_chi2    : %.3f\n', params.baro_gate_chi2);
fprintf('tas_gate_chi2     : %.3f\n', params.tas_gate_chi2);
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

    plot(th, err_hold(1,:), '.');
    hold on;
    plot(th, err_hold(2,:), '.');
    plot(th, err_hold(3,:), '.');

    grid on;
    xlabel('Time [s]');
    ylabel('Position error wrt held-out GPS [m]');
    legend('N','E','D');
    title('GPS Position Holdout Error During GPS Dropout');
end

figure('Name','Innovation Gate Flags');

plot(t, double(log_gate_gps_pos), 'LineWidth', 1.0);
hold on;
plot(t, double(log_gate_gps_vel), 'LineWidth', 1.0);
plot(t, double(log_gate_baro), 'LineWidth', 1.0);
plot(t, double(log_gate_airspeed), 'LineWidth', 1.0);

grid on;
xlabel('Time [s]');
ylabel('Gate rejected flag');
legend('GPS pos','GPS vel','Baro','TAS','Location','best');
title('Innovation Gate Rejections');

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
        "log_gate_gps_pos", "log_gate_gps_vel", "log_gate_baro", "log_gate_airspeed", ...
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
        "accepted_gps_pos_updates", "accepted_gps_vel_updates", ...
        "accepted_baro_updates", "accepted_airspeed_updates", ...
        "rejected_gps_pos_updates", "rejected_gps_vel_updates", ...
        "rejected_baro_updates", "rejected_airspeed_updates", ...
        "USE_BARO", "USE_AIRSPEED", "USE_ATT_INIT_FOR_DEBUG", ...
        "GPS_MEAS_MODE", "GPS_MEAS_WINDOWS", ...
        "GPS_HEIGHT_ANCHOR_TIMEOUT_S", ...
        "ESTIMATE_BARO_BIAS_WITH_GPS", ...
        "GPS_OUTAGE_SIGMA_A_SCALE", "GPS_OUTAGE_SIGMA_G_SCALE", ...
        "GPS_OUTAGE_SIGMA_BA_RW_SCALE", ...
        "USE_COMBINED_PX4_REFERENCE", ...
        "USE_OBSERVABILITY_ANALYSIS", ...
        "params", "px4_cmp", "obs", "-v7.3");

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
function [z_tas, ok] = nearest_airspeed_for_wind_init_local(sim, k, max_dt_s)
%NEAREST_AIRSPEED_FOR_WIND_INIT_LOCAL
% Wind initialization için GPS velocity sample'ına en yakın TAS sample'ını bulur.

z_tas = nan;
ok = false;

if ~isfield(sim, "airspeed_available") || ~isfield(sim, "airspeed")
    return;
end

idx_tas = find(sim.airspeed_available & isfinite(sim.airspeed));

if isempty(idx_tas)
    return;
end

[dt_min, jj] = min(abs(sim.t(idx_tas) - sim.t(k)));

if isempty(jj) || ~isfinite(dt_min) || dt_min > max_dt_s
    return;
end

kk_tas = idx_tas(jj);
z_tas = sim.airspeed(kk_tas);

ok = isfinite(z_tas);
end

function [q_nb, ok] = attitude_for_wind_init_local(sim, k, state)
%ATTITUDE_FOR_WIND_INIT_LOCAL
% Wind initialization için attitude alır.
% Öncelik:
%   1) sim.q_ref(:,k)
%   2) state.q_nb

q_nb = nan(4,1);
ok = false;

if isfield(sim, "q_ref") && size(sim.q_ref,2) >= k
    q_try = sim.q_ref(:,k);

    if numel(q_try) == 4 && all(isfinite(q_try)) && norm(q_try) > 1e-9
        q_nb = q_try(:) / norm(q_try);
        ok = true;
        return;
    end
end

if isfield(state, "q_nb")
    q_try = state.q_nb(:);

    if numel(q_try) == 4 && all(isfinite(q_try)) && norm(q_try) > 1e-9
        q_nb = q_try / norm(q_try);
        ok = true;
        return;
    end
end
end