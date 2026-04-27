clear; clc; close all;

%% ============================================================
%  RAW MRS RECTANGLE MODE
% ============================================================
real_data_file = "rectangle_raw_clean.mat";

USE_BARO = false;

% Kritik:
% true  -> ilk debug koşusunda attitude başlangıcı truth'tan alınır.
%          Bu, filtreye truth measurement vermek değildir ama gerçek sistem testi de değildir.
%
% false -> gerçek sisteme daha yakın: p0 GPS'ten, v0=0, q0 identity.
%          Yaw başlangıcı bilinmezse yaw/konum sapması beklenir.
%
% Önce true ile frame/veri hattını doğrula, sonra false'a geç.
USE_TRUTH_ATT_INIT_FOR_DEBUG = true;

fprintf('[main raw] Raw gercek veri modu: %s yukleniyor...\n', real_data_file);

%% 1) Init filter
[state, P, params] = init_filter();

% ENU convention:
% p = [E; N; U], gravity is down.
params.g_n = [0; 0; -9.81];

% Raw RTK genelde iyi ama yine de ilk koşuda gevşek bırak.
params.sigma_gps_pos = [0.4; 0.4; 0.6];
params.R_gps_pos = diag(params.sigma_gps_pos.^2);

% Baro yükseklik gürültüsü.
params.sigma_baro = 0.8;
params.R_baro = params.sigma_baro^2;

%% 2) Build sim from raw MAT
sim = build_sim_from_raw_mat(real_data_file);

t = sim.t;
dt = sim.dt;
N = numel(t);

%% 3) Start index
% Şimdilik baştan başlıyoruz.
% İleride burada IMU+GPS ile gerçek pre-launch init detection ekleriz.
k_start = 1;

fprintf('[main raw] Baslangic indeksi: k_start = %d\n', k_start);

%% 4) Trim sim
if k_start > 1
    t_offset = sim.t(k_start);

    sim.t = sim.t(k_start:end) - t_offset;
    sim.imu_gyro = sim.imu_gyro(:,k_start:end);
    sim.imu_accel = sim.imu_accel(:,k_start:end);

    sim.p_true = sim.p_true(:,k_start:end);
    sim.v_true = sim.v_true(:,k_start:end);
    sim.q_true = sim.q_true(:,k_start:end);

    sim.gps_pos = sim.gps_pos(:,k_start:end);
    sim.gps_vel = sim.gps_vel(:,k_start:end);
    sim.gps_pos_available = sim.gps_pos_available(k_start:end);
    sim.gps_vel_available = sim.gps_vel_available(k_start:end);

    sim.baro = sim.baro(k_start:end);
    sim.baro_available = sim.baro_available(k_start:end);

    sim.range = sim.range(k_start:end);
    sim.range_available = sim.range_available(k_start:end);
end

t = sim.t;
dt = sim.dt;
N = numel(t);

%% 5) Initial state
% Position: ilk raw GPS/RTK ölçümünden başlat.
idx_gps0 = find(sim.gps_pos_available, 1, 'first');

if isempty(idx_gps0)
    warning('[main raw] GPS ölçümü yok, p0=[0;0;0] kullaniliyor.');
    state.p_n = [0;0;0];
else
    state.p_n = sim.gps_pos(:,idx_gps0);
end

% Velocity: raw RTK velocity vermediği için gerçek sisteme yakın şekilde sıfır.
state.v_n = [0;0;0];

% Attitude initialization:
if USE_TRUTH_ATT_INIT_FOR_DEBUG
    state.q_nb = sim.q_true(:,1) / norm(sim.q_true(:,1));
    fprintf('[main raw] DEBUG: attitude q0 truth/reference üzerinden başlatıldı.\n');
else
    state.q_nb = [1;0;0;0];
    fprintf('[main raw] REAL-LIKE: attitude q0 identity seçildi. Yaw sapması beklenebilir.\n');
end

% Bias başlangıcı:
state.b_g = [0;0;0];
state.b_a = [0;0;0];

fprintf('[main raw] Initial p = [%.3f %.3f %.3f]^T m\n', state.p_n);
fprintf('[main raw] Initial v = [%.3f %.3f %.3f]^T m/s\n', state.v_n);

%% 6) Logs
log_p = zeros(3, N);
log_v = zeros(3, N);
log_q = zeros(4, N);
log_bg = zeros(3, N);
log_ba = zeros(3, N);
log_Pdiag = zeros(15, N);

log_res_gps_pos = nan(3, N);
log_res_gps_vel = nan(3, N);
log_res_baro = nan(1, N);

%% 7) Main EKF loop
for k = 1:N

    imu.gyro_m = sim.imu_gyro(:,k);
    imu.accel_m = sim.imu_accel(:,k);

    % Nominal propagation
    state = propagate_nominal(state, imu, params, dt);

    % Covariance propagation
    [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt);
    P = propagate_covariance(P, F, G, Qd, dt);

    % Raw GPS/RTK position update
    if sim.gps_pos_available(k)
        z_gps_pos = sim.gps_pos(:,k);

        if all(isfinite(z_gps_pos))
            [state, P, residual, ~, ~] = update_gnss_pos(state, P, z_gps_pos, params);
            log_res_gps_pos(:,k) = residual;
        end
    end

    % GPS velocity update disabled, because /rtk_raw does not directly provide velocity.
    % sim.gps_vel_available is intentionally false.

    % Raw barometer update
    if USE_BARO && sim.baro_available(k)
        z_baro = sim.baro(k);

        if isfinite(z_baro)
            [state, P, residual, ~, ~] = update_baro(state, P, z_baro, params);
            log_res_baro(k) = residual;
        end
    end

    % Log
    log_p(:,k) = state.p_n;
    log_v(:,k) = state.v_n;
    log_q(:,k) = state.q_nb;
    log_bg(:,k) = state.b_g;
    log_ba(:,k) = state.b_a;
    log_Pdiag(:,k) = diag(P);
end

%% 8) Results
plot_results_raw(sim, t, log_p, log_v, log_q, log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params);

%% 9) Error metrics
err_p = log_p - sim.p_true;
err_v = log_v - sim.v_true;

fprintf('\n--- RAW RUN FINISHED ---\n');
fprintf('Mean pos error norm [m]    : %.4f\n', mean(vecnorm(err_p,2,1), 'omitnan'));
fprintf('Final pos error norm [m]   : %.4f\n', norm(err_p(:,end)));
fprintf('Mean vel error norm [m/s]  : %.4f\n', mean(vecnorm(err_v,2,1), 'omitnan'));
fprintf('Final vel error norm [m/s] : %.4f\n', norm(err_v(:,end)));
fprintf('GPS pos updates            : %d\n', sum(sim.gps_pos_available));
fprintf('BARO updates               : %d\n', sum(sim.baro_available & isfinite(sim.baro)));