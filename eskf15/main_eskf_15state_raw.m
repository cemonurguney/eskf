clear; clc; close all;

%% ============================================================
%  RAW MRS RECTANGLE MODE - 15 STATE ESKF
% ============================================================
% Inputs:
%   - raw IMU        : imu_gyro, imu_accel
%   - raw RTK/GPS    : gps_pos_enu
%   - pseudo GPS vel : smoothed gps_pos_enu finite difference
%   - raw baro       : baro_alt
%
% Reference / truth:
%   - /rtk_fused_odom -> sim.p_true, sim.v_true, sim.q_true
%
% Important:
%   Truth is used only for plotting/error metrics.
%   Truth is NOT used as measurement update.
% ============================================================

real_data_file = "rectangle_raw_clean.mat";

%% ---------------- USER FLAGS ----------------
USE_BARO = true;
USE_PSEUDO_GPS_VEL = true;

% Debug attitude initialization:
% true:
%   q0 truth/reference üzerinden başlatılır.
%   Bu measurement update değildir, sadece frame/veri hattını debug etmek için.
%
% false:
%   q0 identity alınır.
%   Gerçek sisteme daha yakın ama yaw başlangıcı bilinmezse sonuç bozulabilir.
USE_TRUTH_ATT_INIT_FOR_DEBUG = true;

% Initial velocity:
% true:
%   ilk pseudo GPS velocity varsa v0 olarak kullanır.
%
% false:
%   v0 = [0;0;0]
USE_PSEUDO_GPS_VEL_INIT = true;

fprintf('[main raw] Raw gerçek veri modu: %s yükleniyor...\n', real_data_file);

%% ============================================================
% 1) Init filter
% ============================================================
[state, P, params] = init_filter();

% ENU convention:
% p = [E; N; U]
% gravity is down in ENU
params.g_n = [0; 0; -9.81];

%% ============================================================
% 2) Measurement covariance tuning
% ============================================================

% Raw RTK/GPS position'a daha fazla güven.
% Önceki testte GPS güveni artırınca position error düştü.
params.sigma_gps_pos = [0.4; 0.4; 0.6];
params.R_gps_pos = diag(params.sigma_gps_pos.^2);

% Pseudo GPS velocity, gerçek Doppler velocity değil.
% Raw GPS position türevinden geldiği için fazla güvenme.
% Eğer velocity update sistemi titretirse bunları büyüt:
%   [0.7; 0.7; 1.2]
params.sigma_gps_vel = [0.45; 0.45; 0.90];
params.R_gps_vel = diag(params.sigma_gps_vel.^2);

% Baro güveni azaltıldı.
% RTK altitude daha iyi göründüğü için baro sadece zayıf yardımcı olacak.
params.sigma_baro = 3.0;
params.R_baro = params.sigma_baro^2;

%% ============================================================
% 3) Build sim from raw MAT
% ============================================================
sim = build_sim_from_raw_mat(real_data_file);

t = sim.t;
dt = sim.dt;
N = numel(t);

%% ============================================================
% 4) Start index
% ============================================================
% Şimdilik baştan başlıyoruz.
% İleride gerçek sistem için pre-launch/static detection buraya eklenebilir.
k_start = 1;

fprintf('[main raw] Başlangıç indeksi: k_start = %d\n', k_start);

%% ============================================================
% 5) Trim sim if needed
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

    sim.range = sim.range(k_start:end);
    sim.range_available = sim.range_available(k_start:end);
end

t = sim.t;
dt = sim.dt;
N = numel(t);

%% ============================================================
% 6) Initial state
% ============================================================

% ---------------- Position init ----------------
% İlk raw GPS/RTK ölçümünden başlat.
idx_gps0 = find(sim.gps_pos_available & all(isfinite(sim.gps_pos),1), 1, 'first');

if isempty(idx_gps0)
    warning('[main raw] GPS ölçümü yok, p0 = [0;0;0] kullanılıyor.');
    state.p_n = [0;0;0];
else
    state.p_n = sim.gps_pos(:,idx_gps0);
end

% ---------------- Velocity init ----------------
idx_gps_vel0 = find(sim.gps_vel_available & all(isfinite(sim.gps_vel),1), 1, 'first');

if USE_PSEUDO_GPS_VEL_INIT && ~isempty(idx_gps_vel0)
    state.v_n = sim.gps_vel(:,idx_gps_vel0);
    fprintf('[main raw] Initial velocity pseudo GPS velocity üzerinden başlatıldı.\n');
else
    state.v_n = [0;0;0];
    fprintf('[main raw] Initial velocity v0 = [0;0;0] seçildi.\n');
end

% ---------------- Attitude init ----------------
if USE_TRUTH_ATT_INIT_FOR_DEBUG
    state.q_nb = sim.q_true(:,1);
    state.q_nb = state.q_nb / norm(state.q_nb);
    fprintf('[main raw] DEBUG: attitude q0 truth/reference üzerinden başlatıldı.\n');
else
    state.q_nb = [1;0;0;0];
    fprintf('[main raw] REAL-LIKE: attitude q0 identity seçildi. Yaw sapması beklenebilir.\n');
end

% ---------------- Bias init ----------------
state.b_g = [0;0;0];
state.b_a = [0;0;0];

fprintf('[main raw] Initial p = [%.3f %.3f %.3f]^T m\n', state.p_n);
fprintf('[main raw] Initial v = [%.3f %.3f %.3f]^T m/s\n', state.v_n);
fprintf('[main raw] USE_BARO = %d\n', USE_BARO);
fprintf('[main raw] USE_PSEUDO_GPS_VEL = %d\n', USE_PSEUDO_GPS_VEL);

%% ============================================================
% 7) Logs
% ============================================================
log_p = zeros(3, N);
log_v = zeros(3, N);
log_q = zeros(4, N);
log_bg = zeros(3, N);
log_ba = zeros(3, N);
log_Pdiag = zeros(15, N);

log_res_gps_pos = nan(3, N);
log_res_gps_vel = nan(3, N);
log_res_baro = nan(1, N);

used_gps_pos_updates = 0;
used_gps_vel_updates = 0;
used_baro_updates = 0;

%% ============================================================
% 8) Main ESKF loop
% ============================================================
for k = 1:N

    %% ---------------- IMU measurement ----------------
    imu.gyro_m = sim.imu_gyro(:,k);
    imu.accel_m = sim.imu_accel(:,k);

    %% ---------------- Nominal propagation ----------------
    state = propagate_nominal(state, imu, params, dt);

    %% ---------------- Covariance propagation ----------------
    [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt);
    P = propagate_covariance(P, F, G, Qd, dt);

    %% ---------------- Raw GPS/RTK position update ----------------
    if sim.gps_pos_available(k)
        z_gps_pos = sim.gps_pos(:,k);

        if all(isfinite(z_gps_pos))
            [state, P, residual, ~, ~] = update_gnss_pos(state, P, z_gps_pos, params);
            log_res_gps_pos(:,k) = residual;

            used_gps_pos_updates = used_gps_pos_updates + 1;
        end
    end

    %% ---------------- Pseudo GPS velocity update ----------------
    if USE_PSEUDO_GPS_VEL && sim.gps_vel_available(k)
        z_gps_vel = sim.gps_vel(:,k);

        if all(isfinite(z_gps_vel))
            [state, P, residual, ~, ~] = update_gnss_vel(state, P, z_gps_vel, params);
            log_res_gps_vel(:,k) = residual;

            used_gps_vel_updates = used_gps_vel_updates + 1;
        end
    end

    %% ---------------- Raw barometer update ----------------
    if USE_BARO && sim.baro_available(k)
        z_baro = sim.baro(k);

        if isfinite(z_baro)
            [state, P, residual, ~, ~] = update_baro(state, P, z_baro, params);
            log_res_baro(k) = residual;

            used_baro_updates = used_baro_updates + 1;
        end
    end

    %% ---------------- Log ----------------
    log_p(:,k) = state.p_n;
    log_v(:,k) = state.v_n;
    log_q(:,k) = state.q_nb;
    log_bg(:,k) = state.b_g;
    log_ba(:,k) = state.b_a;
    log_Pdiag(:,k) = diag(P);
end

%% ============================================================
% 9) Plot results
% ============================================================
plot_results_raw(sim, t, log_p, log_v, log_q, log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params);

%% ============================================================
% 10) Error metrics
% ============================================================
err_p = log_p - sim.p_true;
err_v = log_v - sim.v_true;

pos_err_norm = vecnorm(err_p, 2, 1);
vel_err_norm = vecnorm(err_v, 2, 1);

fprintf('\n--- RAW RUN FINISHED ---\n');

fprintf('Mean pos error norm [m]    : %.4f\n', mean(pos_err_norm, 'omitnan'));
fprintf('Final pos error norm [m]   : %.4f\n', norm(err_p(:,end)));
fprintf('Max pos error norm [m]     : %.4f\n', max(pos_err_norm, [], 'omitnan'));

fprintf('Mean vel error norm [m/s]  : %.4f\n', mean(vel_err_norm, 'omitnan'));
fprintf('Final vel error norm [m/s] : %.4f\n', norm(err_v(:,end)));
fprintf('Max vel error norm [m/s]   : %.4f\n', max(vel_err_norm, [], 'omitnan'));

fprintf('\n--- UPDATE COUNTS ---\n');
fprintf('GPS pos updates available   : %d\n', sum(sim.gps_pos_available));
fprintf('GPS pos updates used        : %d\n', used_gps_pos_updates);

fprintf('GPS vel updates available   : %d\n', sum(sim.gps_vel_available));
fprintf('GPS vel updates used        : %d\n', used_gps_vel_updates);

fprintf('BARO updates available      : %d\n', sum(sim.baro_available & isfinite(sim.baro)));
fprintf('BARO updates used           : %d\n', used_baro_updates);

fprintf('\n--- RESIDUAL SUMMARY ---\n');

gps_pos_res_norm = vecnorm(log_res_gps_pos, 2, 1);
gps_vel_res_norm = vecnorm(log_res_gps_vel, 2, 1);

fprintf('Mean GPS pos residual norm [m]    : %.4f\n', mean(gps_pos_res_norm, 'omitnan'));
fprintf('Max  GPS pos residual norm [m]    : %.4f\n', max(gps_pos_res_norm, [], 'omitnan'));

fprintf('Mean GPS vel residual norm [m/s]  : %.4f\n', mean(gps_vel_res_norm, 'omitnan'));
fprintf('Max  GPS vel residual norm [m/s]  : %.4f\n', max(gps_vel_res_norm, [], 'omitnan'));

fprintf('Mean baro residual abs [m]        : %.4f\n', mean(abs(log_res_baro), 'omitnan'));
fprintf('Max  baro residual abs [m]        : %.4f\n', max(abs(log_res_baro), [], 'omitnan'));

fprintf('\n--- FINAL STATES ---\n');
fprintf('Final estimated p [m]     : [%.4f %.4f %.4f]^T\n', log_p(:,end));
fprintf('Final truth p [m]         : [%.4f %.4f %.4f]^T\n', sim.p_true(:,end));
fprintf('Final estimated v [m/s]   : [%.4f %.4f %.4f]^T\n', log_v(:,end));
fprintf('Final truth v [m/s]       : [%.4f %.4f %.4f]^T\n', sim.v_true(:,end));

fprintf('Final gyro bias [rad/s]   : [%.6f %.6f %.6f]^T\n', log_bg(:,end));
fprintf('Final accel bias [m/s^2]  : [%.6f %.6f %.6f]^T\n', log_ba(:,end));

fprintf('\n--- CONFIG USED ---\n');
fprintf('USE_BARO                  : %d\n', USE_BARO);
fprintf('USE_PSEUDO_GPS_VEL        : %d\n', USE_PSEUDO_GPS_VEL);
fprintf('USE_TRUTH_ATT_INIT_DEBUG  : %d\n', USE_TRUTH_ATT_INIT_FOR_DEBUG);

fprintf('sigma_gps_pos             : [%.3f %.3f %.3f]\n', params.sigma_gps_pos);
fprintf('sigma_gps_vel             : [%.3f %.3f %.3f]\n', params.sigma_gps_vel);
fprintf('sigma_baro                : %.3f\n', params.sigma_baro);