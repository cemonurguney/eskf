clear; clc; close all;

%% ============================================================
% Fixed-wing files
% MATLAB Online için tüm dosyaları aynı klasöre koy.
% Dosya adın farklıysa sadece burayı değiştir.
% ============================================================

fixedwing_file = "0002_20.48_fixedwing_eskf.mat";
combined_csv_file = "0002_20.48.csv";

%% ---------------- USER FLAGS ----------------
USE_BARO = true;
USE_GPS_VEL = true;

% İlk testte attitude debug init açık olsun.
% Sistem çalışınca false yapıp identity başlatmayı deneriz.
USE_ATT_INIT_FOR_DEBUG = true;

% PX4 reference kıyası için combined CSV kullan.
USE_COMBINED_PX4_REFERENCE = true;

% Run sonunda sonuçları .mat olarak kaydet.
SAVE_RUN_OUTPUT = true;
save_file = "fixedwing_run_with_px4_compare.mat";

%% ============================================================
% 1) Init filter
% ============================================================

[state, P, params] = init_filter();

%% ============================================================
% 1.5) Fixed-wing tuning
% ============================================================

params.sensor_profile = 'PX4 Fixed-Wing ULog: IMU + GPS + baro';

% PX4 fixed-wing local frame: NED
% p = [N; E; D]
% gravity positive Down
params.g_n = [0; 0; 9.81];

% IMU noise, ilk deneme için çok agresif yapmıyoruz.
params.sigma_g = deg2rad(0.20);
params.sigma_a = 0.080;

params.sigma_bg_rw = deg2rad(0.010);
params.sigma_ba_rw = 0.010;

% GPS 5 Hz fixed-wing log
params.sigma_gps_pos = [2.0; 2.0; 3.0];
params.R_gps_pos = diag(params.sigma_gps_pos.^2);

params.sigma_gps_vel = [0.50; 0.50; 0.80];
params.R_gps_vel = diag(params.sigma_gps_vel.^2);

% Baro offset state yok, manuel offset de yok.
% O yüzden baroyu ilk testte kapalı tutuyoruz veya gevşek bağlıyoruz.
params.sigma_baro = 6.0;
params.R_baro = params.sigma_baro^2;

params.use_joseph_form = true;

%% ============================================================
% 2) Build sim from fixed-wing MAT
% ============================================================

% Aktif uçuş kısmı.
% Bu aralık MAT dosyasının absolute log zamanına göre.
t_start = 600;
t_end   = 975;

sim = build_sim_from_fixedwing_mat(fixedwing_file, t_start, t_end);

t = sim.t;
N = numel(t);

fprintf("\n=== MAIN FIXED-WING ESKF ===\n");
fprintf("N = %d samples\n", N);
fprintf("t range = %.3f to %.3f s\n", t(1), t(end));

%% ============================================================
% 3) Initial state
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

fprintf('[fixedwing main] Initial p = [%.3f %.3f %.3f]^T m\n', state.p_n);
fprintf('[fixedwing main] Initial v = [%.3f %.3f %.3f]^T m/s\n', state.v_n);
fprintf('[fixedwing main] Initial bg = [%.6f %.6f %.6f]^T rad/s\n', state.b_g);
fprintf('[fixedwing main] Initial ba = [%.6f %.6f %.6f]^T m/s^2\n', state.b_a);
fprintf('[fixedwing main] USE_BARO = %d\n', USE_BARO);
fprintf('[fixedwing main] USE_GPS_VEL = %d\n', USE_GPS_VEL);
fprintf('[fixedwing main] USE_ATT_INIT_FOR_DEBUG = %d\n', USE_ATT_INIT_FOR_DEBUG);

%% ============================================================
% 4) Initial covariance
% ============================================================

idx_p  = 1:3;
idx_v  = 4:6;
idx_th = 7:9;
idx_bg = 10:12;
idx_ba = 13:15;

sigma_p0  = [3.0; 3.0; 5.0];
sigma_v0  = [1.0; 1.0; 1.5];

if USE_ATT_INIT_FOR_DEBUG
    sigma_th0 = deg2rad([3.0; 3.0; 8.0]);
else
    sigma_th0 = deg2rad([5.0; 5.0; 45.0]);
end

sigma_bg0 = deg2rad([0.5; 0.5; 0.8]);
sigma_ba0 = [0.30; 0.30; 0.50];

P = zeros(15,15);
P(idx_p, idx_p)   = diag(sigma_p0.^2);
P(idx_v, idx_v)   = diag(sigma_v0.^2);
P(idx_th, idx_th) = diag(sigma_th0.^2);
P(idx_bg, idx_bg) = diag(sigma_bg0.^2);
P(idx_ba, idx_ba) = diag(sigma_ba0.^2);

%% ============================================================
% 5) Allocate logs
% ============================================================

log_p = nan(3,N);
log_v = nan(3,N);
log_q = nan(4,N);
log_bg = nan(3,N);
log_ba = nan(3,N);
log_Pdiag = nan(15,N);

log_res_gps_pos = nan(3,N);
log_res_gps_vel = nan(3,N);
log_res_baro = nan(1,N);

used_gps_pos_updates = 0;
used_gps_vel_updates = 0;
used_baro_updates = 0;

log_p(:,1) = state.p_n;
log_v(:,1) = state.v_n;
log_q(:,1) = state.q_nb;
log_bg(:,1) = state.b_g;
log_ba(:,1) = state.b_a;
log_Pdiag(:,1) = diag(P);

%% ============================================================
% 6) Main ESKF replay loop
% ============================================================

fprintf("\n[fixedwing main] ESKF replay başlıyor...\n");

for k = 2:N

    dt_k = sim.t(k) - sim.t(k-1);

    if dt_k <= 0 || dt_k > 0.1
        warning("Bad dt at k=%d: %.6f s, propagation skipped.", k, dt_k);
        continue;
    end

    imu.gyro_m  = sim.imu_gyro(:,k);
    imu.accel_m = sim.imu_accel(:,k);

    %% ---------------- Propagation ----------------
    state = propagate_nominal(state, imu, params, dt_k);

    [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt_k);
    P = propagate_covariance(P, F, G, Qd, dt_k);

    %% ---------------- GPS position update ----------------
    if sim.gps_pos_available(k)
        z_gps_pos = sim.gps_pos(:,k);

        if all(isfinite(z_gps_pos))
            [state, P, residual, ~, ~] = update_gnss_pos(state, P, z_gps_pos, params);
            log_res_gps_pos(:,k) = residual;
            used_gps_pos_updates = used_gps_pos_updates + 1;
        end
    end

    %% ---------------- GPS velocity update ----------------
    if USE_GPS_VEL && sim.gps_vel_available(k)
        z_gps_vel = sim.gps_vel(:,k);

        if all(isfinite(z_gps_vel))
            [state, P, residual, ~, ~] = update_gnss_vel(state, P, z_gps_vel, params);
            log_res_gps_vel(:,k) = residual;
            used_gps_vel_updates = used_gps_vel_updates + 1;
        end
    end

    %% ---------------- Baro update ----------------
    if USE_BARO && sim.baro_available(k)
        z_baro = sim.baro(k);

        if isfinite(z_baro)
            [state, P, residual, ~, ~] = update_baro(state, P, z_baro, params);
            log_res_baro(k) = residual;
            used_baro_updates = used_baro_updates + 1;
        end
    end

    %% ---------------- Logging ----------------
    log_p(:,k) = state.p_n;
    log_v(:,k) = state.v_n;
    log_q(:,k) = state.q_nb;
    log_bg(:,k) = state.b_g;
    log_ba(:,k) = state.b_a;
    log_Pdiag(:,k) = diag(P);
end

fprintf("[fixedwing main] ESKF replay bitti.\n");

%% ============================================================
% 7) Plot observation/reference results
% ============================================================

plot_results_fixedwing_obs(sim, t, log_p, log_v, log_q, ...
    log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params);

%% ============================================================
% 8) Observation consistency metrics
% ============================================================

idx_gps_pos = sim.gps_pos_available & all(isfinite(sim.gps_pos),1);
idx_gps_vel = sim.gps_vel_available & all(isfinite(sim.gps_vel),1);
idx_baro = sim.baro_available & isfinite(sim.baro);

gps_pos_rmse = nan(3,1);
gps_pos_rmse_norm = nan;
gps_vel_rmse = nan(3,1);
gps_vel_rmse_norm = nan;
baro_rmse = nan;
baro_mean_error = nan;
baro_std_error = nan;

fprintf('\n--- FIXED-WING RUN FINISHED ---\n');

if any(idx_gps_pos)
    gps_pos_err = log_p(:,idx_gps_pos) - sim.gps_pos(:,idx_gps_pos);
    gps_pos_rmse = sqrt(mean(gps_pos_err.^2, 2, 'omitnan'));
    gps_pos_rmse_norm = sqrt(mean(sum(gps_pos_err.^2,1), 'omitnan'));

    fprintf('GPS position consistency RMSE N/E/D [m] : [%.3f %.3f %.3f]\n', gps_pos_rmse);
    fprintf('GPS position consistency RMSE norm [m]  : %.3f\n', gps_pos_rmse_norm);
end

if any(idx_gps_vel)
    gps_vel_err = log_v(:,idx_gps_vel) - sim.gps_vel(:,idx_gps_vel);
    gps_vel_rmse = sqrt(mean(gps_vel_err.^2, 2, 'omitnan'));
    gps_vel_rmse_norm = sqrt(mean(sum(gps_vel_err.^2,1), 'omitnan'));

    fprintf('GPS velocity consistency RMSE N/E/D [m/s] : [%.3f %.3f %.3f]\n', gps_vel_rmse);
    fprintf('GPS velocity consistency RMSE norm [m/s]  : %.3f\n', gps_vel_rmse_norm);
end

if any(idx_baro)
    baro_err = log_p(3,idx_baro) - sim.baro(idx_baro);

    baro_rmse = sqrt(mean(baro_err.^2, 'omitnan'));
    baro_mean_error = mean(baro_err, 'omitnan');
    baro_std_error = std(baro_err, 0, 'omitnan');

    fprintf('Baro Down consistency RMSE [m]       : %.3f\n', baro_rmse);
    fprintf('Baro Down consistency mean error [m] : %.3f\n', baro_mean_error);
    fprintf('Baro Down consistency std error [m]  : %.3f\n', baro_std_error);
end

fprintf('\n--- UPDATE COUNTS ---\n');
fprintf('GPS pos updates available : %d\n', sum(sim.gps_pos_available));
fprintf('GPS pos updates used      : %d\n', used_gps_pos_updates);
fprintf('GPS vel updates available : %d\n', sum(sim.gps_vel_available));
fprintf('GPS vel updates used      : %d\n', used_gps_vel_updates);
fprintf('BARO updates available    : %d\n', sum(sim.baro_available & isfinite(sim.baro)));
fprintf('BARO updates used         : %d\n', used_baro_updates);

fprintf('\n--- FINAL STATES ---\n');
fprintf('Final estimated p NED [m]     : [%.4f %.4f %.4f]^T\n', log_p(:,end));
fprintf('Final estimated v NED [m/s]   : [%.4f %.4f %.4f]^T\n', log_v(:,end));
fprintf('Final gyro bias [rad/s]       : [%.6f %.6f %.6f]^T\n', log_bg(:,end));
fprintf('Final accel bias [m/s^2]      : [%.6f %.6f %.6f]^T\n', log_ba(:,end));

fprintf('\n--- CONFIG USED ---\n');
fprintf('Sensor profile : %s\n', params.sensor_profile);
fprintf('USE_BARO       : %d\n', USE_BARO);
fprintf('USE_GPS_VEL    : %d\n', USE_GPS_VEL);
fprintf('ATT DEBUG INIT : %d\n', USE_ATT_INIT_FOR_DEBUG);
fprintf('sigma_gps_pos  : [%.3f %.3f %.3f]\n', params.sigma_gps_pos);
fprintf('sigma_gps_vel  : [%.3f %.3f %.3f]\n', params.sigma_gps_vel);
fprintf('sigma_baro     : %.3f\n', params.sigma_baro);

%% ============================================================
% 9) PX4 reference comparison from combined CSV
% ============================================================

px4_cmp = struct();
px4_cmp.available = false;

if USE_COMBINED_PX4_REFERENCE
    px4_cmp = compare_with_px4_reference_from_csv_autoalign( ...
        combined_csv_file, t, log_p, log_v, sim);
end

%% ============================================================
% 10) Save run output
% ============================================================

if SAVE_RUN_OUTPUT
    fprintf("\n[fixedwing main] Saving run output to %s ...\n", save_file);

    save(save_file, ...
        "fixedwing_file", "combined_csv_file", ...
        "sim", "t", ...
        "log_p", "log_v", "log_q", "log_bg", "log_ba", "log_Pdiag", ...
        "log_res_gps_pos", "log_res_gps_vel", "log_res_baro", ...
        "gps_pos_rmse", "gps_pos_rmse_norm", ...
        "gps_vel_rmse", "gps_vel_rmse_norm", ...
        "baro_rmse", "baro_mean_error", "baro_std_error", ...
        "used_gps_pos_updates", "used_gps_vel_updates", "used_baro_updates", ...
        "USE_BARO", "USE_GPS_VEL", "USE_ATT_INIT_FOR_DEBUG", ...
        "USE_COMBINED_PX4_REFERENCE", ...
        "params", "px4_cmp", "-v7.3");

    fprintf("[fixedwing main] Saved: %s\n", save_file);
end

%% ============================================================
% LOCAL FUNCTION: PX4 reference comparison with GPS auto-alignment
% ============================================================

function px4_cmp = compare_with_px4_reference_from_csv_autoalign( ...
    combined_csv_file, t, log_p, log_v, sim)

    px4_cmp = struct();
    px4_cmp.available = false;
    px4_cmp.file = combined_csv_file;
    px4_cmp.t = t;

    if ~isfile(combined_csv_file)
        warning("Combined CSV file not found: %s", combined_csv_file);
        return;
    end

    Tref = readtable(combined_csv_file);

    required_cols = ["t_s", ...
        "gps_available", "gps_n_m", "gps_e_m", "gps_d_m", ...
        "px4_local_n_m", "px4_local_e_m", "px4_local_d_m", ...
        "px4_local_vn_m_s", "px4_local_ve_m_s", "px4_local_vd_m_s"];

    if ~all(ismember(required_cols, string(Tref.Properties.VariableNames)))
        warning("Combined CSV does not contain required reference columns.");
        disp(Tref.Properties.VariableNames');
        return;
    end

    %% CSV arrays
    t_csv = Tref.t_s(:)';
    t_csv = t_csv - t_csv(1);

    gps_csv = [
        Tref.gps_n_m(:)'
        Tref.gps_e_m(:)'
        Tref.gps_d_m(:)'
    ];

    gps_csv_valid = ...
        Tref.gps_available(:)' == 1 & ...
        all(isfinite(gps_csv), 1) & ...
        isfinite(t_csv);

    px4_p_ref = [
        Tref.px4_local_n_m(:)'
        Tref.px4_local_e_m(:)'
        Tref.px4_local_d_m(:)'
    ];

    px4_v_ref = [
        Tref.px4_local_vn_m_s(:)'
        Tref.px4_local_ve_m_s(:)'
        Tref.px4_local_vd_m_s(:)'
    ];

    valid_px4_p = all(isfinite(px4_p_ref), 1) & isfinite(t_csv);
    valid_px4_v = all(isfinite(px4_v_ref), 1) & isfinite(t_csv);

    %% sim GPS target
    idx_sim_gps = sim.gps_pos_available & all(isfinite(sim.gps_pos), 1);

    if sum(idx_sim_gps) < 10 || sum(gps_csv_valid) < 10
        warning("Not enough GPS observations for auto-alignment.");
        return;
    end

    t_sim_gps = t(idx_sim_gps);
    gps_sim = sim.gps_pos(:, idx_sim_gps);

    t_csv_gps = t_csv(gps_csv_valid);
    gps_csv_data = gps_csv(:, gps_csv_valid);

    fprintf("\n=== PX4/GPS AUTO ALIGNMENT ===\n");
    fprintf("ESKF/sim t range : %.3f - %.3f s\n", t(1), t(end));
    fprintf("CSV t range      : %.3f - %.3f s\n", t_csv(1), t_csv(end));
    fprintf("sim GPS points   : %d\n", sum(idx_sim_gps));
    fprintf("CSV GPS points   : %d\n", sum(gps_csv_valid));

    %% Search tau robustly using overlap, not strict full containment
    max_tau = max(0, t_csv(end) - t(1));
    min_tau = -max(0, t(end) - t_csv(1));

    % Bu durumda genelde tau = 0 doğru. Yine de +/- 30 s arıyoruz.
    tau_grid = -30:0.5:30;

    best_tau = 0;
    best_score = inf;
    best_offset_gps = [0;0;0];
    best_valid_count = 0;

    for ii = 1:numel(tau_grid)
        tau = tau_grid(ii);

        tq = t_sim_gps + tau;

        gps_csv_i = nan(3, numel(tq));
        for ax = 1:3
            gps_csv_i(ax,:) = interp1(t_csv_gps, gps_csv_data(ax,:), tq, "linear", nan);
        end

        valid_i = all(isfinite(gps_csv_i), 1) & all(isfinite(gps_sim), 1);

        % En az GPS noktalarının %60'ı ortak olsun.
        if sum(valid_i) < 0.60*numel(t_sim_gps)
            continue;
        end

        offset_i = median(gps_sim(:,valid_i) - gps_csv_i(:,valid_i), 2, "omitnan");
        err_i = gps_sim(:,valid_i) - (gps_csv_i(:,valid_i) + offset_i);

        score_i = sqrt(mean(sum(err_i(1:2,:).^2, 1), "omitnan"));

        if score_i < best_score
            best_score = score_i;
            best_tau = tau;
            best_offset_gps = offset_i;
            best_valid_count = sum(valid_i);
        end
    end

    %% Fine search
    tau_fine = best_tau-2:0.05:best_tau+2;

    for ii = 1:numel(tau_fine)
        tau = tau_fine(ii);

        tq = t_sim_gps + tau;

        gps_csv_i = nan(3, numel(tq));
        for ax = 1:3
            gps_csv_i(ax,:) = interp1(t_csv_gps, gps_csv_data(ax,:), tq, "linear", nan);
        end

        valid_i = all(isfinite(gps_csv_i), 1) & all(isfinite(gps_sim), 1);

        if sum(valid_i) < 0.60*numel(t_sim_gps)
            continue;
        end

        offset_i = median(gps_sim(:,valid_i) - gps_csv_i(:,valid_i), 2, "omitnan");
        err_i = gps_sim(:,valid_i) - (gps_csv_i(:,valid_i) + offset_i);

        score_i = sqrt(mean(sum(err_i(1:2,:).^2, 1), "omitnan"));

        if score_i < best_score
            best_score = score_i;
            best_tau = tau;
            best_offset_gps = offset_i;
            best_valid_count = sum(valid_i);
        end
    end

    %% If search failed, force tau = 0 over common overlap
    if ~isfinite(best_score)
        warning("Auto alignment failed. Forcing tau = 0 over common overlap.");

        best_tau = 0;
        tq = t_sim_gps;

        gps_csv_i = nan(3, numel(tq));
        for ax = 1:3
            gps_csv_i(ax,:) = interp1(t_csv_gps, gps_csv_data(ax,:), tq, "linear", nan);
        end

        valid_i = all(isfinite(gps_csv_i), 1) & all(isfinite(gps_sim), 1);

        if sum(valid_i) >= 10
            best_offset_gps = median(gps_sim(:,valid_i) - gps_csv_i(:,valid_i), 2, "omitnan");
            err_i = gps_sim(:,valid_i) - (gps_csv_i(:,valid_i) + best_offset_gps);
            best_score = sqrt(mean(sum(err_i(1:2,:).^2, 1), "omitnan"));
            best_valid_count = sum(valid_i);
        end
    end

    %% Build best aligned GPS segment for plotting
    tq_best = t_sim_gps + best_tau;
    gps_csv_best = nan(3, numel(tq_best));

    for ax = 1:3
        gps_csv_best(ax,:) = interp1(t_csv_gps, gps_csv_data(ax,:), tq_best, "linear", nan);
    end

    valid_best_gps = all(isfinite(gps_csv_best),1) & all(isfinite(gps_sim),1);

    fprintf("Best CSV time offset tau [s]     : %.3f\n", best_tau);
    fprintf("GPS alignment score NE RMSE [m]  : %.3f\n", best_score);
    fprintf("GPS alignment offset N/E/D [m]   : [%.3f %.3f %.3f]\n", best_offset_gps);
    fprintf("GPS alignment valid points       : %d\n", best_valid_count);

    px4_cmp.best_tau = best_tau;
    px4_cmp.gps_alignment_score_ne_rmse = best_score;
    px4_cmp.gps_alignment_offset = best_offset_gps;
    px4_cmp.gps_alignment_valid_points = best_valid_count;

    %% Interpolate PX4 reference at aligned time
    tq_all = t + best_tau;

    px4_p_i = nan(3, numel(t));
    px4_v_i = nan(3, numel(t));

    if any(valid_px4_p)
        for ax = 1:3
            px4_p_i(ax,:) = interp1(t_csv(valid_px4_p), px4_p_ref(ax,valid_px4_p), tq_all, "linear", nan);
        end
    end

    if any(valid_px4_v)
        for ax = 1:3
            px4_v_i(ax,:) = interp1(t_csv(valid_px4_v), px4_v_ref(ax,valid_px4_v), tq_all, "linear", nan);
        end
    end

    %% Position origin alignment
    valid_p_i = all(isfinite(px4_p_i), 1) & all(isfinite(log_p), 1);

    if any(valid_p_i)
        first_valid = find(valid_p_i, 1, "first");

        offset_p_first = log_p(:,first_valid) - px4_p_i(:,first_valid);

        % Daha sağlam hizalama: GPS alignment offset'ini de dene.
        px4_p_i_first_aligned = px4_p_i + offset_p_first;
        px4_p_i_gps_aligned   = px4_p_i + best_offset_gps;

        err_first = log_p - px4_p_i_first_aligned;
        err_gps   = log_p - px4_p_i_gps_aligned;

        valid_err_first = all(isfinite(err_first),1);
        valid_err_gps   = all(isfinite(err_gps),1);

        rmse_first = sqrt(mean(sum(err_first(:,valid_err_first).^2,1), "omitnan"));
        rmse_gps   = sqrt(mean(sum(err_gps(:,valid_err_gps).^2,1), "omitnan"));

        if rmse_gps < rmse_first
            px4_p_i_aligned = px4_p_i_gps_aligned;
            offset_p = best_offset_gps;
            err_p_px4 = err_gps;
            valid_err_p = valid_err_gps;
            align_mode = "gps-median-offset";
        else
            px4_p_i_aligned = px4_p_i_first_aligned;
            offset_p = offset_p_first;
            err_p_px4 = err_first;
            valid_err_p = valid_err_first;
            align_mode = "first-sample-offset";
        end

        err_p_px4_raw = log_p - px4_p_i;

        px4_cmp.available = true;
        px4_cmp.position_alignment_mode = align_mode;
        px4_cmp.p_ref_raw = px4_p_i;
        px4_cmp.p_ref = px4_p_i_aligned;
        px4_cmp.position_alignment_offset = offset_p;

        px4_cmp.err_p_raw = err_p_px4_raw;
        px4_cmp.err_p = err_p_px4;

        px4_cmp.pos_rmse_axis_raw = sqrt(mean(err_p_px4_raw(:,valid_err_p).^2, 2, "omitnan"));
        px4_cmp.pos_rmse_norm_raw = sqrt(mean(sum(err_p_px4_raw(:,valid_err_p).^2, 1), "omitnan"));

        px4_cmp.pos_rmse_axis = sqrt(mean(err_p_px4(:,valid_err_p).^2, 2, "omitnan"));
        px4_cmp.pos_rmse_norm = sqrt(mean(sum(err_p_px4(:,valid_err_p).^2, 1), "omitnan"));
        px4_cmp.pos_mean_axis = mean(err_p_px4(:,valid_err_p), 2, "omitnan");
        px4_cmp.pos_std_axis = std(err_p_px4(:,valid_err_p), 0, 2, "omitnan");

        fprintf("\n=== ESKF vs PX4 LOCAL POSITION REFERENCE ===\n");
        fprintf("Position alignment mode          : %s\n", align_mode);
        fprintf("Raw Position RMSE N/E/D [m]      : [%.3f %.3f %.3f]\n", px4_cmp.pos_rmse_axis_raw);
        fprintf("Raw Position RMSE norm [m]       : %.3f\n", px4_cmp.pos_rmse_norm_raw);
        fprintf("Applied origin offset N/E/D [m]  : [%.3f %.3f %.3f]\n", px4_cmp.position_alignment_offset);
        fprintf("Aligned Position RMSE N/E/D [m]  : [%.3f %.3f %.3f]\n", px4_cmp.pos_rmse_axis);
        fprintf("Aligned Position RMSE norm [m]   : %.3f\n", px4_cmp.pos_rmse_norm);
        fprintf("Aligned mean error N/E/D [m]     : [%.3f %.3f %.3f]\n", px4_cmp.pos_mean_axis);

        figure('Name','PX4 Auto Alignment Check');
        plot(gps_csv(2,gps_csv_valid), gps_csv(1,gps_csv_valid), "Color", [0.7 0.7 0.7]);
        hold on;
        if any(valid_best_gps)
            plot(gps_csv_best(2,valid_best_gps)+best_offset_gps(2), ...
                 gps_csv_best(1,valid_best_gps)+best_offset_gps(1), "b.");
        end
        plot(gps_sim(2,:), gps_sim(1,:), "r.");
        grid on;
        axis equal;
        xlabel("East [m]");
        ylabel("North [m]");
        legend("Full CSV GPS", "Aligned CSV GPS segment", "sim GPS segment");
        title(sprintf("GPS Auto Alignment, tau = %.2f s", best_tau));

        figure('Name','ESKF Position Error w.r.t. PX4 EKF Reference');
        plot(t(valid_err_p), err_p_px4(1,valid_err_p), "r"); hold on;
        plot(t(valid_err_p), err_p_px4(2,valid_err_p), "g");
        plot(t(valid_err_p), err_p_px4(3,valid_err_p), "b");
        grid on;
        xlabel("Time [s]");
        ylabel("Position error wrt aligned PX4 ref [m]");
        legend("N","E","D");
        title("ESKF Position Error w.r.t. Origin-Aligned PX4 EKF Reference");

        figure('Name','Ground Track: ESKF vs PX4 Reference');
        plot(px4_p_i(2,:), px4_p_i(1,:), "Color", [0.5 0.5 0.5], "LineStyle", "--", "LineWidth", 1.0);
        hold on;
        plot(px4_p_i_aligned(2,:), px4_p_i_aligned(1,:), "k--", "LineWidth", 1.2);
        plot(log_p(2,:), log_p(1,:), "r", "LineWidth", 1.2);
        grid on;
        axis equal;
        xlabel("East [m]");
        ylabel("North [m]");
        legend("PX4 raw ref", "PX4 origin-aligned ref", "Our ESKF");
        title("Ground Track: ESKF vs PX4 Reference");

        figure('Name','Altitude: ESKF vs PX4 Reference');
        plot(t, -px4_p_i_aligned(3,:), "k--", "LineWidth", 1.2);
        hold on;
        plot(t, -log_p(3,:), "r", "LineWidth", 1.2);
        grid on;
        xlabel("Time [s]");
        ylabel("Altitude Up [m]");
        legend("PX4 EKF origin-aligned reference", "Our ESKF");
        title("Altitude: ESKF vs PX4 Reference");
    else
        warning("PX4 local position reference could not be interpolated.");
    end

    %% Velocity comparison
    valid_v_i = all(isfinite(px4_v_i), 1) & all(isfinite(log_v), 1);

    if any(valid_v_i)
        err_v_px4 = log_v - px4_v_i;

        px4_cmp.available = true;
        px4_cmp.v_ref = px4_v_i;
        px4_cmp.err_v = err_v_px4;

        px4_cmp.vel_rmse_axis = sqrt(mean(err_v_px4(:,valid_v_i).^2, 2, "omitnan"));
        px4_cmp.vel_rmse_norm = sqrt(mean(sum(err_v_px4(:,valid_v_i).^2, 1), "omitnan"));
        px4_cmp.vel_mean_axis = mean(err_v_px4(:,valid_v_i), 2, "omitnan");
        px4_cmp.vel_std_axis = std(err_v_px4(:,valid_v_i), 0, 2, "omitnan");

        fprintf("\n=== ESKF vs PX4 LOCAL VELOCITY REFERENCE ===\n");
        fprintf("Velocity RMSE N/E/D [m/s] : [%.3f %.3f %.3f]\n", px4_cmp.vel_rmse_axis);
        fprintf("Velocity RMSE norm [m/s]  : %.3f\n", px4_cmp.vel_rmse_norm);
        fprintf("Velocity mean error [m/s] : [%.3f %.3f %.3f]\n", px4_cmp.vel_mean_axis);

        figure('Name','ESKF Velocity Error w.r.t. PX4 EKF Reference');
        plot(t(valid_v_i), err_v_px4(1,valid_v_i), "r"); hold on;
        plot(t(valid_v_i), err_v_px4(2,valid_v_i), "g");
        plot(t(valid_v_i), err_v_px4(3,valid_v_i), "b");
        grid on;
        xlabel("Time [s]");
        ylabel("Velocity error wrt PX4 ref [m/s]");
        legend("V_N","V_E","V_D");
        title("ESKF Velocity Error w.r.t. PX4 EKF Reference");

        figure('Name','Velocity: ESKF vs PX4 Reference');

        subplot(3,1,1);
        plot(t(valid_v_i), px4_v_i(1,valid_v_i), "k--"); hold on;
        plot(t(valid_v_i), log_v(1,valid_v_i), "r");
        grid on;
        ylabel("V_N [m/s]");
        legend("PX4 ref", "ESKF");

        subplot(3,1,2);
        plot(t(valid_v_i), px4_v_i(2,valid_v_i), "k--"); hold on;
        plot(t(valid_v_i), log_v(2,valid_v_i), "r");
        grid on;
        ylabel("V_E [m/s]");
        legend("PX4 ref", "ESKF");

        subplot(3,1,3);
        plot(t(valid_v_i), px4_v_i(3,valid_v_i), "k--"); hold on;
        plot(t(valid_v_i), log_v(3,valid_v_i), "r");
        grid on;
        xlabel("Time [s]");
        ylabel("V_D [m/s]");
        legend("PX4 ref", "ESKF");

        sgtitle("Velocity: ESKF vs PX4 Reference");
    else
        warning("PX4 local velocity reference could not be interpolated.");
    end
end