clear; clc; close all;

%% =========================
% MODE SEÇ
%% =========================
mode = "real";   % "real" | "sin" | "const"

%% =========================
% DATA SEÇ
%% =========================
switch mode
    case "real"
        sim = build_sim_from_mat("hover_clean.mat");

    case "sin"
        sim = build_sim_sine();

    case "const"
        sim = build_sim_const_accel();

    otherwise
        error("Bilinmeyen mode: " + mode);
end

%% =========================
% FILTER INIT
%% =========================
[state, P, params] = init_filter();

% Gerçek veri için başlangıcı truth'a yaklaştır
if mode == "real"
    state.p_n = sim.p_true(:,1);
    state.v_n = sim.v_true(:,1);
    state.q_nb = sim.q_true(:,1);
end

N = length(sim.t);

%% =========================
% LOG
%% =========================
log_p  = zeros(3,N);
log_v  = zeros(3,N);
log_q  = zeros(4,N);
log_bg = zeros(3,N);
log_ba = zeros(3,N);

log_res_gps  = nan(3,N);
log_res_baro = nan(1,N);

%% =========================
% MAIN LOOP
%% =========================
for k = 1:N

    imu.gyro_m  = sim.imu_gyro(:,k);
    imu.accel_m = sim.imu_accel(:,k);

    dt = sim.dt;

    % ---------------------------------
    % 1) PROPAGATION
    % ---------------------------------
    state = propagate_nominal(state, imu, params, dt);
    [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt);
    P = propagate_covariance(P, F, G, Qd, dt);

    % ---------------------------------
    % 2) GPS POSITION UPDATE
    % ---------------------------------
    if sim.gps_pos_available(k)
        z_gps = sim.gps_pos(:,k);

        % BEKLENEN: update fonksiyonu dx_hat döndürmeli
        [state, P, dx_hat_gps, residual_gps] = update_gnss_pos(state, P, z_gps, params);

        state = inject_error_state(state, dx_hat_gps);

        if ~isempty(residual_gps)
            log_res_gps(:,k) = residual_gps;
        end
    end

    % ---------------------------------
    % 3) BARO UPDATE
    % ---------------------------------
    if sim.baro_available(k)
        z_baro = sim.baro(k);

        % BEKLENEN: update fonksiyonu dx_hat döndürmeli
        [state, P, dx_hat_baro, residual_baro] = update_baro(state, P, z_baro, params);

        state = inject_error_state(state, dx_hat_baro);

        if ~isempty(residual_baro)
            log_res_baro(k) = residual_baro;
        end
    end

    % ---------------------------------
    % 4) LOG
    % ---------------------------------
    log_p(:,k)  = state.p_n;
    log_v(:,k)  = state.v_n;
    log_q(:,k)  = state.q_nb;
    log_bg(:,k) = state.b_g;
    log_ba(:,k) = state.b_a;
end

%% =========================
% ERROR
%% =========================
err_p = log_p - sim.p_true;
err_v = log_v - sim.v_true;

%% =========================
% PLOTS
%% =========================
figure('Name','Position Error','Color','w');
plot(sim.t, err_p');
grid on;
xlabel('Time [s]');
ylabel('Position Error [m]');
title('Position Error');
legend('x','y','z');

figure('Name','Velocity Error','Color','w');
plot(sim.t, err_v');
grid on;
xlabel('Time [s]');
ylabel('Velocity Error [m/s]');
title('Velocity Error');
legend('vx','vy','vz');

figure('Name','Estimated Position vs Truth','Color','w');
plot(sim.t, sim.p_true(1,:), '--', sim.t, log_p(1,:), '-'); hold on;
plot(sim.t, sim.p_true(2,:), '--', sim.t, log_p(2,:), '-');
plot(sim.t, sim.p_true(3,:), '--', sim.t, log_p(3,:), '-');
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
title('Estimated vs Truth Position');
legend('x true','x est','y true','y est','z true','z est');

disp('--- RUN FINISHED ---');
disp(['Final position error norm [m]: ', num2str(norm(err_p(:,end)))]);
disp(['Final velocity error norm [m/s]: ', num2str(norm(err_v(:,end)))]);