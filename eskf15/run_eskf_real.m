clear; clc;

% 1) DATA
sim = build_sim_from_mat("hover_clean.mat");

% 2) FILTER INIT
[state, P, params] = init_filter();

N = length(sim.t);

log_p = zeros(3,N);

% 3) MAIN LOOP
for k = 1:N

    imu.gyro_m  = sim.imu_gyro(:,k);
    imu.accel_m = sim.imu_accel(:,k);

    dt = sim.dt;

    % PROPAGATION
    state = propagate_nominal(state, imu, params, dt);
    [F,G,Qd] = compute_F_G_Qd(state, imu, params, dt);
    P = propagate_covariance(P, F, G, Qd, dt);

    % GPS
    if sim.gps_pos_available(k)
        z = sim.gps_pos(:,k);
        [state, P] = update_gnss_pos(state, P, z, params);
    end

    % BARO
    if sim.baro_available(k)
        z = sim.baro(k);
        [state, P] = update_baro(state, P, z, params);
    end

    log_p(:,k) = state.p_n;
end

% 4) ERROR
err = log_p - sim.p_true;

figure
plot(sim.t, err')
grid on
title("Position Error (REAL DATA)")
legend("x","y","z")