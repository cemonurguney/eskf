function plot_results_raw(sim, t, log_p, log_v, log_q, log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params)
%PLOT_RESULTS_RAW
% Raw MRS rectangle testi için ENU uyumlu plot fonksiyonu.

    %% 1) Position
    figure('Name','Position ENU');
    labels_p = {'East E [m]', 'North N [m]', 'Up U [m]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, sim.p_true(i,:), 'b', 'LineWidth', 1.4); hold on;
        plot(t, log_p(i,:), '--r', 'LineWidth', 1.2);

        if isfield(sim, 'gps_pos') && i <= size(sim.gps_pos,1)
            idx = sim.gps_pos_available & isfinite(sim.gps_pos(i,:));
            plot(t(idx), sim.gps_pos(i,idx), 'k.', 'MarkerSize', 6);
            legend('truth','estimate','raw GPS','Location','best');
        else
            legend('truth','estimate','Location','best');
        end

        grid on;
        xlabel('Time [s]');
        ylabel(labels_p{i});
        if i == 3 && isfield(sim, 'baro')
            idx_baro = sim.baro_available & isfinite(sim.baro);
            plot(t(idx_baro), sim.baro(idx_baro), 'm.', 'MarkerSize', 5);
            legend('truth','estimate','raw GPS','baro','Location','best');
        end
    end

    %% 2) Velocity
    figure('Name','Velocity ENU');
    labels_v = {'v_E [m/s]', 'v_N [m/s]', 'v_U [m/s]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, sim.v_true(i,:), 'b', 'LineWidth', 1.4); hold on;
        plot(t, log_v(i,:), '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_v{i});
        legend('truth','estimate','Location','best');
    end

    %% 3) Attitude
    N = numel(t);
    eul_true = zeros(3, N);
    eul_est = zeros(3, N);

    for k = 1:N
        eul_true(:,k) = quat_to_euler321_local(sim.q_true(:,k));
        eul_est(:,k) = quat_to_euler321_local(log_q(:,k));
    end

    eul_true = unwrap(eul_true, [], 2);
    eul_est = unwrap(eul_est, [], 2);

    eul_true_deg = rad2deg(eul_true);
    eul_est_deg = rad2deg(eul_est);

    figure('Name','Attitude');
    labels_att = {'roll [deg]', 'pitch [deg]', 'yaw [deg]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, eul_true_deg(i,:), 'b', 'LineWidth', 1.4); hold on;
        plot(t, eul_est_deg(i,:), '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_att{i});
        legend('truth','estimate','Location','best');
    end

    %% 4) 3D Trajectory ENU
    figure('Name','3D Trajectory ENU');

    x_true = sim.p_true(1,:);
    y_true = sim.p_true(2,:);
    z_true = sim.p_true(3,:);

    x_est = log_p(1,:);
    y_est = log_p(2,:);
    z_est = log_p(3,:);

    plot3(x_true, y_true, z_true, 'b', 'LineWidth', 1.8); hold on;
    plot3(x_est, y_est, z_est, '--r', 'LineWidth', 1.6);

    if isfield(sim, 'gps_pos')
        idx_gps = sim.gps_pos_available & all(isfinite(sim.gps_pos),1);
        plot3(sim.gps_pos(1,idx_gps), sim.gps_pos(2,idx_gps), sim.gps_pos(3,idx_gps), ...
            'k.', 'MarkerSize', 5);
    end

    plot3(x_true(1), y_true(1), z_true(1), 'go', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot3(x_true(end), y_true(end), z_true(end), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);

    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    zlabel('Up [m]');
    title('3D trajectory ENU: truth vs estimate');
    legend('truth','estimate','raw GPS','start','end','Location','best');
    view(3);

    %% 5) XY View
    figure('Name','XY trajectory ENU');
    plot(sim.p_true(1,:), sim.p_true(2,:), 'b', 'LineWidth', 1.6); hold on;
    plot(log_p(1,:), log_p(2,:), '--r', 'LineWidth', 1.4);

    if isfield(sim, 'gps_pos')
        idx_gps = sim.gps_pos_available & all(isfinite(sim.gps_pos),1);
        plot(sim.gps_pos(1,idx_gps), sim.gps_pos(2,idx_gps), 'k.', 'MarkerSize', 5);
    end

    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    title('XY trajectory ENU');
    legend('truth','estimate','raw GPS','Location','best');

    %% 6) Bias estimates
    figure('Name','Gyro Bias Estimate');
    labels_bg = {'b_{g,x} [rad/s]', 'b_{g,y} [rad/s]', 'b_{g,z} [rad/s]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, log_bg(i,:), '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_bg{i});
        legend('estimate','Location','best');
    end

    figure('Name','Accel Bias Estimate');
    labels_ba = {'b_{a,x} [m/s^2]', 'b_{a,y} [m/s^2]', 'b_{a,z} [m/s^2]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, log_ba(i,:), '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_ba{i});
        legend('estimate','Location','best');
    end

    %% 7) Covariance diagonal
    figure('Name','P diagonal');
    plot(t, log_Pdiag.', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    ylabel('diag(P)');
    title('Covariance diagonal elements');

    %% 8) Residuals
    figure('Name','Residuals');

    subplot(3,1,1);
    plot(t, log_res_gps_pos.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('GPS pos residual [m]');
    legend('E','N','U','Location','best');

    subplot(3,1,2);
    plot(t, log_res_gps_vel.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('GPS vel residual');
    legend('E','N','U','Location','best');

    subplot(3,1,3);
    plot(t, log_res_baro.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('Baro residual [m]');
    legend('baro','Location','best');

    %% 9) Position error norm
    err_p = log_p - sim.p_true;
    figure('Name','Position Error Norm');
    plot(t, vecnorm(err_p,2,1), 'LineWidth', 1.3);
    grid on;
    xlabel('Time [s]');
    ylabel('||p_{est} - p_{truth}|| [m]');
    title('Position error norm');

end

%% =============================================================
% LOCAL HELPERS
% =============================================================
function eul = quat_to_euler321_local(q)
    q = q(:);
    q = q / norm(q);

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    roll = atan2(2*(qw*qx + qy*qz), 1 - 2*(qx^2 + qy^2));

    s = 2*(qw*qy - qz*qx);
    s = max(min(s, 1), -1);
    pitch = asin(s);

    yaw = atan2(2*(qw*qz + qx*qy), 1 - 2*(qy^2 + qz^2));

    eul = [roll; pitch; yaw];
end