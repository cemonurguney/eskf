function plot_results_fixedwing_obs(sim, t, log_p, log_v, log_q, ...
    log_bg, log_ba, log_Pdiag, ...
    log_res_gps_pos, log_res_gps_vel, log_res_baro, params)
%PLOT_RESULTS_FIXEDWING_OBS
% Fixed-wing gerçek ULog replay için observation/reference merkezli plotlar.
%
% Frame:
%   NED
%   p = [North; East; Down]
%   v = [Vn; Ve; Vd]
%
% Not:
%   Burada "truth" yoktur.
%   GPS, baro, logged attitude, airspeed ve wind estimate gözlem/reference
%   olarak çizilir. Bunlar absolute ground truth değildir.

    N = numel(t);

    %% Common indices
    idx_gps_pos = false(1,N);
    if isfield(sim, "gps_pos_available")
        idx_gps_pos = sim.gps_pos_available & all(isfinite(sim.gps_pos),1);
    end

    idx_gps_vel = false(1,N);
    if isfield(sim, "gps_vel_available")
        idx_gps_vel = sim.gps_vel_available & all(isfinite(sim.gps_vel),1);
    end

    idx_baro = false(1,N);
    if isfield(sim, "baro_available")
        idx_baro = sim.baro_available & isfinite(sim.baro);
    end

    %% ============================================================
    % 1) Position NED with GPS observations
    % ============================================================
    figure('Name','Fixed-Wing Position NED');
    labels_p = {'North N [m]', 'East E [m]', 'Down D [m]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, log_p(i,:), 'r', 'LineWidth', 1.2);
        hold on;

        if any(idx_gps_pos)
            plot(t(idx_gps_pos), sim.gps_pos(i,idx_gps_pos), 'k.', 'MarkerSize', 5);
            legend('ESKF estimate','GPS observation','Location','best');
        else
            legend('ESKF estimate','Location','best');
        end

        grid on;
        xlabel('Time [s]');
        ylabel(labels_p{i});
    end
    sgtitle('Fixed-Wing Position: ESKF vs GPS Observations');

    %% ============================================================
    % 2) Altitude Up view
    % ============================================================
    figure('Name','Fixed-Wing Altitude Observations');

    plot(t, -log_p(3,:), 'r', 'LineWidth', 1.3);
    hold on;

    if any(idx_gps_pos)
        plot(t(idx_gps_pos), -sim.gps_pos(3,idx_gps_pos), 'k.', 'MarkerSize', 5);
    end

    if any(idx_baro)
        plot(t(idx_baro), -sim.baro(idx_baro), 'm.', 'MarkerSize', 4);
    end

    grid on;
    xlabel('Time [s]');
    ylabel('Altitude Up = -Down [m]');
    title('Altitude: ESKF vs GPS and Baro Observations');
    legend('ESKF altitude','GPS altitude','Baro altitude','Location','best');

    %% ============================================================
    % 3) Velocity NED with GPS velocity observations
    % ============================================================
    figure('Name','Fixed-Wing Velocity NED');
    labels_v = {'V_N [m/s]', 'V_E [m/s]', 'V_D [m/s]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, log_v(i,:), 'r', 'LineWidth', 1.2);
        hold on;

        if any(idx_gps_vel)
            plot(t(idx_gps_vel), sim.gps_vel(i,idx_gps_vel), 'k.', 'MarkerSize', 5);
            legend('ESKF estimate','GPS velocity obs','Location','best');
        else
            legend('ESKF estimate','Location','best');
        end

        grid on;
        xlabel('Time [s]');
        ylabel(labels_v{i});
    end
    sgtitle('Fixed-Wing Velocity: ESKF vs GPS Velocity Observations');

    %% ============================================================
    % 4) Speed magnitude
    % ============================================================
    figure('Name','Speed Magnitude');

    est_speed = vecnorm(log_v,2,1);
    plot(t, est_speed, 'r', 'LineWidth', 1.3);
    hold on;

    leg = {'ESKF ground speed'};

    if any(idx_gps_vel)
        gps_speed = vecnorm(sim.gps_vel(:,idx_gps_vel),2,1);
        plot(t(idx_gps_vel), gps_speed, 'k.', 'MarkerSize', 5);
        leg{end+1} = 'GPS ground speed';
    end

    if isfield(sim, "airspeed_t") && ~isempty(sim.airspeed_t)
        if isfield(sim, "ias_m_s") && ~isempty(sim.ias_m_s)
            plot(sim.airspeed_t, sim.ias_m_s, 'b.', 'MarkerSize', 4);
            leg{end+1} = 'IAS';
        end
        if isfield(sim, "tas_m_s") && ~isempty(sim.tas_m_s)
            plot(sim.airspeed_t, sim.tas_m_s, 'g.', 'MarkerSize', 4);
            leg{end+1} = 'TAS';
        end
    end

    grid on;
    xlabel('Time [s]');
    ylabel('Speed [m/s]');
    title('Speed Observations');
    legend(leg,'Location','best');

    %% ============================================================
    % 5) Attitude: ESKF vs logged attitude reference
    % ============================================================
    eul_est = nan(3,N);
    for k = 1:N
        if all(isfinite(log_q(:,k))) && norm(log_q(:,k)) > 0
            eul_est(:,k) = quat_to_euler321_local(log_q(:,k));
        end
    end

    eul_est = unwrap(eul_est, [], 2);
    eul_est_deg = rad2deg(eul_est);

    figure('Name','Fixed-Wing Attitude');
    labels_att = {'roll [deg]', 'pitch [deg]', 'yaw [deg]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, eul_est_deg(i,:), 'r', 'LineWidth', 1.2);
        hold on;

        if isfield(sim, "att_t") && ~isempty(sim.att_t) && isfield(sim, "rpy_ref") && ~isempty(sim.rpy_ref)
            rpy_ref_deg = rad2deg(unwrap(sim.rpy_ref, [], 2));
            plot(sim.att_t, rpy_ref_deg(i,:), 'k.', 'MarkerSize', 4);
            legend('ESKF attitude','Logged attitude ref','Location','best');
        else
            legend('ESKF attitude','Location','best');
        end

        grid on;
        xlabel('Time [s]');
        ylabel(labels_att{i});
    end
    sgtitle('Attitude: ESKF vs Logged Attitude Reference');

    %% ============================================================
    % 6) Ground track NED
    % ============================================================
    figure('Name','Ground Track NE');

    if any(idx_gps_pos)
        plot(sim.gps_pos(2,idx_gps_pos), sim.gps_pos(1,idx_gps_pos), 'k.', 'MarkerSize', 5);
        hold on;
    end

    plot(log_p(2,:), log_p(1,:), 'r', 'LineWidth', 1.3);

    if any(idx_gps_pos)
        plot(sim.gps_pos(2,find(idx_gps_pos,1,'first')), ...
             sim.gps_pos(1,find(idx_gps_pos,1,'first')), ...
             'go', 'MarkerSize', 8, 'LineWidth', 1.5);

        plot(sim.gps_pos(2,find(idx_gps_pos,1,'last')), ...
             sim.gps_pos(1,find(idx_gps_pos,1,'last')), ...
             'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
    end

    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    title('Ground Track: GPS Observations vs ESKF');
    legend('GPS observation','ESKF estimate','start','end','Location','best');

    %% ============================================================
    % 7) 3D trajectory with altitude-up display
    % ============================================================
    figure('Name','3D Trajectory NEU');

    if any(idx_gps_pos)
        plot3(sim.gps_pos(2,idx_gps_pos), ...
              sim.gps_pos(1,idx_gps_pos), ...
             -sim.gps_pos(3,idx_gps_pos), ...
              'k.', 'MarkerSize', 5);
        hold on;
    end

    plot3(log_p(2,:), log_p(1,:), -log_p(3,:), 'r', 'LineWidth', 1.4);

    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    zlabel('Altitude Up [m]');
    title('3D Trajectory: Observation View');
    legend('GPS observation','ESKF estimate','Location','best');
    view(3);

    %% ============================================================
    % 8) IMU observations
    % ============================================================
    if isfield(sim, "imu_gyro") && isfield(sim, "imu_accel")
        figure('Name','IMU Observations');

        subplot(3,1,1);
        plot(t, sim.imu_gyro.', 'LineWidth', 0.8);
        grid on;
        xlabel('Time [s]');
        ylabel('Gyro [rad/s]');
        legend('\omega_x','\omega_y','\omega_z','Location','best');
        title('Gyroscope Measurements');

        subplot(3,1,2);
        plot(t, sim.imu_accel.', 'LineWidth', 0.8);
        grid on;
        xlabel('Time [s]');
        ylabel('Accel [m/s^2]');
        legend('a_x','a_y','a_z','Location','best');
        title('Accelerometer Measurements');

        subplot(3,1,3);
        plot(t, vecnorm(sim.imu_accel,2,1), 'k', 'LineWidth', 1.0);
        grid on;
        xlabel('Time [s]');
        ylabel('||a|| [m/s^2]');
        title('Accelerometer Norm');
    end

    %% ============================================================
    % 9) Wind estimate observation, if available
    % ============================================================
    if isfield(sim, "wind_t") && ~isempty(sim.wind_t) && isfield(sim, "wind_ned_m_s") && ~isempty(sim.wind_ned_m_s)
        figure('Name','PX4 Wind Estimate Reference');

        labels_w = {'Wind N [m/s]', 'Wind E [m/s]', 'Wind D [m/s]'};

        for i = 1:min(3,size(sim.wind_ned_m_s,1))
            subplot(3,1,i);
            plot(sim.wind_t, sim.wind_ned_m_s(i,:), 'LineWidth', 1.1);
            grid on;
            xlabel('Time [s]');
            ylabel(labels_w{i});
        end

        sgtitle('Logged Wind Estimate Reference');
    end

    %% ============================================================
    % 10) Bias estimates
    % ============================================================
    figure('Name','Gyro Bias Estimate');
    labels_bg = {'b_{g,x} [rad/s]', 'b_{g,y} [rad/s]', 'b_{g,z} [rad/s]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, log_bg(i,:), 'r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_bg{i});
    end
    sgtitle('Estimated Gyro Bias');

    figure('Name','Accel Bias Estimate');
    labels_ba = {'b_{a,x} [m/s^2]', 'b_{a,y} [m/s^2]', 'b_{a,z} [m/s^2]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, log_ba(i,:), 'r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_ba{i});
    end
    sgtitle('Estimated Accelerometer Bias');

    %% ============================================================
    % 11) Covariance diagonal, grouped
    % ============================================================
    figure('Name','Covariance Diagonal Groups');

    subplot(5,1,1);
    plot(t, log_Pdiag(1:3,:).', 'LineWidth', 1.0);
    grid on;
    ylabel('P_p');
    legend('N','E','D','Location','best');
    title('Position covariance diagonal');

    subplot(5,1,2);
    plot(t, log_Pdiag(4:6,:).', 'LineWidth', 1.0);
    grid on;
    ylabel('P_v');
    legend('V_N','V_E','V_D','Location','best');
    title('Velocity covariance diagonal');

    subplot(5,1,3);
    plot(t, log_Pdiag(7:9,:).', 'LineWidth', 1.0);
    grid on;
    ylabel('P_\theta');
    legend('\theta_x','\theta_y','\theta_z','Location','best');
    title('Attitude error covariance diagonal');

    subplot(5,1,4);
    plot(t, log_Pdiag(10:12,:).', 'LineWidth', 1.0);
    grid on;
    ylabel('P_{bg}');
    legend('b_{gx}','b_{gy}','b_{gz}','Location','best');
    title('Gyro bias covariance diagonal');

    subplot(5,1,5);
    plot(t, log_Pdiag(13:15,:).', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    ylabel('P_{ba}');
    legend('b_{ax}','b_{ay}','b_{az}','Location','best');
    title('Accel bias covariance diagonal');

    %% ============================================================
    % 12) Residuals
    % ============================================================
    figure('Name','Measurement Residuals');

    subplot(3,1,1);
    plot(t, log_res_gps_pos.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('GPS pos residual [m]');
    legend('N','E','D','Location','best');
    title('GPS Position Residual');

    subplot(3,1,2);
    plot(t, log_res_gps_vel.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('GPS vel residual [m/s]');
    legend('V_N','V_E','V_D','Location','best');
    title('GPS Velocity Residual');

    subplot(3,1,3);
    plot(t, log_res_baro.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('Baro residual [m]');
    legend('baro','Location','best');
    title('Barometer Residual');

    %% ============================================================
    % 13) Consistency errors w.r.t observations
    % ============================================================
    figure('Name','Observation Consistency Errors');

    subplot(3,1,1);
    if any(idx_gps_pos)
        gps_pos_err = log_p(:,idx_gps_pos) - sim.gps_pos(:,idx_gps_pos);
        plot(t(idx_gps_pos), gps_pos_err.', '.');
        grid on;
        ylabel('p - GPS [m]');
        legend('N','E','D','Location','best');
        title(sprintf('GPS Position Consistency, RMSE norm = %.3f m', ...
            sqrt(mean(sum(gps_pos_err.^2,1), 'omitnan'))));
    end

    subplot(3,1,2);
    if any(idx_gps_vel)
        gps_vel_err = log_v(:,idx_gps_vel) - sim.gps_vel(:,idx_gps_vel);
        plot(t(idx_gps_vel), gps_vel_err.', '.');
        grid on;
        ylabel('v - GPS vel [m/s]');
        legend('V_N','V_E','V_D','Location','best');
        title(sprintf('GPS Velocity Consistency, RMSE norm = %.3f m/s', ...
            sqrt(mean(sum(gps_vel_err.^2,1), 'omitnan'))));
    end

    subplot(3,1,3);
    if any(idx_baro)
        baro_err = log_p(3,idx_baro) - sim.baro(idx_baro);
        plot(t(idx_baro), baro_err, 'k.');
        grid on;
        xlabel('Time [s]');
        ylabel('D - baro D [m]');
        title(sprintf('Baro Down Consistency, RMSE = %.3f m, mean = %.3f m', ...
            sqrt(mean(baro_err.^2, 'omitnan')), mean(baro_err, 'omitnan')));
    end

    %% ============================================================
    % 14) Print numeric observation summary
    % ============================================================
    fprintf('\n--- FIXED-WING OBSERVATION CONSISTENCY SUMMARY ---\n');

    if any(idx_gps_pos)
        gps_pos_err = log_p(:,idx_gps_pos) - sim.gps_pos(:,idx_gps_pos);
        gps_pos_rmse_axis = sqrt(mean(gps_pos_err.^2, 2, 'omitnan'));
        fprintf('GPS pos RMSE N/E/D [m]      : [%.3f %.3f %.3f]\n', gps_pos_rmse_axis);
        fprintf('GPS pos RMSE norm [m]       : %.3f\n', sqrt(mean(sum(gps_pos_err.^2,1), 'omitnan')));
    end

    if any(idx_gps_vel)
        gps_vel_err = log_v(:,idx_gps_vel) - sim.gps_vel(:,idx_gps_vel);
        gps_vel_rmse_axis = sqrt(mean(gps_vel_err.^2, 2, 'omitnan'));
        fprintf('GPS vel RMSE N/E/D [m/s]    : [%.3f %.3f %.3f]\n', gps_vel_rmse_axis);
        fprintf('GPS vel RMSE norm [m/s]     : %.3f\n', sqrt(mean(sum(gps_vel_err.^2,1), 'omitnan')));
    end

    if any(idx_baro)
        baro_err = log_p(3,idx_baro) - sim.baro(idx_baro);
        fprintf('Baro Down RMSE [m]          : %.3f\n', sqrt(mean(baro_err.^2, 'omitnan')));
        fprintf('Baro Down mean error [m]    : %.3f\n', mean(baro_err, 'omitnan'));
    end

    if isfield(params, "sensor_profile")
        fprintf('Sensor profile              : %s\n', params.sensor_profile);
    end
end

%% =============================================================
% LOCAL HELPER
% =============================================================
function eul = quat_to_euler321_local(q)
    q = q(:);
    q = q / norm(q);

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    roll = atan2(2*(qw*qx + qy*qz), ...
                 1 - 2*(qx^2 + qy^2));

    s = 2*(qw*qy - qz*qx);
    s = max(min(s, 1), -1);
    pitch = asin(s);

    yaw = atan2(2*(qw*qz + qx*qy), ...
                1 - 2*(qy^2 + qz^2));

    eul = [roll; pitch; yaw];
end