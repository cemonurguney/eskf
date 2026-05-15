function plot_results(sim, t, log_p, log_v, log_q, log_bg, log_ba, log_Pdiag, ...
                      log_res_gps_pos, log_res_gps_vel, log_res_baro, params)
%PLOT_RESULTS
% Sonuç grafiklerini çizer.

    %% 1) Position
    figure('Name','Position');
    labels_p = {'p_N [m]', 'p_E [m]', 'p_D [m]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, sim.p_true(i,:), 'b', 'LineWidth', 1.4); hold on;
        plot(t, log_p(i,:), '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_p{i});
        legend('truth','estimate','Location','best');
    end

    %% 2) Velocity
    figure('Name','Velocity');
    labels_v = {'v_N [m/s]', 'v_E [m/s]', 'v_D [m/s]'};

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
    eul_est  = zeros(3, N);

    for k = 1:N
        eul_true(:,k) = quat_to_euler321_local(sim.q_true(:,k));
        eul_est(:,k)  = quat_to_euler321_local(log_q(:,k));
    end

    eul_true = unwrap(eul_true, [], 2);
    eul_est  = unwrap(eul_est,  [], 2);

    eul_true_deg = rad2deg(eul_true);
    eul_est_deg  = rad2deg(eul_est);

    figure('Name','Attitude');
    labels_att = {'roll [deg]', 'pitch [deg]', 'yaw [deg]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, eul_true_deg(i,:), 'b', 'LineWidth', 1.4); hold on;
        plot(t, eul_est_deg(i,:),  '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_att{i});
        legend('truth','estimate','Location','best');
    end

    %% 4) 3D Trajectory
    figure('Name','3D Trajectory');

    x_true = sim.p_true(2,:);
    y_true = sim.p_true(1,:);
    z_true = -sim.p_true(3,:);

    x_est = log_p(2,:);
    y_est = log_p(1,:);
    z_est = -log_p(3,:);

    plot3(x_true, y_true, z_true, 'b', 'LineWidth', 1.8); hold on;
    plot3(x_est,  y_est,  z_est,  '--r', 'LineWidth', 1.6);

    plot3(x_true(1),   y_true(1),   z_true(1),   'go', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot3(x_true(end), y_true(end), z_true(end), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);

    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    zlabel('Up [m]');
    title('3D trajectory: truth vs estimate');
    legend('truth','estimate','start','end','Location','best');
    view(3);

    %% 5) 3D Pose Samples
    figure('Name','3D Pose Samples');

    plot3(x_true, y_true, z_true, 'b', 'LineWidth', 1.2); hold on;
    plot3(x_est,  y_est,  z_est,  '--r', 'LineWidth', 1.2);

    idx_samples = round(linspace(1, N, 12));
    idx_samples = unique(max(1, min(N, idx_samples)));

    axis_len_truth = 18;
    axis_len_est   = 14;

    for ii = 1:numel(idx_samples)
        k = idx_samples(ii);

        draw_pose_axes_local(sim.p_true(:,k), sim.q_true(:,k), axis_len_truth, '-');
        draw_pose_axes_local(log_p(:,k),      log_q(:,k),      axis_len_est,   '--');
    end

    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    zlabel('Up [m]');
    title('3D pose samples: truth (solid) vs estimate (dashed)');
    view(3);

    %% 6) Compare Views
    plot_trajectory_compare_views_local(sim, log_p, t);

    %% 7) UAV-like 3D Animation
    %animate_uav_compare_local(sim, log_p, log_q, t, params);

    %% 8) Gyro Bias Estimate
    figure('Name','Gyro Bias Estimate');
    labels_bg = {'b_{g,x} [rad/s]', 'b_{g,y} [rad/s]', 'b_{g,z} [rad/s]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, sim.bg_true(i,:), 'b', 'LineWidth', 1.4); hold on;
        plot(t, log_bg(i,:), '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_bg{i});
        legend('truth','estimate','Location','best');
    end

    %% 9) Accel Bias Estimate
    figure('Name','Accel Bias Estimate');
    labels_ba = {'b_{a,x} [m/s^2]', 'b_{a,y} [m/s^2]', 'b_{a,z} [m/s^2]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, sim.ba_true(i,:), 'b', 'LineWidth', 1.4); hold on;
        plot(t, log_ba(i,:), '--r', 'LineWidth', 1.2);
        grid on;
        xlabel('Time [s]');
        ylabel(labels_ba{i});
        legend('truth','estimate','Location','best');
    end

    %% 10) Covariance diagonal
    figure('Name','P diagonal');
    plot(t, log_Pdiag.', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    ylabel('diag(P)');
    title('Covariance diagonal elements');

    %% 11) Residuals
    figure('Name','Residuals');

    subplot(3,1,1);
    plot(t, log_res_gps_pos.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('GPS pos residual');
    legend('N','E','D','Location','best');

    subplot(3,1,2);
    plot(t, log_res_gps_vel.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('GPS vel residual');
    legend('N','E','D','Location','best');

    subplot(3,1,3);
    plot(t, log_res_baro.', '.');
    grid on;
    xlabel('Time [s]');
    ylabel('Baro residual');
    legend('baro','Location','best');
end

%% =============================================================
% YEREL YARDIMCI FONKSİYONLAR
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

function draw_pose_axes_local(p_n, q_nb, L, line_style)
    R_nb = quat_to_rotmat_local(q_nb);

    origin_plot = [p_n(2); p_n(1); -p_n(3)];

    ex_n = R_nb(:,1);
    ey_n = R_nb(:,2);
    ez_n = R_nb(:,3);

    ex_plot = [ex_n(2); ex_n(1); -ex_n(3)];
    ey_plot = [ey_n(2); ey_n(1); -ey_n(3)];
    ez_plot = [ez_n(2); ez_n(1); -ez_n(3)];

    quiver3(origin_plot(1), origin_plot(2), origin_plot(3), ...
            L*ex_plot(1), L*ex_plot(2), L*ex_plot(3), 0, ...
            'r', 'LineStyle', line_style, 'LineWidth', 1.2); hold on;

    quiver3(origin_plot(1), origin_plot(2), origin_plot(3), ...
            L*ey_plot(1), L*ey_plot(2), L*ey_plot(3), 0, ...
            'g', 'LineStyle', line_style, 'LineWidth', 1.2); hold on;

    quiver3(origin_plot(1), origin_plot(2), origin_plot(3), ...
            L*ez_plot(1), L*ez_plot(2), L*ez_plot(3), 0, ...
            'b', 'LineStyle', line_style, 'LineWidth', 1.2); hold on;
end

function plot_trajectory_compare_views_local(sim, log_p, t)
    p_true = sim.p_true;
    p_est  = log_p;

    x_true = p_true(2,:);
    y_true = p_true(1,:);
    z_true = -p_true(3,:);

    x_est = p_est(2,:);
    y_est = p_est(1,:);
    z_est = -p_est(3,:);

    err = p_est - p_true;
    err_norm = sqrt(sum(err.^2,1));

    figure('Name','Trajectory Compare Views');

    subplot(2,2,1);
    plot3(x_true, y_true, z_true, 'b', 'LineWidth', 1.8); hold on;
    plot3(x_est,  y_est,  z_est,  '--r', 'LineWidth', 1.6);
    plot3(x_true(1),   y_true(1),   z_true(1),   'go', 'MarkerSize', 7, 'LineWidth', 1.5);
    plot3(x_true(end), y_true(end), z_true(end), 'ko', 'MarkerSize', 7, 'LineWidth', 1.5);
    grid on; axis equal;
    xlabel('East [m]'); ylabel('North [m]'); zlabel('Up [m]');
    title('3D trajectory');
    legend('truth','estimate','start','end','Location','best');
    view(3);

    subplot(2,2,2);
    plot(x_true, y_true, 'b', 'LineWidth', 1.8); hold on;
    plot(x_est,  y_est,  '--r', 'LineWidth', 1.6);
    grid on; axis equal;
    xlabel('East [m]'); ylabel('North [m]');
    title('Top view (East-North)');
    legend('truth','estimate','Location','best');

    subplot(2,2,3);
    plot(x_true, z_true, 'b', 'LineWidth', 1.8); hold on;
    plot(x_est,  z_est,  '--r', 'LineWidth', 1.6);
    grid on; axis equal;
    xlabel('East [m]'); ylabel('Up [m]');
    title('Side view (East-Up)');
    legend('truth','estimate','Location','best');

    subplot(2,2,4);
    plot(t, err_norm, 'k', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('||p_{est} - p_{true}|| [m]');
    title('Position error norm');
end

% function animate_uav_compare_local(sim, log_p, log_q, t, params)
% % Animasyon süresi simülasyon süresiyle senkron akar.
% % playback_rate = 1.0 ise gerçek zamanlı oynar.
% 
%     if isfield(params,'anim') && isfield(params.anim,'step')
%         step = params.anim.step;
%     else
%         step = 5;
%     end
% 
%     if isfield(params,'anim') && isfield(params.anim,'trail_len')
%         trail_len = params.anim.trail_len;
%     else
%         trail_len = 120;
%     end
% 
%     if isfield(params,'anim') && isfield(params.anim,'axis_len_truth')
%         axis_len_truth = params.anim.axis_len_truth;
%     else
%         axis_len_truth = 14;
%     end
% 
%     if isfield(params,'anim') && isfield(params.anim,'axis_len_est')
%         axis_len_est = params.anim.axis_len_est;
%     else
%         axis_len_est = 12;
%     end
% 
%     if isfield(params,'anim') && isfield(params.anim,'playback_rate')
%         playback_rate = params.anim.playback_rate;
%     else
%         playback_rate = 1.0;
%     end
% 
%     if playback_rate <= 0
%         playback_rate = 1.0;
%     end
% 
%     p_true = sim.p_true;
%     p_est  = log_p;
% 
%     N = numel(t);
% 
%     x_true = p_true(2,:);
%     y_true = p_true(1,:);
%     z_true = -p_true(3,:);
% 
%     x_est = p_est(2,:);
%     y_est = p_est(1,:);
%     z_est = -p_est(3,:);
% 
%     figAnim = figure('Name','UAV-like 3D Animation');
%     axAnim = axes('Parent', figAnim);
% 
%     hold(axAnim, 'on');
%     grid(axAnim, 'on');
%     axis(axAnim, 'equal');
%     xlabel(axAnim, 'East [m]');
%     ylabel(axAnim, 'North [m]');
%     zlabel(axAnim, 'Up [m]');
%     view(axAnim, 3);
% 
%     xmin = min([x_true x_est]); xmax = max([x_true x_est]);
%     ymin = min([y_true y_est]); ymax = max([y_true y_est]);
%     zmin = min([z_true z_est]); zmax = max([z_true z_est]);
% 
%     dx = xmax - xmin; if dx < 1, dx = 1; end
%     dy = ymax - ymin; if dy < 1, dy = 1; end
%     dz = zmax - zmin; if dz < 1, dz = 1; end
% 
%     xlim(axAnim, [xmin-0.1*dx, xmax+0.1*dx]);
%     ylim(axAnim, [ymin-0.1*dy, ymax+0.1*dy]);
%     zlim(axAnim, [zmin-0.1*dz, zmax+0.1*dz]);
% 
%     hTrailTrue = plot3(axAnim, nan, nan, nan, 'b', 'LineWidth', 1.8);
%     hTrailEst  = plot3(axAnim, nan, nan, nan, '--r', 'LineWidth', 1.6);
% 
%     hBodyTrue = plot3(axAnim, nan, nan, nan, 'bo', 'MarkerSize', 7, 'LineWidth', 1.5);
%     hBodyEst  = plot3(axAnim, nan, nan, nan, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);
% 
%     hTx = quiver3(axAnim, 0,0,0,0,0,0,0,'r','LineWidth',1.8);
%     hTy = quiver3(axAnim, 0,0,0,0,0,0,0,'g','LineWidth',1.8);
%     hTz = quiver3(axAnim, 0,0,0,0,0,0,0,'b','LineWidth',1.8);
% 
%     hEx = quiver3(axAnim, 0,0,0,0,0,0,0,'r','LineWidth',1.2);
%     hEy = quiver3(axAnim, 0,0,0,0,0,0,0,'g','LineWidth',1.2);
%     hEz = quiver3(axAnim, 0,0,0,0,0,0,0,'b','LineWidth',1.2);
% 
%     hNoseTrue = plot3(axAnim, nan, nan, nan, 'b', 'LineWidth', 2.2);
%     hNoseEst  = plot3(axAnim, nan, nan, nan, '--r', 'LineWidth', 1.8);
% 
%     legend(axAnim, 'truth trail','estimate trail','truth pos','estimate pos','Location','best');
% 
%     t0 = tic;
% 
%     for k = 1:step:N
%         i0 = max(1, k-trail_len);
% 
%         set(hTrailTrue, 'XData', x_true(i0:k), 'YData', y_true(i0:k), 'ZData', z_true(i0:k));
%         set(hTrailEst,  'XData', x_est(i0:k),  'YData', y_est(i0:k),  'ZData', z_est(i0:k));
% 
%         set(hBodyTrue, 'XData', x_true(k), 'YData', y_true(k), 'ZData', z_true(k));
%         set(hBodyEst,  'XData', x_est(k),  'YData', y_est(k),  'ZData', z_est(k));
% 
%         [oT, exT, eyT, ezT] = pose_axes_for_plot_local(p_true(:,k), sim.q_true(:,k), axis_len_truth);
% 
%         set(hTx, 'XData', oT(1), 'YData', oT(2), 'ZData', oT(3), ...
%                  'UData', exT(1), 'VData', exT(2), 'WData', exT(3));
%         set(hTy, 'XData', oT(1), 'YData', oT(2), 'ZData', oT(3), ...
%                  'UData', eyT(1), 'VData', eyT(2), 'WData', eyT(3));
%         set(hTz, 'XData', oT(1), 'YData', oT(2), 'ZData', oT(3), ...
%                  'UData', ezT(1), 'VData', ezT(2), 'WData', ezT(3));
% 
%         [oE, exE, eyE, ezE] = pose_axes_for_plot_local(p_est(:,k), log_q(:,k), axis_len_est);
% 
%         set(hEx, 'XData', oE(1), 'YData', oE(2), 'ZData', oE(3), ...
%                  'UData', exE(1), 'VData', exE(2), 'WData', exE(3));
%         set(hEy, 'XData', oE(1), 'YData', oE(2), 'ZData', oE(3), ...
%                  'UData', eyE(1), 'VData', eyE(2), 'WData', eyE(3));
%         set(hEz, 'XData', oE(1), 'YData', oE(2), 'ZData', oE(3), ...
%                  'UData', ezE(1), 'VData', ezE(2), 'WData', ezE(3));
% 
%         set(hNoseTrue, 'XData', [oT(1), oT(1)+exT(1)], ...
%                        'YData', [oT(2), oT(2)+exT(2)], ...
%                        'ZData', [oT(3), oT(3)+exT(3)]);
% 
%         set(hNoseEst,  'XData', [oE(1), oE(1)+exE(1)], ...
%                        'YData', [oE(2), oE(2)+exE(2)], ...
%                        'ZData', [oE(3), oE(3)+exE(3)]);
% 
%         title(axAnim, sprintf('Truth vs estimate UAV motion   |   t = %.2f s', t(k)));
% 
%         drawnow limitrate;
% 
%         target_elapsed = (t(k) - t(1)) / playback_rate;
%         current_elapsed = toc(t0);
% 
%         if current_elapsed < target_elapsed
%             pause(target_elapsed - current_elapsed);
%         end
%     end
% end

function [origin_plot, ex_plot, ey_plot, ez_plot] = pose_axes_for_plot_local(p_n, q_nb, L)
    R_nb = quat_to_rotmat_local(q_nb);

    origin_plot = [p_n(2); p_n(1); -p_n(3)];

    ex_n = R_nb(:,1);
    ey_n = R_nb(:,2);
    ez_n = R_nb(:,3);

    ex_plot = L * [ex_n(2); ex_n(1); -ex_n(3)];
    ey_plot = L * [ey_n(2); ey_n(1); -ey_n(3)];
    ez_plot = L * [ez_n(2); ez_n(1); -ez_n(3)];
end

function R = quat_to_rotmat_local(q)
    q = q(:);
    q = q / norm(q);

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    R = [1 - 2*(qy^2 + qz^2),   2*(qx*qy - qz*qw),   2*(qx*qz + qy*qw);
         2*(qx*qy + qz*qw),     1 - 2*(qx^2 + qz^2), 2*(qy*qz - qx*qw);
         2*(qx*qz - qy*qw),     2*(qy*qz + qx*qw),   1 - 2*(qx^2 + qy^2)];
end