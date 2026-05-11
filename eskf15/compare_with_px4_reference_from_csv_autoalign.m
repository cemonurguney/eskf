function px4_cmp = compare_with_px4_reference_from_csv_autoalign( ...
    combined_csv_file, t, log_p, log_v, sim)
%COMPARE_WITH_PX4_REFERENCE_FROM_CSV_AUTOALIGN
% Combined CSV içindeki PX4 local position/velocity reference ile
% bizim ESKF çıktısını karşılaştırır.
%
% Zaman hizasını GPS izi üzerinden otomatik bulur:
%   CSV time = ESKF reset time + tau
%
% Not:
%   PX4 local output mutlak truth değildir.
%   Sadece onboard PX4 EKF reference olarak kullanılır.

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

    %% ============================================================
    % CSV arrays
    % ============================================================

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

    %% ============================================================
    % Sim GPS target
    % ============================================================

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

    %% ============================================================
    % Search tau robustly using overlap
    % ============================================================

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
            gps_csv_i(ax,:) = interp1( ...
                t_csv_gps, ...
                gps_csv_data(ax,:), ...
                tq, ...
                "linear", ...
                nan);
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

    %% ============================================================
    % Fine search
    % ============================================================

    tau_fine = best_tau-2:0.05:best_tau+2;

    for ii = 1:numel(tau_fine)
        tau = tau_fine(ii);

        tq = t_sim_gps + tau;

        gps_csv_i = nan(3, numel(tq));

        for ax = 1:3
            gps_csv_i(ax,:) = interp1( ...
                t_csv_gps, ...
                gps_csv_data(ax,:), ...
                tq, ...
                "linear", ...
                nan);
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

    %% ============================================================
    % If search failed, force tau = 0
    % ============================================================

    if ~isfinite(best_score)
        warning("Auto alignment failed. Forcing tau = 0 over common overlap.");

        best_tau = 0;
        tq = t_sim_gps;

        gps_csv_i = nan(3, numel(tq));

        for ax = 1:3
            gps_csv_i(ax,:) = interp1( ...
                t_csv_gps, ...
                gps_csv_data(ax,:), ...
                tq, ...
                "linear", ...
                nan);
        end

        valid_i = all(isfinite(gps_csv_i), 1) & all(isfinite(gps_sim), 1);

        if sum(valid_i) >= 10
            best_offset_gps = median(gps_sim(:,valid_i) - gps_csv_i(:,valid_i), 2, "omitnan");
            err_i = gps_sim(:,valid_i) - (gps_csv_i(:,valid_i) + best_offset_gps);
            best_score = sqrt(mean(sum(err_i(1:2,:).^2, 1), "omitnan"));
            best_valid_count = sum(valid_i);
        end
    end

    %% ============================================================
    % Build best aligned GPS segment for plotting
    % ============================================================

    tq_best = t_sim_gps + best_tau;
    gps_csv_best = nan(3, numel(tq_best));

    for ax = 1:3
        gps_csv_best(ax,:) = interp1( ...
            t_csv_gps, ...
            gps_csv_data(ax,:), ...
            tq_best, ...
            "linear", ...
            nan);
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

    %% ============================================================
    % Interpolate PX4 reference at aligned time
    % ============================================================

    tq_all = t + best_tau;

    px4_p_i = nan(3, numel(t));
    px4_v_i = nan(3, numel(t));

    if any(valid_px4_p)
        for ax = 1:3
            px4_p_i(ax,:) = interp1( ...
                t_csv(valid_px4_p), ...
                px4_p_ref(ax,valid_px4_p), ...
                tq_all, ...
                "linear", ...
                nan);
        end
    end

    if any(valid_px4_v)
        for ax = 1:3
            px4_v_i(ax,:) = interp1( ...
                t_csv(valid_px4_v), ...
                px4_v_ref(ax,valid_px4_v), ...
                tq_all, ...
                "linear", ...
                nan);
        end
    end

    %% ============================================================
    % Position origin alignment
    % ============================================================

    valid_p_i = all(isfinite(px4_p_i), 1) & all(isfinite(log_p), 1);

    if any(valid_p_i)
        first_valid = find(valid_p_i, 1, "first");

        offset_p_first = log_p(:,first_valid) - px4_p_i(:,first_valid);

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
        plot(px4_p_i(2,:), px4_p_i(1,:), ...
            "Color", [0.5 0.5 0.5], ...
            "LineStyle", "--", ...
            "LineWidth", 1.0);
        hold on;
        plot(px4_p_i_aligned(2,:), px4_p_i_aligned(1,:), ...
            "k--", "LineWidth", 1.2);
        plot(log_p(2,:), log_p(1,:), ...
            "r", "LineWidth", 1.2);
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

    %% ============================================================
    % Velocity comparison
    % ============================================================

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