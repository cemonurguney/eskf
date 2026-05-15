function px4_cmp = compare_with_px4_reference_from_csv_autoalign( ...
    combined_csv_file, t, log_p, log_v, sim)
%COMPARE_WITH_PX4_REFERENCE_FROM_CSV_AUTOALIGN
% Combined CSV içindeki PX4 local position/velocity reference ile
% bizim ESKF çıktısını karşılaştırır.
%
% Zaman hizasını GPS izi üzerinden otomatik bulur:
%   CSV time = ESKF reset time + tau
%
% Position alignment:
%   initial-window-median-offset kullanılır.
%   Bu sadece local-frame origin farkını düzeltir.
%   Tüm trajectory'ye least-squares fit yapılmaz.
%
% Not:
%   PX4 local output mutlak truth değildir.
%   Estimator-to-estimator consistency reference olarak kullanılır.

    px4_cmp = struct();
    px4_cmp.available = false;
    px4_cmp.file = combined_csv_file;
    px4_cmp.t = t;

    if ~isfile(combined_csv_file)
        warning("Combined CSV file not found: %s", combined_csv_file);
        return;
    end

    %% ============================================================
    % 0) User config
    % ============================================================

    position_alignment_mode = "initial-window-median-offset";
    alignment_window_s = 10.0;

    % GPS auto alignment search ranges
    tau_grid_coarse = -30:0.5:30;
    tau_fine_half_width_s = 2.0;
    tau_fine_step_s = 0.05;

    %% ============================================================
    % 1) Read CSV and check columns
    % ============================================================

    Tref = readtable(combined_csv_file, "VariableNamingRule", "preserve");

    var_names = string(Tref.Properties.VariableNames);

    required_cols = [
        "t_s"
        "gps_available"
        "gps_n_m"
        "gps_e_m"
        "gps_d_m"
        "px4_local_n_m"
        "px4_local_e_m"
        "px4_local_d_m"
        "px4_local_vn_m_s"
        "px4_local_ve_m_s"
        "px4_local_vd_m_s"
    ];

    if ~all(ismember(required_cols, var_names))
        warning("Combined CSV does not contain required reference columns.");
        fprintf("Required columns:\n");
        disp(required_cols);
        fprintf("Available columns:\n");
        disp(var_names.');
        return;
    end

    %% ============================================================
    % 2) CSV arrays
    % ============================================================

    t_csv = Tref.t_s(:).';
    t_csv = t_csv - t_csv(1);

    gps_csv = [
        Tref.gps_n_m(:).'
        Tref.gps_e_m(:).'
        Tref.gps_d_m(:).'
    ];

    gps_csv_valid = ...
        Tref.gps_available(:).' == 1 & ...
        all(isfinite(gps_csv), 1) & ...
        isfinite(t_csv);

    px4_p_ref = [
        Tref.px4_local_n_m(:).'
        Tref.px4_local_e_m(:).'
        Tref.px4_local_d_m(:).'
    ];

    px4_v_ref = [
        Tref.px4_local_vn_m_s(:).'
        Tref.px4_local_ve_m_s(:).'
        Tref.px4_local_vd_m_s(:).'
    ];

    valid_px4_p = all(isfinite(px4_p_ref), 1) & isfinite(t_csv);
    valid_px4_v = all(isfinite(px4_v_ref), 1) & isfinite(t_csv);

    %% ============================================================
    % 3) Sim GPS target for time alignment
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
    % 4) Coarse tau search using GPS NE track
    % ============================================================

    best_tau = 0;
    best_score = inf;
    best_offset_gps = [0;0;0];
    best_valid_count = 0;

    for ii = 1:numel(tau_grid_coarse)
        tau = tau_grid_coarse(ii);

        [score_i, offset_i, valid_count_i] = score_gps_alignment_local( ...
            tau, t_sim_gps, gps_sim, t_csv_gps, gps_csv_data);

        if isfinite(score_i) && score_i < best_score
            best_score = score_i;
            best_tau = tau;
            best_offset_gps = offset_i;
            best_valid_count = valid_count_i;
        end
    end

    %% ============================================================
    % 5) Fine tau search
    % ============================================================

    tau_fine = ...
        (best_tau - tau_fine_half_width_s) : ...
        tau_fine_step_s : ...
        (best_tau + tau_fine_half_width_s);

    for ii = 1:numel(tau_fine)
        tau = tau_fine(ii);

        [score_i, offset_i, valid_count_i] = score_gps_alignment_local( ...
            tau, t_sim_gps, gps_sim, t_csv_gps, gps_csv_data);

        if isfinite(score_i) && score_i < best_score
            best_score = score_i;
            best_tau = tau;
            best_offset_gps = offset_i;
            best_valid_count = valid_count_i;
        end
    end

    %% ============================================================
    % 6) If search failed, force tau = 0
    % ============================================================

    if ~isfinite(best_score)
        warning("Auto alignment failed. Forcing tau = 0 over common overlap.");

        best_tau = 0;

        [best_score, best_offset_gps, best_valid_count] = score_gps_alignment_local( ...
            best_tau, t_sim_gps, gps_sim, t_csv_gps, gps_csv_data);
    end

    %% ============================================================
    % 7) Build best aligned GPS segment for plotting
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
    % 8) Interpolate PX4 reference at aligned time
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
    % 9) GPS schedule for split metrics
    % ============================================================

    gps_sched = true(1, numel(t));

    if isfield(sim, "gps_meas_scheduled")
        gps_sched = logical(sim.gps_meas_scheduled(:).');
    end

    if numel(gps_sched) ~= numel(t)
        gps_sched = true(1, numel(t));
    end

    %% ============================================================
    % 10) Position origin alignment, robust initial-window median
    % ============================================================

    valid_p_i = all(isfinite(px4_p_i), 1) & all(isfinite(log_p), 1);

    if any(valid_p_i)

        first_valid = find(valid_p_i, 1, "first");
        t0_align = t(first_valid);

        init_align = valid_p_i & ...
            (t >= t0_align) & ...
            (t <= t0_align + alignment_window_s);

        if sum(init_align) < 5
            warning("[PX4 compare] Initial alignment window has too few samples. Falling back to first valid sample.");

            offset_p = log_p(:,first_valid) - px4_p_i(:,first_valid);
            used_alignment_mode = "first-valid-sample-offset";
            alignment_sample_count = 1;

        else
            switch position_alignment_mode

                case "initial-window-median-offset"
                    offset_p = median( ...
                        log_p(:,init_align) - px4_p_i(:,init_align), ...
                        2, "omitnan");

                case "initial-window-mean-offset"
                    offset_p = mean( ...
                        log_p(:,init_align) - px4_p_i(:,init_align), ...
                        2, "omitnan");

                case "first-sample-offset"
                    offset_p = log_p(:,first_valid) - px4_p_i(:,first_valid);

                case "gps-median-offset"
                    offset_p = best_offset_gps;

                otherwise
                    error("Unknown position_alignment_mode: %s", position_alignment_mode);
            end

            used_alignment_mode = position_alignment_mode;
            alignment_sample_count = sum(init_align);
        end

        px4_p_i_aligned = px4_p_i + offset_p;

        err_p_px4_raw = log_p - px4_p_i;
        err_p_px4 = log_p - px4_p_i_aligned;

        valid_err_p = all(isfinite(err_p_px4), 1);

        valid_err_p_on  = valid_err_p & gps_sched;
        valid_err_p_off = valid_err_p & ~gps_sched;

        %% All-run position metrics
        pos_rmse_axis_raw = rmse_axis_local(err_p_px4_raw, valid_err_p);
        pos_rmse_norm_raw = rmse_norm_local(err_p_px4_raw, valid_err_p);

        pos_rmse_axis = rmse_axis_local(err_p_px4, valid_err_p);
        pos_rmse_norm = rmse_norm_local(err_p_px4, valid_err_p);
        pos_mean_axis = mean_axis_local(err_p_px4, valid_err_p);
        pos_std_axis  = std_axis_local(err_p_px4, valid_err_p);

        %% GPS ON/OFF position metrics
        pos_rmse_axis_on = rmse_axis_local(err_p_px4, valid_err_p_on);
        pos_rmse_norm_on = rmse_norm_local(err_p_px4, valid_err_p_on);
        pos_mean_axis_on = mean_axis_local(err_p_px4, valid_err_p_on);

        pos_rmse_axis_off = rmse_axis_local(err_p_px4, valid_err_p_off);
        pos_rmse_norm_off = rmse_norm_local(err_p_px4, valid_err_p_off);
        pos_mean_axis_off = mean_axis_local(err_p_px4, valid_err_p_off);

        %% Store
        px4_cmp.available = true;

        px4_cmp.position_alignment_mode = used_alignment_mode;
        px4_cmp.position_alignment_window_s = alignment_window_s;
        px4_cmp.position_alignment_sample_count = alignment_sample_count;
        px4_cmp.position_alignment_offset = offset_p;

        px4_cmp.p_ref_raw = px4_p_i;
        px4_cmp.p_ref = px4_p_i_aligned;

        px4_cmp.err_p_raw = err_p_px4_raw;
        px4_cmp.err_p = err_p_px4;

        px4_cmp.valid_err_p = valid_err_p;
        px4_cmp.valid_err_p_on = valid_err_p_on;
        px4_cmp.valid_err_p_off = valid_err_p_off;

        px4_cmp.pos_rmse_axis_raw = pos_rmse_axis_raw;
        px4_cmp.pos_rmse_norm_raw = pos_rmse_norm_raw;

        px4_cmp.pos_rmse_axis = pos_rmse_axis;
        px4_cmp.pos_rmse_norm = pos_rmse_norm;
        px4_cmp.pos_mean_axis = pos_mean_axis;
        px4_cmp.pos_std_axis = pos_std_axis;

        px4_cmp.pos_rmse_axis_on = pos_rmse_axis_on;
        px4_cmp.pos_rmse_norm_on = pos_rmse_norm_on;
        px4_cmp.pos_mean_axis_on = pos_mean_axis_on;

        px4_cmp.pos_rmse_axis_off = pos_rmse_axis_off;
        px4_cmp.pos_rmse_norm_off = pos_rmse_norm_off;
        px4_cmp.pos_mean_axis_off = pos_mean_axis_off;

        %% Print
        fprintf("\n=== ESKF vs PX4 LOCAL POSITION REFERENCE ===\n");
        fprintf("Position alignment mode          : %s\n", used_alignment_mode);
        fprintf("Alignment window [s]             : %.3f\n", alignment_window_s);
        fprintf("Alignment samples                : %d\n", alignment_sample_count);

        fprintf("Raw Position RMSE N/E/D [m]      : [%.3f %.3f %.3f]\n", pos_rmse_axis_raw);
        fprintf("Raw Position RMSE norm [m]       : %.3f\n", pos_rmse_norm_raw);

        fprintf("Applied origin offset N/E/D [m]  : [%.3f %.3f %.3f]\n", offset_p);

        fprintf("Aligned Position RMSE N/E/D [m]  : [%.3f %.3f %.3f]\n", pos_rmse_axis);
        fprintf("Aligned Position RMSE norm [m]   : %.3f\n", pos_rmse_norm);
        fprintf("Aligned mean error N/E/D [m]     : [%.3f %.3f %.3f]\n", pos_mean_axis);

        if any(valid_err_p_on)
            fprintf("Aligned Position RMSE GPS ON N/E/D [m] : [%.3f %.3f %.3f]\n", pos_rmse_axis_on);
            fprintf("Aligned Position RMSE GPS ON norm [m]  : %.3f\n", pos_rmse_norm_on);
        end

        if any(valid_err_p_off)
            fprintf("Aligned Position RMSE GPS OFF N/E/D [m]: [%.3f %.3f %.3f]\n", pos_rmse_axis_off);
            fprintf("Aligned Position RMSE GPS OFF norm [m] : %.3f\n", pos_rmse_norm_off);
        end

        %% Plots
        figure('Name','PX4 Auto Alignment Check');
        plot(gps_csv(2,gps_csv_valid), gps_csv(1,gps_csv_valid), ...
            "Color", [0.7 0.7 0.7]);
        hold on;

        if any(valid_best_gps)
            plot( ...
                gps_csv_best(2,valid_best_gps) + best_offset_gps(2), ...
                gps_csv_best(1,valid_best_gps) + best_offset_gps(1), ...
                "b.");
        end

        plot(gps_sim(2,:), gps_sim(1,:), "r.");

        grid on;
        axis equal;
        xlabel("East [m]");
        ylabel("North [m]");
        legend("Full CSV GPS", "Aligned CSV GPS segment", "sim GPS segment", ...
            "Location", "best");
        title(sprintf("GPS Auto Alignment, tau = %.2f s", best_tau));

        figure('Name','ESKF Position Error w.r.t. PX4 EKF Reference');
        hold on;

        shade_gps_off_regions_local(t, gps_sched);

        plot(t(valid_err_p), err_p_px4(1,valid_err_p), "r");
        plot(t(valid_err_p), err_p_px4(2,valid_err_p), "g");
        plot(t(valid_err_p), err_p_px4(3,valid_err_p), "b");

        grid on;
        xlabel("Time [s]");
        ylabel("Position error wrt aligned PX4 ref [m]");
        legend("N", "E", "D", "Location", "best");
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
        legend("PX4 raw ref", "PX4 origin-aligned ref", "Our ESKF", ...
            "Location", "best");
        title("Ground Track: ESKF vs PX4 Reference");

        figure('Name','Altitude: ESKF vs PX4 Reference');
        hold on;

        shade_gps_off_regions_local(t, gps_sched);

        plot(t, -px4_p_i_aligned(3,:), "k--", "LineWidth", 1.2);
        plot(t, -log_p(3,:), "r", "LineWidth", 1.2);

        grid on;
        xlabel("Time [s]");
        ylabel("Altitude Up [m]");
        legend("PX4 EKF origin-aligned reference", "Our ESKF", ...
            "Location", "best");
        title("Altitude: ESKF vs PX4 Reference");

    else
        warning("PX4 local position reference could not be interpolated.");
    end

    %% ============================================================
    % 11) Velocity comparison
    % ============================================================

    valid_v_i = all(isfinite(px4_v_i), 1) & all(isfinite(log_v), 1);

    if any(valid_v_i)

        err_v_px4 = log_v - px4_v_i;

        valid_err_v = all(isfinite(err_v_px4), 1);

        valid_err_v_on  = valid_err_v & gps_sched;
        valid_err_v_off = valid_err_v & ~gps_sched;

        vel_rmse_axis = rmse_axis_local(err_v_px4, valid_err_v);
        vel_rmse_norm = rmse_norm_local(err_v_px4, valid_err_v);
        vel_mean_axis = mean_axis_local(err_v_px4, valid_err_v);
        vel_std_axis  = std_axis_local(err_v_px4, valid_err_v);

        vel_rmse_axis_on = rmse_axis_local(err_v_px4, valid_err_v_on);
        vel_rmse_norm_on = rmse_norm_local(err_v_px4, valid_err_v_on);
        vel_mean_axis_on = mean_axis_local(err_v_px4, valid_err_v_on);

        vel_rmse_axis_off = rmse_axis_local(err_v_px4, valid_err_v_off);
        vel_rmse_norm_off = rmse_norm_local(err_v_px4, valid_err_v_off);
        vel_mean_axis_off = mean_axis_local(err_v_px4, valid_err_v_off);

        px4_cmp.available = true;

        px4_cmp.v_ref = px4_v_i;
        px4_cmp.err_v = err_v_px4;

        px4_cmp.valid_err_v = valid_err_v;
        px4_cmp.valid_err_v_on = valid_err_v_on;
        px4_cmp.valid_err_v_off = valid_err_v_off;

        px4_cmp.vel_rmse_axis = vel_rmse_axis;
        px4_cmp.vel_rmse_norm = vel_rmse_norm;
        px4_cmp.vel_mean_axis = vel_mean_axis;
        px4_cmp.vel_std_axis = vel_std_axis;

        px4_cmp.vel_rmse_axis_on = vel_rmse_axis_on;
        px4_cmp.vel_rmse_norm_on = vel_rmse_norm_on;
        px4_cmp.vel_mean_axis_on = vel_mean_axis_on;

        px4_cmp.vel_rmse_axis_off = vel_rmse_axis_off;
        px4_cmp.vel_rmse_norm_off = vel_rmse_norm_off;
        px4_cmp.vel_mean_axis_off = vel_mean_axis_off;

        fprintf("\n=== ESKF vs PX4 LOCAL VELOCITY REFERENCE ===\n");
        fprintf("Velocity RMSE N/E/D [m/s]        : [%.3f %.3f %.3f]\n", vel_rmse_axis);
        fprintf("Velocity RMSE norm [m/s]         : %.3f\n", vel_rmse_norm);
        fprintf("Velocity mean error [m/s]        : [%.3f %.3f %.3f]\n", vel_mean_axis);

        if any(valid_err_v_on)
            fprintf("Velocity RMSE GPS ON N/E/D [m/s] : [%.3f %.3f %.3f]\n", vel_rmse_axis_on);
            fprintf("Velocity RMSE GPS ON norm [m/s]  : %.3f\n", vel_rmse_norm_on);
        end

        if any(valid_err_v_off)
            fprintf("Velocity RMSE GPS OFF N/E/D [m/s]: [%.3f %.3f %.3f]\n", vel_rmse_axis_off);
            fprintf("Velocity RMSE GPS OFF norm [m/s] : %.3f\n", vel_rmse_norm_off);
        end

        figure('Name','ESKF Velocity Error w.r.t. PX4 EKF Reference');
        hold on;

        shade_gps_off_regions_local(t, gps_sched);

        plot(t(valid_v_i), err_v_px4(1,valid_v_i), "r");
        plot(t(valid_v_i), err_v_px4(2,valid_v_i), "g");
        plot(t(valid_v_i), err_v_px4(3,valid_v_i), "b");

        grid on;
        xlabel("Time [s]");
        ylabel("Velocity error wrt PX4 ref [m/s]");
        legend("V_N", "V_E", "V_D", "Location", "best");
        title("ESKF Velocity Error w.r.t. PX4 EKF Reference");

        figure('Name','Velocity: ESKF vs PX4 Reference');

        subplot(3,1,1);
        hold on;
        shade_gps_off_regions_local(t, gps_sched);
        plot(t(valid_v_i), px4_v_i(1,valid_v_i), "k--");
        plot(t(valid_v_i), log_v(1,valid_v_i), "r");
        grid on;
        ylabel("V_N [m/s]");
        legend("PX4 ref", "ESKF", "Location", "best");

        subplot(3,1,2);
        hold on;
        shade_gps_off_regions_local(t, gps_sched);
        plot(t(valid_v_i), px4_v_i(2,valid_v_i), "k--");
        plot(t(valid_v_i), log_v(2,valid_v_i), "r");
        grid on;
        ylabel("V_E [m/s]");
        legend("PX4 ref", "ESKF", "Location", "best");

        subplot(3,1,3);
        hold on;
        shade_gps_off_regions_local(t, gps_sched);
        plot(t(valid_v_i), px4_v_i(3,valid_v_i), "k--");
        plot(t(valid_v_i), log_v(3,valid_v_i), "r");
        grid on;
        xlabel("Time [s]");
        ylabel("V_D [m/s]");
        legend("PX4 ref", "ESKF", "Location", "best");

        sgtitle("Velocity: ESKF vs PX4 Reference");

    else
        warning("PX4 local velocity reference could not be interpolated.");
    end
end

%% ============================================================
% Local helper: GPS alignment score
% ============================================================
function [score, offset, valid_count] = score_gps_alignment_local( ...
    tau, t_sim_gps, gps_sim, t_csv_gps, gps_csv_data)

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

    if sum(valid_i) < max(10, 0.60*numel(t_sim_gps))
        score = inf;
        offset = [0;0;0];
        valid_count = sum(valid_i);
        return;
    end

    % offset = sim GPS - CSV GPS
    offset = median(gps_sim(:,valid_i) - gps_csv_i(:,valid_i), 2, "omitnan");

    err_i = gps_sim(:,valid_i) - (gps_csv_i(:,valid_i) + offset);

    % Time alignment score only on horizontal NE.
    score = sqrt(mean(sum(err_i(1:2,:).^2, 1), "omitnan"));

    valid_count = sum(valid_i);
end

%% ============================================================
% Local helper: RMSE / mean / std
% ============================================================
function r = rmse_axis_local(err, valid)
    if nargin < 2 || isempty(valid)
        valid = all(isfinite(err), 1);
    end

    if ~any(valid)
        r = nan(size(err,1),1);
        return;
    end

    r = sqrt(mean(err(:,valid).^2, 2, "omitnan"));
end

function r = rmse_norm_local(err, valid)
    if nargin < 2 || isempty(valid)
        valid = all(isfinite(err), 1);
    end

    if ~any(valid)
        r = nan;
        return;
    end

    r = sqrt(mean(sum(err(:,valid).^2, 1), "omitnan"));
end

function m = mean_axis_local(err, valid)
    if nargin < 2 || isempty(valid)
        valid = all(isfinite(err), 1);
    end

    if ~any(valid)
        m = nan(size(err,1),1);
        return;
    end

    m = mean(err(:,valid), 2, "omitnan");
end

function s = std_axis_local(err, valid)
    if nargin < 2 || isempty(valid)
        valid = all(isfinite(err), 1);
    end

    if ~any(valid)
        s = nan(size(err,1),1);
        return;
    end

    s = std(err(:,valid), 0, 2, "omitnan");
end

%% ============================================================
% Local helper: shade GPS OFF regions
% ============================================================
function shade_gps_off_regions_local(t, gps_flag)
%SHADE_GPS_OFF_REGIONS_LOCAL
% GPS OFF bölgelerini mevcut eksende gri arka plan olarak gösterir.
% Legend'a dahil edilmez. MATLAB'ın küçük legend krizi böyle önlenir.

    gps_flag = logical(gps_flag(:).');
    t = t(:).';

    if numel(gps_flag) ~= numel(t)
        return;
    end

    off_flag = ~gps_flag;

    if ~any(off_flag)
        return;
    end

    edges = diff([false, off_flag, false]);
    starts = find(edges == 1);
    stops  = find(edges == -1) - 1;

    yl = ylim;

    for i = 1:numel(starts)
        xs = t(starts(i));
        xe = t(stops(i));

        patch([xs xe xe xs], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.90 0.90 0.90], ...
            "EdgeColor", "none", ...
            "FaceAlpha", 0.45, ...
            "HandleVisibility", "off");
    end

    h_lines = findobj(gca, "Type", "line");

    if ~isempty(h_lines)
        uistack(h_lines, "top");
    end
end