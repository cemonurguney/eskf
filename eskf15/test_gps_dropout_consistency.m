clear; clc; close all;

%% ============================================================
% Load latest dropout run
% ============================================================

load("fixedwing_run_gps_dropout_observability.mat");

fprintf("\n=== LOADED RUN CHECK ===\n");

if exist("GPS_OUTAGE_POS_INFLATE_NE","var")
    fprintf("GPS_OUTAGE_POS_INFLATE_NE   = [%.3f %.3f]\n", GPS_OUTAGE_POS_INFLATE_NE);
else
    warning("GPS_OUTAGE_POS_INFLATE_NE variable not found. MAT file may be old.");
end

if exist("GPS_OUTAGE_VEL_INFLATE_NE","var")
    fprintf("GPS_OUTAGE_VEL_INFLATE_NE   = [%.3f %.3f]\n", GPS_OUTAGE_VEL_INFLATE_NE);
else
    warning("GPS_OUTAGE_VEL_INFLATE_NE variable not found. MAT file may be old.");
end

if exist("GPS_OUTAGE_SIGMA_A_SCALE","var")
    fprintf("GPS_OUTAGE_SIGMA_A_SCALE     = %.3f\n", GPS_OUTAGE_SIGMA_A_SCALE);
end

if exist("GPS_OUTAGE_SIGMA_G_SCALE","var")
    fprintf("GPS_OUTAGE_SIGMA_G_SCALE     = %.3f\n", GPS_OUTAGE_SIGMA_G_SCALE);
end

if exist("GPS_OUTAGE_SIGMA_BA_RW_SCALE","var")
    fprintf("GPS_OUTAGE_SIGMA_BA_RW_SCALE = %.3f\n", GPS_OUTAGE_SIGMA_BA_RW_SCALE);
end

fprintf("GPS pos updates used : %d\n", used_gps_pos_updates);
fprintf("GPS vel updates used : %d\n", used_gps_vel_updates);
fprintf("BARO updates used    : %d\n", used_baro_updates);
fprintf("AIRS updates used    : %d\n", used_airspeed_updates);

fprintf("GPS scheduled IMU samples       : %d / %d\n", sum(log_use_gps_meas), numel(t));
fprintf("Baro bias estimated IMU samples : %d / %d\n", sum(log_estimate_baro_bias), numel(t));

%% ============================================================
% GPS holdout position consistency
% ============================================================

idx_hold_p = sim.gps_pos_available & all(isfinite(sim.gps_pos),1) & ~log_use_gps_meas;

err_p = log_p(:,idx_hold_p) - sim.gps_pos(:,idx_hold_p);

sig_p = sqrt(log_Pdiag(1:3,idx_hold_p));
sig3_p = 3 * sig_p;

fprintf("\n=== GPS HOLDOUT POSITION CONSISTENCY ===\n");

fprintf("Holdout samples: %d\n", sum(idx_hold_p));

fprintf("Holdout RMSE N/E/D [m] = [%.3f %.3f %.3f]\n", ...
    sqrt(mean(err_p.^2,2,"omitnan")));

fprintf("Mean 1-sigma N/E/D [m] = [%.3f %.3f %.3f]\n", ...
    mean(sig_p,2,"omitnan"));

fprintf("Mean 3-sigma N/E/D [m] = [%.3f %.3f %.3f]\n", ...
    mean(sig3_p,2,"omitnan"));

outside_3sig_p = abs(err_p) > sig3_p;

fprintf("Outside 3-sigma ratio N/E/D [%%] = [%.1f %.1f %.1f]\n", ...
    100*mean(outside_3sig_p,2,"omitnan"));

%% ============================================================
% GPS holdout velocity consistency
% ============================================================

idx_hold_v = sim.gps_vel_available & all(isfinite(sim.gps_vel),1) & ~log_use_gps_meas;

err_v = log_v(:,idx_hold_v) - sim.gps_vel(:,idx_hold_v);

sig_v = sqrt(log_Pdiag(4:6,idx_hold_v));
sig3_v = 3 * sig_v;

fprintf("\n=== GPS HOLDOUT VELOCITY CONSISTENCY ===\n");

fprintf("Holdout velocity samples: %d\n", sum(idx_hold_v));

fprintf("Holdout velocity RMSE N/E/D [m/s] = [%.3f %.3f %.3f]\n", ...
    sqrt(mean(err_v.^2,2,"omitnan")));

fprintf("Mean 1-sigma velocity N/E/D [m/s] = [%.3f %.3f %.3f]\n", ...
    mean(sig_v,2,"omitnan"));

fprintf("Mean 3-sigma velocity N/E/D [m/s] = [%.3f %.3f %.3f]\n", ...
    mean(sig3_v,2,"omitnan"));

outside_3sig_v = abs(err_v) > sig3_v;

fprintf("Outside 3-sigma ratio velocity N/E/D [%%] = [%.1f %.1f %.1f]\n", ...
    100*mean(outside_3sig_v,2,"omitnan"));

%% ============================================================
% Used GPS consistency
% ============================================================

idx_used_p = sim.gps_pos_available & all(isfinite(sim.gps_pos),1) & log_use_gps_meas;
idx_used_v = sim.gps_vel_available & all(isfinite(sim.gps_vel),1) & log_use_gps_meas;

fprintf("\n=== GPS USED MEASUREMENT CONSISTENCY ===\n");

if any(idx_used_p)
    err_used_p = log_p(:,idx_used_p) - sim.gps_pos(:,idx_used_p);

    fprintf("Used GPS position RMSE N/E/D [m] = [%.3f %.3f %.3f]\n", ...
        sqrt(mean(err_used_p.^2,2,"omitnan")));

    fprintf("Used GPS position RMSE norm [m]  = %.3f\n", ...
        sqrt(mean(sum(err_used_p.^2,1),"omitnan")));
end

if any(idx_used_v)
    err_used_v = log_v(:,idx_used_v) - sim.gps_vel(:,idx_used_v);

    fprintf("Used GPS velocity RMSE N/E/D [m/s] = [%.3f %.3f %.3f]\n", ...
        sqrt(mean(err_used_v.^2,2,"omitnan")));

    fprintf("Used GPS velocity RMSE norm [m/s]  = %.3f\n", ...
        sqrt(mean(sum(err_used_v.^2,1),"omitnan")));
end

%% ============================================================
% Baro bias observability / freeze check
% ============================================================

fprintf("\n=== BARO BIAS OBSERVABILITY CHECK ===\n");

fprintf("Initial b_baro [m] = %.6f\n", log_bbaro(1));
fprintf("Final b_baro [m]   = %.6f\n", log_bbaro(end));

idx_gps_off = ~log_use_gps_meas;
idx_gps_on  = log_use_gps_meas;

if any(idx_gps_off)
    fprintf("b_baro std while GPS OFF [m] = %.6f\n", ...
        std(log_bbaro(idx_gps_off),0,"omitnan"));
end

if any(idx_gps_on)
    fprintf("b_baro std while GPS ON  [m] = %.6f\n", ...
        std(log_bbaro(idx_gps_on),0,"omitnan"));
end

fprintf("Baro bias estimated samples: %d\n", sum(log_estimate_baro_bias));

%% ============================================================
% Baro modeled residual
% ============================================================

idx_baro = sim.baro_available & isfinite(sim.baro);

baro_raw_err = log_p(3,idx_baro) - sim.baro(idx_baro);
baro_model_err = log_p(3,idx_baro) + log_bbaro(idx_baro) - sim.baro(idx_baro);

fprintf("\n=== BARO RESIDUAL CHECK ===\n");

fprintf("Baro raw RMSE [m]       = %.3f\n", sqrt(mean(baro_raw_err.^2,"omitnan")));
fprintf("Baro raw mean [m]       = %.3f\n", mean(baro_raw_err,"omitnan"));

fprintf("Baro modeled RMSE [m]   = %.3f\n", sqrt(mean(baro_model_err.^2,"omitnan")));
fprintf("Baro modeled mean [m]   = %.3f\n", mean(baro_model_err,"omitnan"));
fprintf("Baro modeled std [m]    = %.3f\n", std(baro_model_err,0,"omitnan"));

%% ============================================================
% Wind reference check
% ============================================================

if isfield(sim,"wind_ref_available") && any(sim.wind_ref_available)
    idx_wind = sim.wind_ref_available & all(isfinite(sim.wind_ref),1);

    wind_err = log_wind_ne(:,idx_wind) - sim.wind_ref(:,idx_wind);

    fprintf("\n=== WIND VS PX4 REFERENCE CHECK ===\n");

    fprintf("Wind RMSE N/E [m/s]  = [%.3f %.3f]\n", ...
        sqrt(mean(wind_err.^2,2,"omitnan")));

    fprintf("Wind RMSE norm [m/s] = %.3f\n", ...
        sqrt(mean(sum(wind_err.^2,1),"omitnan")));

    fprintf("Wind mean N/E [m/s]  = [%.3f %.3f]\n", ...
        mean(wind_err,2,"omitnan"));
end

%% ============================================================
% Airspeed residual check
% ============================================================

idx_air = isfinite(log_res_airspeed);

fprintf("\n=== AIRSPEED RESIDUAL CHECK ===\n");

fprintf("Airspeed residual RMSE [m/s] = %.3f\n", ...
    sqrt(mean(log_res_airspeed(idx_air).^2,"omitnan")));

fprintf("Airspeed residual mean [m/s] = %.3f\n", ...
    mean(log_res_airspeed(idx_air),"omitnan"));

fprintf("Airspeed residual std [m/s]  = %.3f\n", ...
    std(log_res_airspeed(idx_air),0,"omitnan"));

%% ============================================================
% Plots
% ============================================================

labels = ["N","E","D"];

if any(idx_hold_p)
    th = t(idx_hold_p);

    figure("Name","GPS Holdout Position Error vs 3-sigma");

    for ax = 1:3
        subplot(3,1,ax);
        plot(th, err_p(ax,:), "k.");
        hold on;
        plot(th, sig3_p(ax,:), "r", "LineWidth", 1.1);
        plot(th, -sig3_p(ax,:), "r", "LineWidth", 1.1);

        grid on;
        ylabel(labels(ax) + " error [m]");
        legend("holdout error", "+/-3 sigma");
    end

    xlabel("Time [s]");
    sgtitle("GPS Holdout Position Error vs ESKF 3-sigma");
end

if any(idx_hold_v)
    thv = t(idx_hold_v);

    figure("Name","GPS Holdout Velocity Error vs 3-sigma");

    for ax = 1:3
        subplot(3,1,ax);
        plot(thv, err_v(ax,:), "k.");
        hold on;
        plot(thv, sig3_v(ax,:), "r", "LineWidth", 1.1);
        plot(thv, -sig3_v(ax,:), "r", "LineWidth", 1.1);

        grid on;
        ylabel(labels(ax) + " error [m/s]");
        legend("holdout error", "+/-3 sigma");
    end

    xlabel("Time [s]");
    sgtitle("GPS Holdout Velocity Error vs ESKF 3-sigma");
end

figure("Name","GPS Schedule and Baro Bias");
yyaxis left;
plot(t, log_bbaro, "LineWidth", 1.3);
ylabel("b_{baro} [m]");
grid on;

yyaxis right;
stairs(t, double(log_use_gps_meas), "LineWidth", 1.1);
hold on;
stairs(t, double(log_estimate_baro_bias), "--", "LineWidth", 1.1);
ylabel("Flag");
ylim([-0.1 1.1]);

xlabel("Time [s]");
legend("b_{baro}", "GPS pos+vel scheduled", "baro bias estimated", "Location","best");
title("GPS Schedule and Barometer Bias Observability");

figure("Name","Position Sigma");
plot(t, sqrt(log_Pdiag(1,:)), "LineWidth", 1.2); hold on;
plot(t, sqrt(log_Pdiag(2,:)), "LineWidth", 1.2);
plot(t, sqrt(log_Pdiag(3,:)), "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Position 1-sigma [m]");
legend("sigma_N","sigma_E","sigma_D");
title("Position Covariance Growth");

figure("Name","Velocity Sigma");
plot(t, sqrt(log_Pdiag(4,:)), "LineWidth", 1.2); hold on;
plot(t, sqrt(log_Pdiag(5,:)), "LineWidth", 1.2);
plot(t, sqrt(log_Pdiag(6,:)), "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("Velocity 1-sigma [m/s]");
legend("sigma_VN","sigma_VE","sigma_VD");
title("Velocity Covariance Growth");