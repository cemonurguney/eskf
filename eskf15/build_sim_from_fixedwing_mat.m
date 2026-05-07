function sim = build_sim_from_fixedwing_mat(mat_file, t_start, t_end)
%BUILD_SIM_FROM_FIXEDWING_MAT
% Fixed-wing PX4 ULog'dan üretilmiş *_fixedwing_eskf.mat dosyasını,
% mevcut 15-state ESKF main kodunun beklediği sim struct formatına çevirir.
%
% Frame:
%   NED
%   p = [N; E; D]
%   v = [Vn; Ve; Vd]
%
% Baro:
%   baro_z_down_m doğrudan p_D ölçümü gibi kullanılır.
%   Bu versiyonda manuel baro offset uygulanmaz.

    if nargin < 2
        t_start = -inf;
    end

    if nargin < 3
        t_end = inf;
    end

    fw = load(mat_file);

    %% ---------------- IMU trim ----------------
    imu_t = fw.imu_t(:);

    imu_keep = imu_t >= t_start & imu_t <= t_end ...
        & all(isfinite(fw.gyro_rad_s),2) ...
        & all(isfinite(fw.acc_m_s2),2);

    imu_t = imu_t(imu_keep);
    gyro_rad_s = fw.gyro_rad_s(imu_keep,:);
    acc_m_s2   = fw.acc_m_s2(imu_keep,:);

    if isempty(imu_t)
        error("Trim sonrası IMU verisi boş kaldı. t_start/t_end aralığını kontrol et.");
    end

    %% ---------------- Time reset ----------------
    t0 = imu_t(1);

    imu_t = imu_t - t0;
    t = imu_t(:).';

    N = numel(t);

    if N < 2
        error("IMU sample sayısı yetersiz.");
    end

    dt_vec = [median(diff(t)), diff(t)];
    dt_nom = median(diff(t));

    %% ---------------- Base sim ----------------
    sim.t = t;
    sim.dt = dt_nom;
    sim.dt_vec = dt_vec;

    sim.imu_gyro  = gyro_rad_s.';  % 3xN
    sim.imu_accel = acc_m_s2.';    % 3xN

    %% ---------------- GPS position / velocity onto IMU grid ----------------
    gps_t = fw.gps_t(:);

    gps_keep = gps_t >= t_start & gps_t <= t_end ...
        & all(isfinite(fw.gps_pos_ned_m),2);

    gps_t = gps_t(gps_keep) - t0;
    gps_pos = fw.gps_pos_ned_m(gps_keep,:);

    gps_vel = fw.gps_vel_ned_m_s(gps_keep,:);
    gps_vel_valid = all(isfinite(gps_vel),2);

    sim.gps_pos = nan(3,N);
    sim.gps_vel = nan(3,N);

    sim.gps_pos_available = false(1,N);
    sim.gps_vel_available = false(1,N);

    for i = 1:numel(gps_t)
        [~, k] = min(abs(t - gps_t(i)));

        if k >= 1 && k <= N
            sim.gps_pos(:,k) = gps_pos(i,:).';
            sim.gps_pos_available(k) = true;

            if gps_vel_valid(i)
                sim.gps_vel(:,k) = gps_vel(i,:).';
                sim.gps_vel_available(k) = true;
            end
        end
    end

    %% ---------------- Baro onto IMU grid ----------------
    if isfield(fw, "baro_t") && isfield(fw, "baro_z_down_m")
        baro_t = fw.baro_t(:);
        baro_z_down_m = fw.baro_z_down_m(:);

        baro_keep = baro_t >= t_start & baro_t <= t_end ...
            & isfinite(baro_z_down_m);

        baro_t = baro_t(baro_keep) - t0;
        baro_z_down_m = baro_z_down_m(baro_keep);

        sim.baro = nan(1,N);
        sim.baro_available = false(1,N);

        for i = 1:numel(baro_t)
            [~, k] = min(abs(t - baro_t(i)));

            if k >= 1 && k <= N
                sim.baro(k) = baro_z_down_m(i);
                sim.baro_available(k) = true;
            end
        end
    else
        sim.baro = nan(1,N);
        sim.baro_available = false(1,N);
    end

    %% ---------------- Attitude reference/debug ----------------
    if isfield(fw, "att_t") && ~isempty(fw.att_t) && isfield(fw, "q_nb")
        att_t = fw.att_t(:);

        att_keep = att_t >= t_start & att_t <= t_end ...
            & all(isfinite(fw.q_nb),2);

        sim.att_t = att_t(att_keep).' - t0;
        sim.q_ref = fw.q_nb(att_keep,:).';      % 4xM

        if isfield(fw, "rpy_rad")
            sim.rpy_ref = fw.rpy_rad(att_keep,:).'; % 3xM
        else
            sim.rpy_ref = [];
        end
    else
        sim.att_t = [];
        sim.q_ref = [];
        sim.rpy_ref = [];
    end

    %% ---------------- Extra observations for plotting ----------------
    % Airspeed
    if isfield(fw, "airspeed_t") && ~isempty(fw.airspeed_t)
        airspeed_t = fw.airspeed_t(:);
        air_keep = airspeed_t >= t_start & airspeed_t <= t_end;

        sim.airspeed_t = airspeed_t(air_keep).' - t0;

        if isfield(fw, "ias_m_s") && ~isempty(fw.ias_m_s)
            ias = fw.ias_m_s(:);
            sim.ias_m_s = ias(air_keep).';
        else
            sim.ias_m_s = [];
        end

        if isfield(fw, "tas_m_s") && ~isempty(fw.tas_m_s)
            tas = fw.tas_m_s(:);
            sim.tas_m_s = tas(air_keep).';
        else
            sim.tas_m_s = [];
        end
    else
        sim.airspeed_t = [];
        sim.ias_m_s = [];
        sim.tas_m_s = [];
    end

    % Wind estimate from PX4, if available
    if isfield(fw, "wind_t") && ~isempty(fw.wind_t) && isfield(fw, "wind_ned_m_s")
        wind_t = fw.wind_t(:);

        wind_keep = wind_t >= t_start & wind_t <= t_end ...
            & all(isfinite(fw.wind_ned_m_s),2);

        sim.wind_t = wind_t(wind_keep).' - t0;
        sim.wind_ned_m_s = fw.wind_ned_m_s(wind_keep,:).';
    else
        sim.wind_t = [];
        sim.wind_ned_m_s = [];
    end

    sim.imu_acc_norm = vecnorm(sim.imu_accel, 2, 1);
    sim.imu_gyro_norm = vecnorm(sim.imu_gyro, 2, 1);

    %% ---------------- Reference placeholders ----------------
    % Gerçek truth yok. Bunlar sadece eski plotların patlamaması için.
    sim.p_true = nan(3,N);
    sim.v_true = nan(3,N);
    sim.q_true = nan(4,N);

    if any(sim.gps_pos_available)
        idx_gps = find(sim.gps_pos_available);

        for ax = 1:3
            sim.p_true(ax,:) = interp1(t(idx_gps), sim.gps_pos(ax,idx_gps), t, "linear", "extrap");
        end
    end

    if any(sim.gps_vel_available)
        idx_vel = find(sim.gps_vel_available);

        for ax = 1:3
            sim.v_true(ax,:) = interp1(t(idx_vel), sim.gps_vel(ax,idx_vel), t, "linear", "extrap");
        end
    end

    if ~isempty(sim.q_ref)
        for j = 1:N
            [~, ia] = min(abs(sim.att_t - t(j)));
            sim.q_true(:,j) = sim.q_ref(:,ia);
        end
    else
        sim.q_true = repmat([1;0;0;0], 1, N);
    end

    %% ---------------- Unused fields expected by old code ----------------
    sim.range = nan(1,N);
    sim.range_available = false(1,N);

    %% ---------------- Summary ----------------
    fprintf("=== FIXED-WING SIM READY ===\n");
    fprintf("MAT file: %s\n", mat_file);
    fprintf("Original time span selected: %.2f to %.2f s\n", t_start, t_end);
    fprintf("Reset duration: %.2f s\n", t(end));
    fprintf("IMU samples: %d, approx rate %.2f Hz\n", N, 1/median(diff(t)));
    fprintf("GPS pos updates: %d\n", sum(sim.gps_pos_available));
    fprintf("GPS vel updates: %d\n", sum(sim.gps_vel_available));
    fprintf("Baro updates: %d\n", sum(sim.baro_available));
end