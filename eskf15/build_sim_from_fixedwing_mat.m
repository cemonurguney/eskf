function sim = build_sim_from_fixedwing_mat(mat_file, t_start, t_end)
%BUILD_SIM_FROM_FIXEDWING_MAT
% Fixed-wing PX4 ULog'dan uretilmis *_fixedwing_eskf.mat dosyasini,
% mevcut ESKF main kodunun bekledigi sim struct formatina cevirir.
%
% Frame convention:
%   NED
%   p = [N; E; D]
%   v = [Vn; Ve; Vd]
%   g_n = [0; 0; +9.81]
%
% Baro convention:
%   sim.baro = baro_z_down_m
%   16/18-state modelde:
%       z_baro = p_D + b_baro + noise
%
% Airspeed:
%   Logged TAS kullanilir.
%   TAS yoksa IAS fallback kullanilir.
%
% Wind:
%   Logged PX4 wind estimate update olarak kullanilmaz.
%   Sadece sim.wind_ref ile karsilastirma referansi olarak saklanir.
%
% Usage:
%   sim = build_sim_from_fixedwing_mat("0002_20.48_fixedwing_eskf.mat", 600, 975);

    if nargin < 2 || isempty(t_start)
        t_start = -inf;
    end
    if nargin < 3 || isempty(t_end)
        t_end = inf;
    end

    fw = load(mat_file);

    required_fields = {"imu_t", "gyro_rad_s", "acc_m_s2", ...
                       "gps_t", "gps_pos_ned_m", "gps_vel_ned_m_s", ...
                       "baro_t", "baro_z_down_m"};

    for i = 1:numel(required_fields)
        if ~isfield(fw, required_fields{i})
            error("MAT file icinde gerekli alan yok: %s", required_fields{i});
        end
    end

    %% ============================================================
    % 1) IMU trim and cleanup
    % ============================================================
    imu_t = fw.imu_t(:);
    gyro_rad_s = fw.gyro_rad_s;
    acc_m_s2   = fw.acc_m_s2;

    valid_imu = isfinite(imu_t) & ...
        all(isfinite(gyro_rad_s), 2) & ...
        all(isfinite(acc_m_s2), 2) & ...
        imu_t >= t_start & imu_t <= t_end;

    imu_t = imu_t(valid_imu);
    gyro_rad_s = gyro_rad_s(valid_imu, :);
    acc_m_s2   = acc_m_s2(valid_imu, :);

    if numel(imu_t) < 5
        error("Trim sonrasi IMU verisi bos veya cok az kaldi. t_start/t_end araligini kontrol et.");
    end

    [imu_t, sort_idx] = sort(imu_t);
    gyro_rad_s = gyro_rad_s(sort_idx, :);
    acc_m_s2   = acc_m_s2(sort_idx, :);

    %% ============================================================
    % 2) Reset time to zero
    % ============================================================
    t0 = imu_t(1);
    imu_t = imu_t - t0;
    t = imu_t(:).';
    N = numel(t);

    dt_all = diff(t);
    dt_nom = median(dt_all(dt_all > 0), "omitnan");
    if ~isfinite(dt_nom) || dt_nom <= 0
        error("IMU timestamp diff gecersiz gorunuyor.");
    end

    sim.t = t;
    sim.dt = dt_nom;
    sim.dt_vec = [dt_nom, diff(t)];

    sim.imu_gyro  = gyro_rad_s.';  % 3xN
    sim.imu_accel = acc_m_s2.';    % 3xN

    %% ============================================================
    % 3) GPS position / velocity mapped to nearest IMU sample
    % ============================================================
    gps_t = fw.gps_t(:);
    gps_pos = fw.gps_pos_ned_m;
    gps_vel = fw.gps_vel_ned_m_s;

    valid_gps_pos = isfinite(gps_t) & ...
        gps_t >= t_start & gps_t <= t_end & ...
        all(isfinite(gps_pos), 2);

    gps_t_abs = gps_t(valid_gps_pos);
    gps_t_rel = gps_t_abs - t0;
    gps_pos = gps_pos(valid_gps_pos, :);
    gps_vel = gps_vel(valid_gps_pos, :);

    in_range = gps_t_rel >= t(1) & gps_t_rel <= t(end);
    gps_t_rel = gps_t_rel(in_range);
    gps_pos = gps_pos(in_range, :);
    gps_vel = gps_vel(in_range, :);

    sim.gps_pos = nan(3, N);
    sim.gps_vel = nan(3, N);
    sim.gps_pos_available = false(1, N);
    sim.gps_vel_available = false(1, N);

    if ~isempty(gps_t_rel)
        gps_k = nearest_time_index(t, gps_t_rel);

        for i = 1:numel(gps_k)
            k = gps_k(i);
            sim.gps_pos(:, k) = gps_pos(i, :).';
            sim.gps_pos_available(k) = true;

            if all(isfinite(gps_vel(i, :)))
                sim.gps_vel(:, k) = gps_vel(i, :).';
                sim.gps_vel_available(k) = true;
            end
        end
    end

    %% ============================================================
    % 4) Barometer mapped to nearest IMU sample
    % ============================================================
    baro_t = fw.baro_t(:);
    baro_z_down_m = fw.baro_z_down_m(:);

    valid_baro = isfinite(baro_t) & isfinite(baro_z_down_m) & ...
        baro_t >= t_start & baro_t <= t_end;

    baro_t_abs = baro_t(valid_baro);
    baro_t_rel = baro_t_abs - t0;
    baro_z_down_m = baro_z_down_m(valid_baro);

    in_range = baro_t_rel >= t(1) & baro_t_rel <= t(end);
    baro_t_rel = baro_t_rel(in_range);
    baro_z_down_m = baro_z_down_m(in_range);

    sim.baro = nan(1, N);
    sim.baro_available = false(1, N);

    if ~isempty(baro_t_rel)
        baro_k = nearest_time_index(t, baro_t_rel);

        for i = 1:numel(baro_k)
            k = baro_k(i);
            sim.baro(k) = baro_z_down_m(i);
            sim.baro_available(k) = true;
        end
    end

    %% ============================================================
    % 5) Attitude reference/debug, not used as measurement update
    % ============================================================
    sim.att_t = [];
    sim.q_ref = [];
    sim.rpy_ref = [];

    if isfield(fw, "att_t") && isfield(fw, "q_nb") && ~isempty(fw.att_t) && ~isempty(fw.q_nb)
        att_t_abs = fw.att_t(:);
        q_nb = fw.q_nb;

        valid_att = isfinite(att_t_abs) & ...
            att_t_abs >= t_start & att_t_abs <= t_end & ...
            all(isfinite(q_nb), 2);

        att_t_rel = att_t_abs(valid_att) - t0;
        q_nb = q_nb(valid_att, :);

        in_range = att_t_rel >= t(1) & att_t_rel <= t(end);
        sim.att_t = att_t_rel(in_range).';
        sim.q_ref = q_nb(in_range, :).';  % 4xM

        for i = 1:size(sim.q_ref, 2)
            nq = norm(sim.q_ref(:, i));
            if nq > 0
                sim.q_ref(:, i) = sim.q_ref(:, i) / nq;
            end
        end

        if isfield(fw, "rpy_rad") && ~isempty(fw.rpy_rad)
            rpy = fw.rpy_rad(valid_att, :);
            rpy = rpy(in_range, :);
            sim.rpy_ref = rpy.'; % 3xM
        end
    end

    %% ============================================================
    % 6) Logged airspeed mapped to nearest IMU sample
    % ============================================================
    sim.airspeed = nan(1, N);
    sim.airspeed_available = false(1, N);
    sim.airspeed_source = "none";

    sim.airspeed_t = [];
    sim.ias_m_s = [];
    sim.tas_m_s = [];

    if isfield(fw, "airspeed_t") && ~isempty(fw.airspeed_t)

        air_t_abs_all = fw.airspeed_t(:);

        tas_all = [];
        ias_all = [];

        if isfield(fw, "tas_m_s") && ~isempty(fw.tas_m_s)
            tas_all = fw.tas_m_s(:);
        end

        if isfield(fw, "ias_m_s") && ~isempty(fw.ias_m_s)
            ias_all = fw.ias_m_s(:);
        end

        air_data_all = [];

        if ~isempty(tas_all) && numel(tas_all) == numel(air_t_abs_all)
            air_data_all = tas_all;
            sim.airspeed_source = "TAS";
        elseif ~isempty(ias_all) && numel(ias_all) == numel(air_t_abs_all)
            air_data_all = ias_all;
            sim.airspeed_source = "IAS";
        end

        if ~isempty(air_data_all)
            valid_air_time = isfinite(air_t_abs_all) & ...
                air_t_abs_all >= t_start & air_t_abs_all <= t_end;

            air_t_abs = air_t_abs_all(valid_air_time);
            air_t_rel = air_t_abs - t0;
            air_data = air_data_all(valid_air_time);

            if ~isempty(tas_all) && numel(tas_all) == numel(air_t_abs_all)
                tas_trim = tas_all(valid_air_time);
            else
                tas_trim = nan(size(air_t_abs));
            end

            if ~isempty(ias_all) && numel(ias_all) == numel(air_t_abs_all)
                ias_trim = ias_all(valid_air_time);
            else
                ias_trim = nan(size(air_t_abs));
            end

            in_range = air_t_rel >= t(1) & air_t_rel <= t(end);
            air_t_rel = air_t_rel(in_range);
            air_data = air_data(in_range);
            tas_trim = tas_trim(in_range);
            ias_trim = ias_trim(in_range);

            valid_air = isfinite(air_t_rel) & ...
                isfinite(air_data) & ...
                air_data > 3.0 & air_data < 80.0;

            sim.airspeed_t = air_t_rel(valid_air).';
            sim.airspeed = nan(1, N);
            sim.airspeed_available = false(1, N);

            if any(valid_air)
                air_t_use = air_t_rel(valid_air);
                air_data_use = air_data(valid_air);

                air_k = nearest_time_index(t, air_t_use);

                for i = 1:numel(air_k)
                    k = air_k(i);
                    sim.airspeed(k) = air_data_use(i);
                    sim.airspeed_available(k) = true;
                end

                sim.tas_m_s = tas_trim(valid_air).';
                sim.ias_m_s = ias_trim(valid_air).';
            end
        end
    end

    %% ============================================================
    % 7) PX4 wind estimate mapped to nearest IMU sample, comparison only
    % ============================================================
    sim.wind_ref = nan(2, N);
    sim.wind_ref_available = false(1, N);

    sim.wind_t = [];
    sim.wind_ned_m_s = [];

    if isfield(fw, "wind_t") && ~isempty(fw.wind_t) && ...
            isfield(fw, "wind_ned_m_s") && ~isempty(fw.wind_ned_m_s)

        wind_t_abs_all = fw.wind_t(:);
        wind_data = fw.wind_ned_m_s;

        % Accept both Mx3 and 3xM formats.
        if size(wind_data, 2) == 3
            wind_data = wind_data.';     % 3xM
        elseif size(wind_data, 1) == 3
            % already 3xM
        else
            error("wind_ned_m_s boyutu 3xM veya Mx3 olmalıdır.");
        end

        if size(wind_data, 2) ~= numel(wind_t_abs_all)
            warning("wind_t ve wind_ned_m_s uzunluklari uyusmuyor. Wind ref map edilmeyecek.");
        else
            valid_wind_time = isfinite(wind_t_abs_all) & ...
                wind_t_abs_all >= t_start & wind_t_abs_all <= t_end & ...
                all(isfinite(wind_data), 1).';

            wind_t_abs = wind_t_abs_all(valid_wind_time);
            wind_t_rel = wind_t_abs - t0;
            wind_data = wind_data(:, valid_wind_time);

            in_range = wind_t_rel >= t(1) & wind_t_rel <= t(end);
            wind_t_rel = wind_t_rel(in_range);
            wind_data = wind_data(:, in_range);

            sim.wind_t = wind_t_rel.';
            sim.wind_ned_m_s = wind_data;

            if ~isempty(wind_t_rel)
                wind_k = nearest_time_index(t, wind_t_rel);

                for i = 1:numel(wind_k)
                    k = wind_k(i);
                    sim.wind_ref(:, k) = wind_data(1:2, i);
                    sim.wind_ref_available(k) = true;
                end
            end
        end
    end

    %% ============================================================
    % 8) Reference-like arrays for plotting/debug only
    % ============================================================
    sim.p_true = nan(3, N);
    sim.v_true = nan(3, N);
    sim.q_true = nan(4, N);

    idx_gps = find(sim.gps_pos_available);
    if numel(idx_gps) >= 2
        [tu, ia] = unique(t(idx_gps), "stable");
        idxu = idx_gps(ia);
        for ax = 1:3
            sim.p_true(ax, :) = interp1(tu, sim.gps_pos(ax, idxu), t, "linear", "extrap");
        end
    end

    idx_vel = find(sim.gps_vel_available);
    if numel(idx_vel) >= 2
        [tu, ia] = unique(t(idx_vel), "stable");
        idxu = idx_vel(ia);
        for ax = 1:3
            sim.v_true(ax, :) = interp1(tu, sim.gps_vel(ax, idxu), t, "linear", "extrap");
        end
    end

    if ~isempty(sim.q_ref) && ~isempty(sim.att_t)
        qk = nearest_time_index(sim.att_t, t(:));
        sim.q_true = sim.q_ref(:, qk);
    else
        sim.q_true = repmat([1; 0; 0; 0], 1, N);
    end

    %% ============================================================
    % 9) Unused fields expected by older code
    % ============================================================
    sim.range = nan(1, N);
    sim.range_available = false(1, N);

    %% ============================================================
    % 10) Summary
    % ============================================================
    fprintf("=== FIXED-WING SIM READY ===\n");
    fprintf("MAT file: %s\n", mat_file);
    fprintf("Requested absolute time window: %.3f to %.3f s\n", t_start, t_end);
    fprintf("Reset duration: %.3f s\n", t(end));
    fprintf("IMU samples: %d, approx rate %.2f Hz\n", N, 1 / median(diff(t), "omitnan"));
    fprintf("GPS pos updates: %d\n", sum(sim.gps_pos_available));
    fprintf("GPS vel updates: %d\n", sum(sim.gps_vel_available));
    fprintf("Baro updates: %d\n", sum(sim.baro_available & isfinite(sim.baro)));
    fprintf("Attitude ref samples: %d\n", size(sim.q_ref, 2));
    fprintf("Airspeed source: %s\n", sim.airspeed_source);
    fprintf("Airspeed mapped samples: %d\n", sum(sim.airspeed_available));

    if any(sim.airspeed_available)
        fprintf("Airspeed mean/min/max [m/s]: %.3f / %.3f / %.3f\n", ...
            mean(sim.airspeed(sim.airspeed_available), "omitnan"), ...
            min(sim.airspeed(sim.airspeed_available)), ...
            max(sim.airspeed(sim.airspeed_available)));
    end

    fprintf("PX4 wind ref mapped samples: %d\n", sum(sim.wind_ref_available));
end

function idx = nearest_time_index(t_grid, t_query)
%NEAREST_TIME_INDEX Fast nearest index mapping for monotonically increasing time vectors.
    t_grid = t_grid(:);
    t_query = t_query(:);

    if isempty(t_grid)
        idx = [];
        return;
    end

    idx = interp1(t_grid, (1:numel(t_grid)).', t_query, "nearest", "extrap");
    idx = round(idx);
    idx = max(1, min(numel(t_grid), idx));
end