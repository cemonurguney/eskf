function sim = build_sim_from_raw_mat(filename)
%BUILD_SIM_FROM_RAW_MAT
% MRS rectangle_raw_clean.mat dosyasını mevcut 15-state ESKF sim formatına çevirir.
%
% Bu sürümde:
%   - IMU raw olarak kullanılır.
%   - GPS/RTK raw LLA'dan üretilmiş gps_pos_enu ölçümü kullanılır.
%   - Baro raw pressure'dan üretilmiş baro_alt ölçümü kullanılır.
%   - Truth/reference yalnızca kıyas ve plot için kullanılır.
%
% Convention:
%   Nav frame: ENU
%   p = [E; N; U]
%   g_n = [0; 0; -9.81] main içinde set edilmeli.

    %% 1) Load
    data = load(filename);

    required_fields = { ...
        't_imu', 'imu_gyro', 'imu_accel', ...
        't_gps', 'gps_pos_enu', ...
        't_baro', 'baro_alt', ...
        't_truth', 'p_truth_rel', 'v_truth', 'q_truth' ...
    };

    for i = 1:numel(required_fields)
        if ~isfield(data, required_fields{i})
            error('MAT dosyasinda eksik alan: %s', required_fields{i});
        end
    end

    %% 2) Raw arrays
    t_imu = data.t_imu(:)';
    imu_gyro = data.imu_gyro;
    imu_accel = data.imu_accel;

    t_gps = data.t_gps(:)';
    gps_pos_enu = data.gps_pos_enu;

    t_baro = data.t_baro(:)';
    baro_alt = data.baro_alt(:)';

    t_truth = data.t_truth(:)';
    p_truth = data.p_truth_rel;
    v_truth = data.v_truth;
    q_truth = data.q_truth;

    %% 3) Dimension sanity
    if size(imu_gyro,1) ~= 3
        imu_gyro = imu_gyro';
    end

    if size(imu_accel,1) ~= 3
        imu_accel = imu_accel';
    end

    if size(gps_pos_enu,1) ~= 3
        gps_pos_enu = gps_pos_enu';
    end

    if size(p_truth,1) ~= 3
        p_truth = p_truth';
    end

    if size(v_truth,1) ~= 3
        v_truth = v_truth';
    end

    if size(q_truth,1) ~= 4
        q_truth = q_truth';
    end

    %% 4) Time normalize
    % Colab extraction zaten normalize ediyor ama burada yine güvene alıyoruz.
    t0 = t_imu(1);

    t_imu = t_imu - t0;
    t_gps = t_gps - t0;
    t_baro = t_baro - t0;
    t_truth = t_truth - t0;

    %% 5) IMU timeline = master timeline
    t = t_imu;
    N = numel(t);

    sim.t = t;
    sim.dt = mean(diff(t));

    sim.imu_gyro = imu_gyro(:,1:N);
    sim.imu_accel = imu_accel(:,1:N);

    %% 6) Truth/reference interpolation to IMU timeline
    % Truth sadece kıyas için kullanılır, filtreye measurement olarak verilmez.
    sim.p_true = interp1(t_truth', p_truth', t, 'linear', 'extrap')';
    sim.v_true = interp1(t_truth', v_truth', t, 'linear', 'extrap')';

    q_interp = interp1(t_truth', q_truth', t, 'linear', 'extrap')';

    for k = 1:N
        nq = norm(q_interp(:,k));
        if nq < 1e-12
            q_interp(:,k) = [1;0;0;0];
        else
            q_interp(:,k) = q_interp(:,k) / nq;
        end
    end

    sim.q_true = q_interp;

    %% 7) Raw GPS/RTK measurement channel
    % Burada GPS'i interpolate edip her adımda update yapmıyoruz.
    % Sadece gerçek GPS timestamp'ine en yakın IMU indexinde update açıyoruz.
    sim.gps_pos = nan(3,N);
    sim.gps_pos_available = false(1,N);

    for i = 1:numel(t_gps)
        [~, idx] = min(abs(t - t_gps(i)));

        if idx >= 1 && idx <= N
            sim.gps_pos(:,idx) = gps_pos_enu(:,i);
            sim.gps_pos_available(idx) = true;
        end
    end

    % Raw /rtk_raw velocity vermediği için GPS velocity update kapalı.
    sim.gps_vel = nan(3,N);
    sim.gps_vel_available = false(1,N);

    %% 8) Raw baro measurement channel
    sim.baro = nan(1,N);
    sim.baro_available = false(1,N);

    for i = 1:numel(t_baro)
        [~, idx] = min(abs(t - t_baro(i)));

        if idx >= 1 && idx <= N
            sim.baro(idx) = baro_alt(i);
            sim.baro_available(idx) = true;
        end
    end

    %% 9) Range placeholder
    sim.range = zeros(1,N);
    sim.range_available = false(1,N);

    %% 10) Metadata
    sim.is_raw_mrs = true;
    sim.nav_frame = 'ENU';
    sim.gps_source = 'raw /rtk_raw converted LLA -> ENU';
    sim.truth_source = '/rtk_fused_odom reference only';

    if isfield(data, 'gps_lla')
        sim.gps_lla = data.gps_lla;
    end

    if isfield(data, 'lat0')
        sim.lat0 = data.lat0;
        sim.lon0 = data.lon0;
        sim.alt0 = data.alt0;
    end

    %% 11) Debug info
    disp("=== RAW SIM READY ===")
    disp(["IMU rate:", num2str(1/sim.dt), "Hz"])
    disp(["GPS raw updates:", num2str(sum(sim.gps_pos_available))])
    disp(["BARO raw updates:", num2str(sum(sim.baro_available))])
    disp(["GPS velocity updates:", num2str(sum(sim.gps_vel_available)), "(disabled)"])
end