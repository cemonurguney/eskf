function sim = build_sim_from_mat(filename)
%BUILD_SIM_FROM_MAT
% hover_clean.mat gibi gercek veriyi mevcut ESKF sim formatina cevirir.
%
% Beklenen MAT alanlari:
%   t_imu, imu_gyro, imu_accel
%   t_baro, baro
%   t_range, range
%   t_odom, p_true, v_true, q_true
%
% Cikti:
%   sim struct:
%       t, dt
%       imu_gyro, imu_accel
%       gps_pos, gps_vel
%       baro
%       range
%       gps_pos_available, gps_vel_available, baro_available, range_available
%       p_true, v_true, q_true
%
% CONVENTION:
%   Nav frame ENU (z up positive) — main_eskf_15state icinde params.g_n bu
%   dosyaya uyacak sekilde [0;0;-9.81] olarak set edilmeli.

    %% =========================
    % 1) Load safely
    %% =========================
    data = load(filename);%s

    % Zorunlu alanlar
    required_fields = {'t_imu','imu_gyro','imu_accel','t_baro','baro','t_odom','p_true','v_true','q_true'};
    for i = 1:numel(required_fields)
        if ~isfield(data, required_fields{i})
            error('MAT dosyasinda eksik alan: %s', required_fields{i});
        end
    end

    % IMU
    t_imu     = data.t_imu(:)';
    imu_gyro  = data.imu_gyro;
    imu_accel = data.imu_accel;

    % BARO
    t_baro = data.t_baro(:)';
    baro   = data.baro(:)';   % pressure [Pa]

    % RANGE opsiyonel
    has_range = isfield(data,'t_range') && (isfield(data,'range') || isfield(data,'range_data'));
    if has_range
        t_range = data.t_range(:)';
        if isfield(data,'range')
            range_data = data.range(:)';
        else
            range_data = data.range_data(:)';
        end
    else
        t_range = [];
        range_data = [];
    end

    % ODOM / truth
    t_odom = data.t_odom(:)';
    p_true = data.p_true;
    v_true = data.v_true;
    q_true = data.q_true;

    %% =========================
    % 2) Normalize times
    %% =========================
    t0 = t_imu(1);

    t_imu  = t_imu  - t0;
    t_baro = t_baro - t0;
    t_odom = t_odom - t0;

    if has_range
        t_range = t_range - t0;
    end

    %% =========================
    % 3) IMU timeline = master timeline
    %% =========================
    t = t_imu;
    N = numel(t);

    sim.t = t;
    sim.dt = mean(diff(t));

    sim.imu_gyro  = imu_gyro;
    sim.imu_accel = imu_accel;

    %% =========================
    % 4) Interpolate truth/odom to IMU timeline
    %% =========================
    p_interp = interp1(t_odom', p_true', t, 'linear', 'extrap')';
    v_interp = interp1(t_odom', v_true', t, 'linear', 'extrap')';
    q_interp = interp1(t_odom', q_true', t, 'linear', 'extrap')';

    % Quaternion normalize
    for k = 1:N
        nq = norm(q_interp(:,k));
        if nq < 1e-12
            q_interp(:,k) = [1;0;0;0];
        else
            q_interp(:,k) = q_interp(:,k) / nq;
        end
    end

    sim.p_true = p_interp;
    sim.v_true = v_interp;
    sim.q_true = q_interp;

    %% =========================
    % 5) GPS measurement channels
    %% =========================
    % Ilk asamada RTK fused odom'u hem measurement hem truth gibi kullaniyoruz.
    % Bu "pipeline calisiyor mu?" testi icin yeterli.
    sim.gps_pos = p_interp;
    sim.gps_vel = v_interp;

    gps_mask = false(1,N);

    % Odom zamanindaki en yakin IMU indexlerine measurement isareti koy
    for i = 1:numel(t_odom)
        [~, idx] = min(abs(t - t_odom(i)));
        gps_mask(idx) = true;
    end

    sim.gps_pos_available = gps_mask;
    sim.gps_vel_available = gps_mask;

    %% =========================
    % 6) BARO: pressure -> relative altitude -> ENU up
    %% =========================
    % Standart atmosfer yaklasimi:
    %   h_up = 44330 * (1 - (P / Pref)^(1/5.255))
    %
    % Pref = ilk birkac orneginin ortalamasi
    %
    % Nav frame ENU oldugu icin z = up positive, yani baro = +alt_up.
    % (NED olsaydi baro = -alt_up olacakti.)
    %
    n_ref = min(50, numel(baro));
    P_ref = mean(baro(1:n_ref));

    alt_up = 44330 * (1 - (baro / P_ref).^(1/5.255));   % up-positive relative altitude [m]
    baro_z = alt_up;                                    % ENU: z up positive

    baro_interp = interp1(t_baro', baro_z', t, 'linear', 'extrap')';
    sim.baro = baro_interp;

    baro_mask = false(1,N);
    for i = 1:numel(t_baro)
        [~, idx] = min(abs(t - t_baro(i)));
        baro_mask(idx) = true;
    end
    sim.baro_available = baro_mask;

    %% =========================
    % 7) RANGE optional
    %% =========================
    if has_range && ~isempty(t_range) && ~isempty(range_data)
        range_interp = interp1(t_range', range_data', t, 'linear', 'extrap')';
        sim.range = range_interp;

        range_mask = false(1,N);
        for i = 1:numel(t_range)
            [~, idx] = min(abs(t - t_range(i)));
            range_mask(idx) = true;
        end
        sim.range_available = range_mask;
    else
        sim.range = zeros(1,N);
        sim.range_available = false(1,N);
    end

    %% =========================
    % 8) Debug info
    %% =========================
    disp("=== SIM READY ===")
    disp(["IMU rate:", num2str(1/sim.dt), "Hz"])
    disp(["GPS updates:", num2str(sum(sim.gps_pos_available))])
    disp(["BARO updates:", num2str(sum(sim.baro_available))])
    disp(["RANGE updates:", num2str(sum(sim.range_available))])

end