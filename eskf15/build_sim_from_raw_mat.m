function sim = build_sim_from_raw_mat(filename)
%BUILD_SIM_FROM_RAW_MAT
% MRS rectangle_raw_clean.mat dosyasını mevcut 15-state ESKF sim formatına çevirir.
%
% Bu sürümde:
%   - IMU raw olarak kullanılır.
%   - GPS/RTK raw LLA'dan üretilmiş gps_pos_enu ölçümü kullanılır.
%   - GPS/RTK pozisyonundan smoothing + finite difference ile pseudo velocity üretilir.
%   - Baro raw pressure'dan üretilmiş baro_alt ölçümü kullanılır.
%   - Truth/reference yalnızca kıyas ve plot için kullanılır.
%
% Convention:
%   Nav frame: ENU
%   p = [E; N; U]
%   v = [v_E; v_N; v_U]
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
    % IMU timeline master olduğu için t0 = t_imu(1) seçiyoruz.
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

    %% 7) Raw GPS/RTK position measurement channel
    % Burada GPS'i interpolate edip her adımda update yapmıyoruz.
    % Sadece gerçek GPS timestamp'ine en yakın IMU indexinde update açıyoruz.
    sim.gps_pos = nan(3,N);
    sim.gps_pos_available = false(1,N);

    for i = 1:numel(t_gps)
        if ~isfinite(t_gps(i)) || any(~isfinite(gps_pos_enu(:,i)))
            continue;
        end

        [~, idx] = min(abs(t - t_gps(i)));

        if idx >= 1 && idx <= N
            sim.gps_pos(:,idx) = gps_pos_enu(:,i);
            sim.gps_pos_available(idx) = true;
        end
    end

    %% 7.5) Pseudo GPS velocity from raw RTK/GPS position
    % /rtk_raw doğrudan velocity vermiyor.
    % Bu yüzden raw RTK/GPS pozisyonundan:
    %
    %   1) NaN/Inf temizliği
    %   2) Duplicate timestamp temizliği
    %   3) Moving average smoothing
    %   4) Central finite difference / gradient
    %
    % ile pseudo GPS velocity üretiyoruz.
    %
    % Dikkat:
    % Bu gerçek Doppler GPS velocity değildir.
    % Pozisyon türevinden elde edilen yardımcı velocity ölçümüdür.
    % Bu yüzden main içinde R_gps_vel fazla küçük seçilmemeli.
    
    sim.gps_vel = nan(3,N);
    sim.gps_vel_available = false(1,N);

    valid_gps = isfinite(t_gps) & all(isfinite(gps_pos_enu),1);

    t_gps_valid = t_gps(valid_gps);
    gps_pos_valid = gps_pos_enu(:,valid_gps);

    % Sadece IMU zaman aralığına düşen GPS örneklerini kullan.
    valid_time = (t_gps_valid >= t(1)) & (t_gps_valid <= t(end));

    t_gps_valid = t_gps_valid(valid_time);
    gps_pos_valid = gps_pos_valid(:,valid_time);

    % Duplicate timestamp temizliği.
    [t_gps_u, ia] = unique(t_gps_valid, 'stable');
    gps_pos_u = gps_pos_valid(:,ia);

    gps_vel_enu = [];
    gps_pos_smooth = [];

    if numel(t_gps_u) >= 7

        % GPS yaklaşık 10 Hz.
        % 7 örnek yaklaşık 0.7 saniyelik smoothing demek.
        % Çok küçük olursa velocity gürültülü olur.
        % Çok büyük olursa velocity gecikir.
        win = 7;

        gps_pos_smooth = gps_pos_u;

        for ax = 1:3
            gps_pos_smooth(ax,:) = movmean(gps_pos_u(ax,:), win);
        end

        gps_vel_enu = zeros(3, numel(t_gps_u));

        % Central finite difference.
        % gradient(Y, X) nonuniform timestamp için de çalışır.
        for ax = 1:3
            gps_vel_enu(ax,:) = gradient(gps_pos_smooth(ax,:), t_gps_u);
        end

        % İlk ve son birkaç türev örneği güvenilmez.
        valid_vel = true(1, numel(t_gps_u));
        valid_vel(1:3) = false;
        valid_vel(end-2:end) = false;

        % Fiziksel olmayan hız zıplamalarını ele.
        % Rectangle uçuşu için 5 m/s üstü çoğu zaman türev gürültüsüdür.
        speed = vecnorm(gps_vel_enu, 2, 1);
        valid_vel = valid_vel & isfinite(speed) & speed < 5.0;

        for i = 1:numel(t_gps_u)
            if ~valid_vel(i)
                continue;
            end

            [~, idx] = min(abs(t - t_gps_u(i)));

            if idx >= 1 && idx <= N
                sim.gps_vel(:,idx) = gps_vel_enu(:,i);
                sim.gps_vel_available(idx) = true;
            end
        end

        sim.gps_vel_source = 'pseudo velocity from smoothed raw RTK/GPS position';
        sim.gps_vel_raw_t = t_gps_u;
        sim.gps_vel_raw = gps_vel_enu;
        sim.gps_pos_smooth_for_vel = gps_pos_smooth;
        sim.gps_vel_smoothing_window = win;

    else
        sim.gps_vel_source = 'not available';
        sim.gps_vel_raw_t = [];
        sim.gps_vel_raw = [];
        sim.gps_pos_smooth_for_vel = [];
        sim.gps_vel_smoothing_window = nan;
    end

    %% 8) Raw baro measurement channel
    sim.baro = nan(1,N);
    sim.baro_available = false(1,N);

    for i = 1:numel(t_baro)
        if ~isfinite(t_baro(i)) || ~isfinite(baro_alt(i))
            continue;
        end

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
    sim.gps_vel_source = sim.gps_vel_source;
    sim.truth_source = '/rtk_fused_odom reference only';

    if isfield(data, 'gps_lla')
        sim.gps_lla = data.gps_lla;
    end

    if isfield(data, 'gps_fix_type')
        sim.gps_fix_type = data.gps_fix_type;
    end

    if isfield(data, 'lat0')
        sim.lat0 = data.lat0;
        sim.lon0 = data.lon0;
        sim.alt0 = data.alt0;
    end

    if isfield(data, 'baro_pressure')
        sim.baro_pressure = data.baro_pressure;
    end

    if isfield(data, 'baro_p0')
        sim.baro_p0 = data.baro_p0;
    end

    %% 11) Debug info
    disp("=== RAW SIM READY ===")
    disp(["IMU rate:", num2str(1/sim.dt), "Hz"])
    disp(["GPS raw position updates:", num2str(sum(sim.gps_pos_available))])
    disp(["GPS pseudo velocity updates:", num2str(sum(sim.gps_vel_available))])
    disp(["BARO raw updates:", num2str(sum(sim.baro_available))])

    if ~isempty(gps_vel_enu)
        speed = vecnorm(gps_vel_enu, 2, 1);
        disp(["Pseudo GPS velocity mean speed:", num2str(mean(speed, 'omitnan')), "m/s"])
        disp(["Pseudo GPS velocity max speed:", num2str(max(speed, [], 'omitnan')), "m/s"])
    end
end