function sim = simulate_truth_and_sensors(params)
%SIMULATE_TRUTH_AND_SENSORS
% Senaryo ayarlı truth ve sensör verisi üretir.
%
% Bu sürümde quaternion/gyro konvansiyonu:
% - q_nb : body -> nav
% - IMU gyro : body-frame angular velocity
% - dolayısıyla truth gyro üretiminde relative quaternion:
%       dq_b = q_{k-1}^* ⊗ q_k
%   kullanılır.

    %% Zaman ayarları
    if isfield(params,'sim') && isfield(params.sim,'dt')
        sim.dt = params.sim.dt;
    else
        sim.dt = 0.01;
    end

    if isfield(params,'sim') && isfield(params.sim,'T_end')
        T_end = params.sim.T_end;
    else
        T_end = 20;
    end

    sim.t = 0:sim.dt:T_end;

    t  = sim.t;
    dt = sim.dt;
    N  = numel(t);

    %% 1) Attitude truth
    [roll_true, pitch_true, yaw_true, q_true] = generate_attitude_truth(t, params);

    %% 2) Translasyon truth (body frame)
    [v_body, a_body_cmd] = generate_translation_truth_body(t, dt, params);

    %% 3) Body -> NED hız / konum
    p_true = zeros(3, N);
    v_true = zeros(3, N);

    for k = 1:N
        R_nb = quat_to_rotmat_local(q_true(:,k));
        v_true(:,k) = R_nb * v_body(:,k);

        if k > 1
            p_true(:,k) = p_true(:,k-1) + 0.5 * (v_true(:,k-1) + v_true(:,k)) * dt;
        end
    end

    %% 4) NED ivme
    a_true_n = zeros(3, N);

    for k = 2:N-1
        a_true_n(:,k) = (v_true(:,k+1) - v_true(:,k-1)) / (2*dt);
    end
    a_true_n(:,1) = (v_true(:,2) - v_true(:,1)) / dt;
    a_true_n(:,N) = (v_true(:,N) - v_true(:,N-1)) / dt;

    %% 5) Truth angular rate (BODY frame)
    omega_true = zeros(3, N);

    for k = 2:N
        % body-frame incremental rotation:
        % q_k = q_{k-1} ⊗ dq_b
        % => dq_b = q_{k-1}^* ⊗ q_k
        dq_b = quat_multiply_local(quat_conjugate_local(q_true(:,k-1)), q_true(:,k));

        rotvec_b = quat_to_rotvec_local(dq_b);
        omega_true(:,k) = rotvec_b / dt;
    end

    omega_true(:,1) = omega_true(:,2);

    %% 6) Gerçek biaslar (random walk)
    bg_true = zeros(3, N);
    ba_true = zeros(3, N);

    bg_true(:,1) = params.truth.bg0(:);
    ba_true(:,1) = params.truth.ba0(:);

    for k = 2:N
        bg_true(:,k) = bg_true(:,k-1) + params.sigma_bg_rw * sqrt(dt) * randn(3,1);
        ba_true(:,k) = ba_true(:,k-1) + params.sigma_ba_rw * sqrt(dt) * randn(3,1);
    end

    %% 7) IMU truth specific force ve ölçüm
    f_true_b  = zeros(3, N);
    imu_gyro  = zeros(3, N);
    imu_accel = zeros(3, N);

    for k = 1:N
        R_nb = quat_to_rotmat_local(q_true(:,k));
        R_bn = R_nb.';

        % a_n = R_nb * f_b + g_n
        % --> f_b = R_bn * (a_n - g_n)
        f_true_b(:,k) = R_bn * (a_true_n(:,k) - params.g_n);

        imu_gyro(:,k) = omega_true(:,k) ...
                      + bg_true(:,k) ...
                      + params.sigma_g * randn(3,1);

        imu_accel(:,k) = f_true_b(:,k) ...
                       + ba_true(:,k) ...
                       + params.sigma_a * randn(3,1);
    end

    %% 8) GNSS position (5 Hz)
    gps_pos_available = false(1, N);
    gps_pos = nan(3, N);

    gps_pos_step = round(0.20 / dt);

    for k = 1:gps_pos_step:N
        gps_pos_available(k) = true;
        gps_pos(:,k) = p_true(:,k) + params.sigma_gps_pos .* randn(3,1);
    end

    %% 9) GNSS velocity (5 Hz)
    gps_vel_available = false(1, N);
    gps_vel = nan(3, N);

    gps_vel_step = round(0.20 / dt);

    for k = 1:gps_vel_step:N
        gps_vel_available(k) = true;
        gps_vel(:,k) = v_true(:,k) + params.sigma_gps_vel .* randn(3,1);
    end

    %% 10) Baro (5 Hz)
    baro_available = false(1, N);
    baro = nan(1, N);

    baro_step = round(0.20 / dt);

    for k = 1:baro_step:N
        baro_available(k) = true;
        baro(k) = p_true(3,k) + params.sigma_baro * randn;
    end

    %% 11) Çıkış
    sim.p_true = p_true;
    sim.v_true = v_true;
    sim.q_true = q_true;

    sim.euler_true = [roll_true; pitch_true; yaw_true];
    sim.omega_true = omega_true;
    sim.a_true_n   = a_true_n;
    sim.f_true_b   = f_true_b;

    sim.v_body     = v_body;
    sim.a_body_cmd = a_body_cmd;

    sim.bg_true = bg_true;
    sim.ba_true = ba_true;

    sim.imu_gyro  = imu_gyro;
    sim.imu_accel = imu_accel;

    sim.gps_pos_available = gps_pos_available;
    sim.gps_pos = gps_pos;

    sim.gps_vel_available = gps_vel_available;
    sim.gps_vel = gps_vel;

    sim.baro_available = baro_available;
    sim.baro = baro;
end

%% =============================================================
% ATTITUDE SENARYO ÜRETİCİSİ
% =============================================================
function [roll_true, pitch_true, yaw_true, q_true] = generate_attitude_truth(t, params)

    N = numel(t);
    mode = params.truth.attitude_mode;

    switch mode
        case 'fixed'
            eul_deg = params.truth.attitude.fixed_euler_deg(:);
            roll_true  = deg2rad(eul_deg(1)) * ones(1,N);
            pitch_true = deg2rad(eul_deg(2)) * ones(1,N);
            yaw_true   = deg2rad(eul_deg(3)) * ones(1,N);

        case 'sinusoid'
            s = params.truth.attitude.sin;

            roll_true  = deg2rad(s.roll_amp_deg)  * sin(s.roll_freq  * t);
            pitch_true = deg2rad(s.pitch_amp_deg) * sin(s.pitch_freq * t);
            yaw_true   = deg2rad(s.yaw_bias_rate_deg) * t ...
                       + deg2rad(s.yaw_amp_deg) * sin(s.yaw_freq * t);

        case 'piecewise'
            tb = params.truth.attitude.piecewise.t_breaks(:).';
            E  = params.truth.attitude.piecewise.euler_deg;

            if size(E,1) ~= 3
                error('attitude.piecewise.euler_deg boyutu 3xM olmalıdır.');
            end

            if size(E,2) == numel(tb)
                roll_true  = deg2rad(interp1(tb, E(1,:), t, 'previous', 'extrap'));
                pitch_true = deg2rad(interp1(tb, E(2,:), t, 'previous', 'extrap'));
                yaw_true   = deg2rad(interp1(tb, E(3,:), t, 'previous', 'extrap'));

            elseif size(E,2) == numel(tb)-1
                roll_true  = piecewise_interval_series_deg(tb, E(1,:), t);
                pitch_true = piecewise_interval_series_deg(tb, E(2,:), t);
                yaw_true   = piecewise_interval_series_deg(tb, E(3,:), t);

            else
                error(['attitude.piecewise için boyut uyumsuz. ', ...
                       'Ya numel(t_breaks)=size(euler_deg,2) olmalı ', ...
                       '(node-based), ya da numel(t_breaks)=size(euler_deg,2)+1 olmalı ', ...
                       '(interval-based).']);
            end

        otherwise
            error('Bilinmeyen attitude_mode: %s', mode);
    end

    q_true = zeros(4, N);
    for k = 1:N
        q_true(:,k) = euler321_to_quat_local([roll_true(k); pitch_true(k); yaw_true(k)]);
    end
end

%% =============================================================
% TRANSLASYON SENARYO ÜRETİCİSİ
% =============================================================
function [v_body, a_body_cmd] = generate_translation_truth_body(t, dt, params)

    N = numel(t);
    mode = params.truth.translation_mode;

    v_body = zeros(3, N);
    a_body_cmd = zeros(3, N);

    v_body(:,1) = params.truth.translation.v0_body(:);

    switch mode
        case 'constant_velocity'
            v_const = params.truth.translation.const_vel_body(:);
            v_body = repmat(v_const, 1, N);
            a_body_cmd(:,:) = 0;

        case 'constant_accel'
            a_const = params.truth.translation.const_acc_body(:);

            for k = 2:N
                a_body_cmd(:,k-1) = a_const;
                v_body(:,k) = v_body(:,k-1) + a_const * dt;
            end
            a_body_cmd(:,N) = a_const;

        case 'random_smooth'
            r = params.truth.translation.random;

            if isfield(r, 'seed') && ~isempty(r.seed)
                rng(r.seed);
            end

            nom_v = r.nom_v_body(:);
            tau_a = r.tau_a(:);
            sigma = r.sigma_a_rand(:);
            tau_v = r.tau_v(:);

            a_rand = zeros(3,N);

            for k = 2:N
                w = randn(3,1);
                a_rand(:,k) = a_rand(:,k-1) ...
                            + dt * (-a_rand(:,k-1) ./ tau_a) ...
                            + sqrt(dt) * sigma .* w;
            end

            a_body_cmd = a_rand;
            v_body(:,1) = nom_v;

            for k = 2:N
                dv = a_body_cmd(:,k-1) + (nom_v - v_body(:,k-1)) ./ tau_v;
                v_body(:,k) = v_body(:,k-1) + dv * dt;
            end

        case 'piecewise_accel'
            tb = params.truth.translation.piecewise.t_breaks(:).';
            A  = params.truth.translation.piecewise.acc_body;

            if size(A,1) ~= 3
                error('translation.piecewise.acc_body boyutu 3xM olmalıdır.');
            end

            if size(A,2) ~= numel(tb)-1
                error('translation.piecewise için size(acc_body,2) = numel(t_breaks)-1 olmalıdır.');
            end

            for k = 2:N
                tk = t(k-1);
                idx = find(tk >= tb(1:end-1) & tk < tb(2:end), 1, 'last');

                if isempty(idx)
                    idx = size(A,2);
                end

                a_now = A(:,idx);
                a_body_cmd(:,k-1) = a_now;
                v_body(:,k) = v_body(:,k-1) + a_now * dt;
            end
            a_body_cmd(:,N) = a_body_cmd(:,N-1);

        otherwise
            error('Bilinmeyen translation_mode: %s', mode);
    end
end

%% =============================================================
% YARDIMCI FONKSİYONLAR
% =============================================================
function series_rad = piecewise_interval_series_deg(tb, values_deg, t)
    N = numel(t);
    series_rad = zeros(1, N);

    for k = 1:N
        tk = t(k);
        idx = find(tk >= tb(1:end-1) & tk < tb(2:end), 1, 'last');

        if isempty(idx)
            idx = numel(values_deg);
        end

        series_rad(k) = deg2rad(values_deg(idx));
    end
end

function q = euler321_to_quat_local(eul)
    roll  = eul(1);
    pitch = eul(2);
    yaw   = eul(3);

    cr = cos(roll/2);  sr = sin(roll/2);
    cp = cos(pitch/2); sp = sin(pitch/2);
    cy = cos(yaw/2);   sy = sin(yaw/2);

    q = [ cr*cp*cy + sr*sp*sy;
          sr*cp*cy - cr*sp*sy;
          cr*sp*cy + sr*cp*sy;
          cr*cp*sy - sr*sp*cy ];

    q = quat_normalize_local(q);
end

function qc = quat_conjugate_local(q)
    qc = [q(1); -q(2:4)];
end

function q = quat_normalize_local(q)
    q = q(:) / norm(q);
end

function q = quat_multiply_local(q1, q2)
    w1 = q1(1); x1 = q1(2); y1 = q1(3); z1 = q1(4);
    w2 = q2(1); x2 = q2(2); y2 = q2(3); z2 = q2(4);

    q = [ w1*w2 - x1*x2 - y1*y2 - z1*z2;
          w1*x2 + x1*w2 + y1*z2 - z1*y2;
          w1*y2 - x1*z2 + y1*w2 + z1*x2;
          w1*z2 + x1*y2 - y1*x2 + z1*w2 ];

    q = quat_normalize_local(q);
end

function R = quat_to_rotmat_local(q)
    q = quat_normalize_local(q);

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    R = [1 - 2*(qy^2 + qz^2),   2*(qx*qy - qz*qw),   2*(qx*qz + qy*qw);
         2*(qx*qy + qz*qw),     1 - 2*(qx^2 + qz^2), 2*(qy*qz - qx*qw);
         2*(qx*qz - qy*qw),     2*(qy*qz + qx*qw),   1 - 2*(qx^2 + qy^2)];
end

function rotvec = quat_to_rotvec_local(q)
    q = quat_normalize_local(q);

    if q(1) < 0
        q = -q;
    end

    v = q(2:4);
    nv = norm(v);

    if nv < 1e-12
        rotvec = 2 * v;
        return;
    end

    angle = 2 * atan2(nv, q(1));
    axis  = v / nv;
    rotvec = angle * axis;
end