function sim = simulate_truth_and_sensors(params)
%SIMULATE_TRUTH_AND_SENSORS
% Senaryo ayarlı truth ve sensör verisi üretir.
%
% Quaternion / gyro konvansiyonu:
% - q_nb : body -> nav
% - IMU gyro : body-frame angular velocity
% - truth gyro üretiminde:
%       dq_b = q_{k-1}^* ⊗ q_k
%
% Desteklenen senaryolar:
% - legacy attitude/translation modes
% - fixedwing_like scenario
%
% Çıktı alanları:
%   sim.p_true, sim.v_true, sim.q_true
%   sim.imu_gyro, sim.imu_accel
%   sim.gps_pos, sim.gps_vel, sim.baro
%   sim.wind_ref, sim.wind_ref_available
%   sim.airspeed, sim.airspeed_available

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
t = sim.t;
dt = sim.dt;
N = numel(t);

%% ============================================================
% 1-5) Truth generation
% ============================================================

wind_ref_from_truth = [];

scenario_mode = "legacy";
if isfield(params, 'truth') && isfield(params.truth, 'scenario_mode')
    scenario_mode = string(params.truth.scenario_mode);
end

if scenario_mode == "fixedwing_like"

    [roll_true, pitch_true, yaw_true, q_true, ...
     p_true, v_true, a_true_n, omega_true, ...
     v_body, a_body_cmd, wind_ref_from_truth] = ...
        generate_fixedwing_like_truth(t, dt, params);

else

    %% 1) Attitude truth
    [roll_true, pitch_true, yaw_true, q_true] = generate_attitude_truth(t, params);

    %% 2) Translasyon truth, body frame
    [v_body, a_body_cmd] = generate_translation_truth_body(t, dt, params);

    %% 3) Body -> nav hız / konum
    p_true = zeros(3, N);
    v_true = zeros(3, N);

    for k = 1:N
        R_nb = quat_to_rotmat_local(q_true(:,k));
        v_true(:,k) = R_nb * v_body(:,k);

        if k > 1
            p_true(:,k) = p_true(:,k-1) + 0.5 * (v_true(:,k-1) + v_true(:,k)) * dt;
        end
    end

    %% 4) Nav-frame acceleration truth
    a_true_n = zeros(3, N);

    for k = 2:N-1
        a_true_n(:,k) = (v_true(:,k+1) - v_true(:,k-1)) / (2*dt);
    end

    a_true_n(:,1) = (v_true(:,2) - v_true(:,1)) / dt;
    a_true_n(:,N) = (v_true(:,N) - v_true(:,N-1)) / dt;

    %% 5) Truth angular rate, body frame
    omega_true = zeros(3, N);

    for k = 2:N
        dq_b = quat_multiply_local(quat_conjugate_local(q_true(:,k-1)), q_true(:,k));
        rotvec_b = quat_to_rotvec_local(dq_b);
        omega_true(:,k) = rotvec_b / dt;
    end

    omega_true(:,1) = omega_true(:,2);
end

%% ============================================================
% 6) True IMU biases, random walk
% ============================================================

bg_true = zeros(3, N);
ba_true = zeros(3, N);

if isfield(params, 'truth') && isfield(params.truth, 'bg0')
    bg_true(:,1) = params.truth.bg0(:);
else
    bg_true(:,1) = [0;0;0];
end

if isfield(params, 'truth') && isfield(params.truth, 'ba0')
    ba_true(:,1) = params.truth.ba0(:);
else
    ba_true(:,1) = [0;0;0];
end

for k = 2:N
    bg_true(:,k) = bg_true(:,k-1) + params.sigma_bg_rw * sqrt(dt) * randn(3,1);
    ba_true(:,k) = ba_true(:,k-1) + params.sigma_ba_rw * sqrt(dt) * randn(3,1);
end

%% ============================================================
% 7) IMU truth specific force and measurement
% ============================================================

f_true_b = zeros(3, N);
imu_gyro = zeros(3, N);
imu_accel = zeros(3, N);

for k = 1:N
    R_nb = quat_to_rotmat_local(q_true(:,k));
    R_bn = R_nb.';

    % a_n = R_nb * f_b + g_n
    % f_b = R_bn * (a_n - g_n)
    f_true_b(:,k) = R_bn * (a_true_n(:,k) - params.g_n);

    imu_gyro(:,k) = omega_true(:,k) + bg_true(:,k) + params.sigma_g * randn(3,1);
    imu_accel(:,k) = f_true_b(:,k) + ba_true(:,k) + params.sigma_a * randn(3,1);
end

%% ============================================================
% 8) GNSS position
% ============================================================

gps_pos_available = false(1, N);
gps_pos = nan(3, N);

if isfield(params, 'synthetic_gps_rate_hz')
    gps_rate_hz = params.synthetic_gps_rate_hz;
else
    gps_rate_hz = 5;
end

gps_pos_step = max(1, round((1/gps_rate_hz) / dt));

for k = 1:gps_pos_step:N
    gps_pos_available(k) = true;
    gps_pos(:,k) = p_true(:,k) + params.sigma_gps_pos .* randn(3,1);
end

%% ============================================================
% 9) GNSS velocity
% ============================================================

gps_vel_available = false(1, N);
gps_vel = nan(3, N);

gps_vel_step = gps_pos_step;

for k = 1:gps_vel_step:N
    gps_vel_available(k) = true;
    gps_vel(:,k) = v_true(:,k) + params.sigma_gps_vel .* randn(3,1);
end

%% ============================================================
% 10) Baro measurement
% ============================================================

baro_available = false(1, N);
baro = nan(1, N);

if isfield(params, 'synthetic_baro_rate_hz')
    baro_rate_hz = params.synthetic_baro_rate_hz;
else
    baro_rate_hz = 5;
end

baro_step = max(1, round((1/baro_rate_hz) / dt));

if isfield(params, 'truth') && isfield(params.truth, 'b_baro0')
    b_baro0_true = params.truth.b_baro0;
else
    b_baro0_true = 0;
end

if isfield(params, 'truth') && isfield(params.truth, 'b_baro_rw')
    b_baro_rw = params.truth.b_baro_rw;
else
    b_baro_rw = 0;
end

b_baro_true = zeros(1,N);
b_baro_true(1) = b_baro0_true;

for k = 2:N
    b_baro_true(k) = b_baro_true(k-1) + b_baro_rw * sqrt(dt) * randn;
end

for k = 1:baro_step:N
    baro_available(k) = true;
    baro(k) = p_true(3,k) + b_baro_true(k) + params.sigma_baro * randn;
end

%% ============================================================
% 11) Synthetic horizontal wind and TAS
% ============================================================

wind_ref = nan(2,N);
wind_ref_available = false(1,N);

airspeed = nan(1,N);
airspeed_available = false(1,N);

enable_wind_tas = false;

if isfield(params, 'synthetic') && isfield(params.synthetic, 'enable_wind_tas')
    enable_wind_tas = params.synthetic.enable_wind_tas;
end

if enable_wind_tas

    % Fixed-wing-like truth trajectory içinde wind kullandıysa aynı wind'i kullan.
    % Legacy senaryoda burada wind üret.
    if ~isempty(wind_ref_from_truth)
        wind_ref = wind_ref_from_truth;
        wN = wind_ref(1,:);
        wE = wind_ref(2,:);
    else
        if isfield(params.synthetic, 'wind_seed') && ~isempty(params.synthetic.wind_seed)
            rng(params.synthetic.wind_seed);
        end

        if isfield(params.synthetic, 'wind_mode')
            wind_mode = params.synthetic.wind_mode;
        else
            wind_mode = "slow_sine";
        end

        switch string(wind_mode)
            case "constant"
                if isfield(params.synthetic, 'wind_constant_ne')
                    w_const = params.synthetic.wind_constant_ne(:);
                else
                    w_const = [2.0; -1.2];
                end

                wN = w_const(1) * ones(1,N);
                wE = w_const(2) * ones(1,N);

            case "slow_sine"
                if isfield(params.synthetic, 'wind_base_ne')
                    w_base = params.synthetic.wind_base_ne(:);
                else
                    w_base = [2.0; -1.2];
                end

                if isfield(params.synthetic, 'wind_amp_ne')
                    w_amp = params.synthetic.wind_amp_ne(:);
                else
                    w_amp = [0.4; 0.3];
                end

                if isfield(params.synthetic, 'wind_freq_ne')
                    w_freq = params.synthetic.wind_freq_ne(:);
                else
                    w_freq = [0.030; 0.025];
                end

                wN = w_base(1) + w_amp(1) * sin(w_freq(1) * t);
                wE = w_base(2) + w_amp(2) * cos(w_freq(2) * t + 0.4);

            otherwise
                error('Unknown synthetic wind_mode: %s', wind_mode);
        end

        wind_ref = [wN; wE];
    end

    wind_ref_available = true(1,N);

    wind_3d = [wind_ref(1,:); wind_ref(2,:); zeros(1,N)];
    v_air = v_true - wind_3d;
    tas_true = vecnorm(v_air, 2, 1);

    if isfield(params, 'sigma_tas')
        sigma_tas = params.sigma_tas;
    else
        sigma_tas = 0.30;
    end

    tas_meas_full = tas_true + sigma_tas * randn(1,N);

    if isfield(params, 'synthetic_airspeed_rate_hz')
        tas_rate_hz = params.synthetic_airspeed_rate_hz;
    else
        tas_rate_hz = 20;
    end

    tas_step = max(1, round((1/tas_rate_hz) / dt));

    for k = 1:tas_step:N
        airspeed_available(k) = true;
        airspeed(k) = tas_meas_full(k);
    end

    fprintf('[simulate] Synthetic wind/TAS enabled.\n');
    fprintf('[simulate] Wind N mean/min/max [m/s]: %.3f / %.3f / %.3f\n', ...
        mean(wN), min(wN), max(wN));
    fprintf('[simulate] Wind E mean/min/max [m/s]: %.3f / %.3f / %.3f\n', ...
        mean(wE), min(wE), max(wE));
    fprintf('[simulate] TAS samples: %d\n', sum(airspeed_available));
    fprintf('[simulate] TAS mean/min/max [m/s]: %.3f / %.3f / %.3f\n', ...
        mean(airspeed(airspeed_available), "omitnan"), ...
        min(airspeed(airspeed_available)), ...
        max(airspeed(airspeed_available)));
end

%% ============================================================
% 12) Output
% ============================================================

sim.p_true = p_true;
sim.v_true = v_true;
sim.q_true = q_true;
sim.euler_true = [roll_true; pitch_true; yaw_true];

sim.omega_true = omega_true;
sim.a_true_n = a_true_n;
sim.f_true_b = f_true_b;

sim.v_body = v_body;
sim.a_body_cmd = a_body_cmd;

sim.bg_true = bg_true;
sim.ba_true = ba_true;
sim.b_baro_true = b_baro_true;

sim.imu_gyro = imu_gyro;
sim.imu_accel = imu_accel;

sim.gps_pos_available = gps_pos_available;
sim.gps_pos = gps_pos;

sim.gps_vel_available = gps_vel_available;
sim.gps_vel = gps_vel;

sim.baro_available = baro_available;
sim.baro = baro;

sim.wind_ref = wind_ref;
sim.wind_ref_available = wind_ref_available;

sim.airspeed = airspeed;
sim.airspeed_available = airspeed_available;
end

%% =============================================================
% FIXED-WING-LIKE TRUTH GENERATOR
% =============================================================
function [roll_true, pitch_true, yaw_true, q_true, ...
          p_true, v_true, a_true_n, omega_true, ...
          v_body, a_body_cmd, wind_ref] = ...
          generate_fixedwing_like_truth(t, dt, params)
%GENERATE_FIXEDWING_LIKE_TRUTH
% Fixed-wing benzeri sentetik truth üretir.
%
% NED frame:
%   p = [N; E; D]
%   D positive down
%
% Mantık:
%   - Önce air-relative velocity tanımlanır.
%   - Wind eklenerek ground velocity oluşturulur.
%   - Position ground velocity ile entegre edilir.
%   - Attitude yaw air-relative heading'e yakın tutulur.
%   - Roll coordinated-turn benzeri yaw-rate'ten türetilir.
%   - Pitch küçük climb/descent açısıdır.

    N = numel(t);

    fw = params.truth.fixedwing;

    if isfield(fw, 'p0_ned')
        p0 = fw.p0_ned(:);
    else
        p0 = [0;0;-30];
    end

    Va_trim = getfield_with_default_local(fw, 'Va_trim', 14.0);
    Va_amp1 = getfield_with_default_local(fw, 'Va_amp1', 1.8);
    Va_amp2 = getfield_with_default_local(fw, 'Va_amp2', 0.7);

    yaw0 = deg2rad(getfield_with_default_local(fw, 'yaw0_deg', 30.0));
    yaw_rate_base = deg2rad(getfield_with_default_local(fw, 'yaw_rate_base_deg_s', 0.35));

    yaw_amp1 = deg2rad(getfield_with_default_local(fw, 'yaw_amp1_deg', 22.0));
    yaw_freq1 = getfield_with_default_local(fw, 'yaw_freq1', 0.020);

    yaw_amp2 = deg2rad(getfield_with_default_local(fw, 'yaw_amp2_deg', 7.0));
    yaw_freq2 = getfield_with_default_local(fw, 'yaw_freq2', 0.060);

    gamma_amp = deg2rad(getfield_with_default_local(fw, 'gamma_amp_deg', 3.0));
    gamma_freq = getfield_with_default_local(fw, 'gamma_freq', 0.018);

    max_roll = deg2rad(getfield_with_default_local(fw, 'max_roll_deg', 35.0));

    %% Wind profile, same source as synthetic TAS
    if isfield(params, 'synthetic') && isfield(params.synthetic, 'wind_base_ne')
        w_base = params.synthetic.wind_base_ne(:);
    else
        w_base = [2.0; -1.2];
    end

    if isfield(params, 'synthetic') && isfield(params.synthetic, 'wind_amp_ne')
        w_amp = params.synthetic.wind_amp_ne(:);
    else
        w_amp = [0.4; 0.3];
    end

    if isfield(params, 'synthetic') && isfield(params.synthetic, 'wind_freq_ne')
        w_freq = params.synthetic.wind_freq_ne(:);
    else
        w_freq = [0.030; 0.025];
    end

    wN = w_base(1) + w_amp(1) * sin(w_freq(1) * t);
    wE = w_base(2) + w_amp(2) * cos(w_freq(2) * t + 0.4);

    wind_ref = [wN; wE];

    %% Air-relative speed and attitude path
    Va = Va_trim ...
        + Va_amp1 * sin(0.014 * t + 0.2) ...
        + Va_amp2 * sin(0.057 * t + 1.1);

    Va = max(Va, 7.0);

    yaw_true = yaw0 ...
        + yaw_rate_base * t ...
        + yaw_amp1 * sin(yaw_freq1 * t) ...
        + yaw_amp2 * sin(yaw_freq2 * t + 0.8);

    yaw_rate = gradient(yaw_true, dt);

    % gamma > 0 means climbing, so D velocity is negative.
    gamma = gamma_amp * sin(gamma_freq * t + 0.5);

    Va_h = Va .* cos(gamma);

    v_air_N = Va_h .* cos(yaw_true);
    v_air_E = Va_h .* sin(yaw_true);
    v_air_D = -Va .* sin(gamma);

    %% Ground velocity = air-relative velocity + wind
    v_true = zeros(3,N);
    v_true(1,:) = v_air_N + wN;
    v_true(2,:) = v_air_E + wE;
    v_true(3,:) = v_air_D;

    %% Position integration
    p_true = zeros(3,N);
    p_true(:,1) = p0;

    for k = 2:N
        p_true(:,k) = p_true(:,k-1) + 0.5 * (v_true(:,k-1) + v_true(:,k)) * dt;
    end

    %% Attitude
    g_abs = abs(params.g_n(3));

    roll_true = atan2(Va_h .* yaw_rate, g_abs);
    roll_true = min(max(roll_true, -max_roll), max_roll);

    % Small pitch following climb angle plus slight oscillation.
    pitch_true = gamma + deg2rad(1.0) * sin(0.033 * t + 0.3);

    q_true = zeros(4,N);

    for k = 1:N
        q_true(:,k) = euler321_to_quat_local([roll_true(k); pitch_true(k); yaw_true(k)]);
    end

    %% Nav-frame acceleration truth
    a_true_n = zeros(3,N);

    for k = 2:N-1
        a_true_n(:,k) = (v_true(:,k+1) - v_true(:,k-1)) / (2*dt);
    end

    a_true_n(:,1) = (v_true(:,2) - v_true(:,1)) / dt;
    a_true_n(:,N) = (v_true(:,N) - v_true(:,N-1)) / dt;

    %% Body-frame velocity and approximate body acceleration command
    v_body = zeros(3,N);
    a_body_cmd = zeros(3,N);

    for k = 1:N
        R_nb = quat_to_rotmat_local(q_true(:,k));
        R_bn = R_nb.';
        v_body(:,k) = R_bn * v_true(:,k);
        a_body_cmd(:,k) = R_bn * a_true_n(:,k);
    end

    %% Body angular rate from quaternion difference
    omega_true = zeros(3,N);

    for k = 2:N
        dq_b = quat_multiply_local(quat_conjugate_local(q_true(:,k-1)), q_true(:,k));
        rotvec_b = quat_to_rotvec_local(dq_b);
        omega_true(:,k) = rotvec_b / dt;
    end

    omega_true(:,1) = omega_true(:,2);

    fprintf('[simulate] Fixed-wing-like truth enabled.\n');
    fprintf('[simulate] Va true mean/min/max [m/s]: %.3f / %.3f / %.3f\n', ...
        mean(Va), min(Va), max(Va));
    fprintf('[simulate] Ground speed mean/min/max [m/s]: %.3f / %.3f / %.3f\n', ...
        mean(vecnorm(v_true,2,1)), min(vecnorm(v_true,2,1)), max(vecnorm(v_true,2,1)));
    fprintf('[simulate] Roll mean/min/max [deg]: %.3f / %.3f / %.3f\n', ...
        mean(rad2deg(roll_true)), min(rad2deg(roll_true)), max(rad2deg(roll_true)));
end

%% =============================================================
% ATTITUDE SENARYO ÜRETİCİSİ, LEGACY
% =============================================================
function [roll_true, pitch_true, yaw_true, q_true] = generate_attitude_truth(t, params)
    N = numel(t);
    mode = params.truth.attitude_mode;

    switch mode
        case 'fixed'
            eul_deg = params.truth.attitude.fixed_euler_deg(:);
            roll_true = deg2rad(eul_deg(1)) * ones(1,N);
            pitch_true = deg2rad(eul_deg(2)) * ones(1,N);
            yaw_true = deg2rad(eul_deg(3)) * ones(1,N);

        case 'sinusoid'
            s = params.truth.attitude.sin;

            roll_true = deg2rad(s.roll_amp_deg) * sin(s.roll_freq * t);
            pitch_true = deg2rad(s.pitch_amp_deg) * sin(s.pitch_freq * t);
            yaw_true = deg2rad(s.yaw_bias_rate_deg) * t + ...
                       deg2rad(s.yaw_amp_deg) * sin(s.yaw_freq * t);

        case 'piecewise'
            tb = params.truth.attitude.piecewise.t_breaks(:).';
            E = params.truth.attitude.piecewise.euler_deg;

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
                error('attitude.piecewise boyut uyumsuz.');
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
% TRANSLASYON SENARYO ÜRETİCİSİ, LEGACY
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
            A = params.truth.translation.piecewise.acc_body;

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
function val = getfield_with_default_local(s, field_name, default_val)
    if isfield(s, field_name)
        val = s.(field_name);
    else
        val = default_val;
    end
end

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
    roll = eul(1);
    pitch = eul(2);
    yaw = eul(3);

    cr = cos(roll/2); sr = sin(roll/2);
    cp = cos(pitch/2); sp = sin(pitch/2);
    cy = cos(yaw/2); sy = sin(yaw/2);

    q = [
        cr*cp*cy + sr*sp*sy;
        sr*cp*cy - cr*sp*sy;
        cr*sp*cy + sr*cp*sy;
        cr*cp*sy - sr*sp*cy
    ];

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

    q = [
        w1*w2 - x1*x2 - y1*y2 - z1*z2;
        w1*x2 + x1*w2 + y1*z2 - z1*y2;
        w1*y2 - x1*z2 + y1*w2 + z1*x2;
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ];

    q = quat_normalize_local(q);
end

function R = quat_to_rotmat_local(q)
    q = quat_normalize_local(q);

    qw = q(1); qx = q(2); qy = q(3); qz = q(4);

    R = [
        1 - 2*(qy^2 + qz^2), 2*(qx*qy - qz*qw), 2*(qx*qz + qy*qw);
        2*(qx*qy + qz*qw), 1 - 2*(qx^2 + qz^2), 2*(qy*qz - qx*qw);
        2*(qx*qz - qy*qw), 2*(qy*qz + qx*qw), 1 - 2*(qx^2 + qy^2)
    ];
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
    axis = v / nv;
    rotvec = angle * axis;
end