function [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt)
%COMPUTE_F_G_QD
% 18-state error-state EKF için F, G ve ayrık süreç gürültüsü Qd üretir.
%
% Error-state sırası:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%
% Nominal state:
%   state.p_n      : 3x1 position
%   state.v_n      : 3x1 velocity
%   state.q_nb     : 4x1 quaternion, body -> nav
%   state.b_g      : 3x1 gyro bias
%   state.b_a      : 3x1 accel bias
%   state.b_baro   : scalar barometer offset
%   state.wind_ne  : 2x1 horizontal wind estimate [wind_N; wind_E]
%
% Noise vector:
%   w = [n_g; n_a; n_bg; n_ba; n_bbaro; n_wind_N; n_wind_E]
%
% Boyutlar:
%   F  : 18x18
%   G  : 18x15
%   Qc : 15x15
%   Qd : 18x18
%
% Not:
%   Wind state deterministic olarak sabit kabul edilir:
%       wind_dot = noise
%   Yani rüzgarı asıl airspeed/TAS update'i gözlenebilir yapacak.
%   Burada sadece random-walk process modeli eklenir.

    %% 0) Input checks
    if nargin < 4
        error('compute_F_G_Qd requires state, imu, params, dt.');
    end

    if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error('dt pozitif ve sonlu scalar olmalıdır.');
    end

    %% 1) Required nominal quantities
    q_nb = state.q_nb(:);
    b_g  = state.b_g(:);
    b_a  = state.b_a(:);

    if numel(q_nb) ~= 4
        error('state.q_nb 4x1 quaternion olmalıdır.');
    end

    if numel(b_g) ~= 3 || numel(b_a) ~= 3
        error('state.b_g ve state.b_a 3x1 olmalıdır.');
    end

    R_nb = quat_to_rotmat(q_nb);

    omega_m = imu.gyro_m(:);
    accel_m = imu.accel_m(:);

    if numel(omega_m) ~= 3 || numel(accel_m) ~= 3
        error('imu.gyro_m ve imu.accel_m 3x1 olmalıdır.');
    end

    omega_c = omega_m - b_g;
    f_c     = accel_m - b_a;

    %% 2) Basic blocks
    I3 = eye(3);

    Omega_skew = skew3(omega_c);
    Fc_skew    = skew3(f_c);

    %% 3) Continuous-time F matrix
    %
    % dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
    %
    % p_dot      = v
    % v_dot      = R_nb * (accel_m - b_a) + g
    % theta_dot  = attitude error dynamics
    % b_g_dot    = gyro bias random walk
    % b_a_dot    = accel bias random walk
    % b_baro_dot = baro offset random walk
    % wind_dot   = horizontal wind random walk
    %
    % Deterministic F içinde:
    %   db_baro, dw_N, dw_E için ek coupling yoktur.
    %   Bu state'ler process noise ve airspeed/baro measurement update ile öğrenilir.

    F = zeros(18,18);

    % dp_dot = dv
    F(1:3, 4:6) = I3;

    % dv_dot error terms
    % a_n = R_nb * (accel_m - b_a) + g_n
    %
    % delta_v_dot wrt attitude error:
    %   -R_nb * skew(f_c) * dtheta
    %
    % delta_v_dot wrt accel bias:
    %   -R_nb * dba
    F(4:6, 7:9)   = -R_nb * Fc_skew;
    F(4:6, 13:15) = -R_nb;

    % attitude error dynamics
    F(7:9, 7:9)   = -Omega_skew;
    F(7:9, 10:12) = -I3;

    % Bias, baro offset and wind deterministic dynamics are zero.
    % F(10:18,:) remains zero except terms above.

    %% 4) Noise input matrix G
    %
    % Noise vector:
    %   w = [
    %       n_g        1:3
    %       n_a        4:6
    %       n_bg       7:9
    %       n_ba       10:12
    %       n_bbaro    13
    %       n_wind_NE  14:15
    %   ]

    G = zeros(18,15);

    % Gyro measurement noise affects attitude error
    G(7:9, 1:3) = -I3;

    % Accelerometer measurement noise affects velocity error
    G(4:6, 4:6) = -R_nb;

    % Gyro bias random walk
    G(10:12, 7:9) = I3;

    % Accelerometer bias random walk
    G(13:15, 10:12) = I3;

    % Barometer offset random walk
    G(16, 13) = 1;

    % Horizontal wind random walk
    G(17:18, 14:15) = eye(2);

    %% 5) Continuous-time noise covariance Qc

    if ~isfield(params, 'sigma_baro_bias_rw')
        params.sigma_baro_bias_rw = 0.01;   % [m/sqrt(s)]
    end

    if ~isfield(params, 'sigma_wind_rw')
        params.sigma_wind_rw = 0.05;        % [m/s/sqrt(s)]
    end

    Qc = diag([
        (params.sigma_g^2)           * ones(3,1);   % gyro white noise
        (params.sigma_a^2)           * ones(3,1);   % accel white noise
        (params.sigma_bg_rw^2)       * ones(3,1);   % gyro bias random walk
        (params.sigma_ba_rw^2)       * ones(3,1);   % accel bias random walk
        params.sigma_baro_bias_rw^2;                % baro offset random walk
        (params.sigma_wind_rw^2)     * ones(2,1)    % horizontal wind random walk
    ]);

    %% 6) Discrete process noise Qd
    %
    % Continuous noise mapped to state:
    %   Qx_c = G Qc G'
    %
    % Discrete approximation:
    %   Qd ≈ Qx_c dt + 0.5 * (F Qx_c + Qx_c F') dt^2

    Qx_c = G * Qc * G.';

    Qd = Qx_c * dt + 0.5 * (F * Qx_c + Qx_c * F.') * dt^2;

    Qd = 0.5 * (Qd + Qd.');
end