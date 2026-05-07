function [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt)
%COMPUTE_F_G_QD
% 16-state error-state EKF için F, G ve ayrık süreç gürültüsü Qd üretir.
%
% Error-state sırası:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro]
%
% Nominal state:
%   state.p_n     : 3x1 position
%   state.v_n     : 3x1 velocity
%   state.q_nb    : 4x1 quaternion, body -> nav
%   state.b_g     : 3x1 gyro bias
%   state.b_a     : 3x1 accel bias
%   state.b_baro  : scalar barometer offset
%
% Noise vector:
%   w = [n_g; n_a; n_bg; n_ba; n_bbaro]
%
% Boyutlar:
%   F  : 16x16
%   G  : 16x13
%   Qc : 13x13
%   Qd : 16x16

    %% 0) Input checks
    if nargin < 4
        error('compute_F_G_Qd requires state, imu, params, dt.');
    end

    if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error('dt pozitif ve sonlu scalar olmalıdır.');
    end

    %% 1) Gerekli büyüklükler
    q_nb = state.q_nb(:);
    b_g  = state.b_g(:);
    b_a  = state.b_a(:);

    R_nb = quat_to_rotmat(q_nb);

    omega_m = imu.gyro_m(:);
    accel_m = imu.accel_m(:);

    omega_c = omega_m - b_g;
    f_c     = accel_m - b_a;

    %% 2) Temel bloklar
    I3 = eye(3);

    Omega_skew = skew3(omega_c);
    Fc_skew    = skew3(f_c);

    %% 3) F matrisi
    %
    % dx = [dp; dv; dtheta; dbg; dba; db_baro]
    %
    % Baro offset nominal propagation:
    %   b_baro_dot = w_bbaro
    %
    % Bu yüzden deterministic F içinde 16. state için ek bağ yok:
    %   F(16,:) = 0
    %
    % Yani baro offset sabit kabul edilir, ama random walk noise ile
    % yavaş değişmesine izin verilir.

    F = zeros(16,16);

    % dp_dot = dv
    F(1:3, 4:6) = I3;

    % dv_dot error terms
    % a_n = R_nb * (accel_m - b_a) + g_n
    % delta_v_dot wrt dtheta and dba
    F(4:6, 7:9)   = -R_nb * Fc_skew;
    F(4:6, 13:15) = -R_nb;

    % attitude error dynamics
    F(7:9, 7:9)   = -Omega_skew;
    F(7:9, 10:12) = -I3;

    % gyro bias random walk deterministic part: zero
    % accel bias random walk deterministic part: zero
    % baro bias random walk deterministic part: zero

    %% 4) G matrisi
    %
    % w = [n_g; n_a; n_bg; n_ba; n_bbaro]
    %
    % n_g      : gyro white noise, 3x1
    % n_a      : accel white noise, 3x1
    % n_bg     : gyro bias random walk, 3x1
    % n_ba     : accel bias random walk, 3x1
    % n_bbaro  : baro offset random walk, 1x1

    G = zeros(16,13);

    % accel measurement noise affects velocity
    G(4:6, 4:6) = -R_nb;

    % gyro measurement noise affects attitude
    G(7:9, 1:3) = -I3;

    % gyro bias random walk
    G(10:12, 7:9) = I3;

    % accel bias random walk
    G(13:15, 10:12) = I3;

    % baro offset random walk
    G(16, 13) = 1;

    %% 5) Sürekli zaman noise covariance Qc

    if ~isfield(params, 'sigma_baro_bias_rw')
        params.sigma_baro_bias_rw = 0.01; % [m/sqrt(s)] default
    end

    Qc = diag([
        (params.sigma_g^2)  * ones(3,1);          % gyro white noise
        (params.sigma_a^2)  * ones(3,1);          % accel white noise
        (params.sigma_bg_rw^2) * ones(3,1);       % gyro bias RW
        (params.sigma_ba_rw^2) * ones(3,1);       % accel bias RW
        params.sigma_baro_bias_rw^2              % baro offset RW
    ]);

    %% 6) Ayrık süreç gürültüsü Qd
    %
    % Mevcut 15-state kodla aynı yaklaşım:
    %   Qx_c = G Qc G'
    %   Qd ≈ Qx_c dt + 0.5 (F Qx_c + Qx_c F') dt^2

    Qx_c = G * Qc * G.';

    Qd = Qx_c * dt + 0.5 * (F * Qx_c + Qx_c * F.') * dt^2;

    Qd = 0.5 * (Qd + Qd.');
end