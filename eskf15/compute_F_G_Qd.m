function [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt)
%COMPUTE_F_G_QD
% 15-state error-state EKF için F, G ve ayrık süreç gürültüsü Qd üretir.

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
    % delta x = [dp; dv; dtheta; dbg; dba]
    F = zeros(15,15);

    F(1:3, 4:6) = I3;

    F(4:6, 7:9)   = -R_nb * Fc_skew;
    F(4:6, 13:15) = -R_nb;

    F(7:9, 7:9)   = -Omega_skew;
    F(7:9, 10:12) = -I3;

    %% 4) G matrisi
    % w = [n_g; n_a; n_bg; n_ba]
    G = zeros(15,12);

    G(4:6, 4:6) = -R_nb;
    G(7:9, 1:3) = -I3;
    G(10:12, 7:9) = I3;
    G(13:15, 10:12) = I3;

    %% 5) Sürekli zaman noise covariance
    Qc = diag([
        (params.sigma_g^2)     * ones(3,1);
        (params.sigma_a^2)     * ones(3,1);
        (params.sigma_bg_rw^2) * ones(3,1);
        (params.sigma_ba_rw^2) * ones(3,1)
    ]);

    %% 6) Ayrık süreç gürültüsü Qd
    Qx_c = G * Qc * G.';
    Qd = Qx_c * dt + 0.5 * (F * Qx_c + Qx_c * F.') * dt^2;

    Qd = 0.5 * (Qd + Qd.');
end