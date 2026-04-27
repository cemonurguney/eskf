function state = propagate_nominal(state, imu, params, dt)
%PROPAGATE_NOMINAL
% 15-state error-state EKF için nominal state propagation.
%
% Bu sürümde kullanılan konvansiyon:
% - q_nb : body -> nav dönüşünü temsil eder
% - gyro ölçümü body frame'dedir
% - bu yüzden quaternion artımı sağdan uygulanır:
%       q_k+1 = q_k ⊗ dq_b

    %% 1) IMU ölçümleri
    omega_m = imu.gyro_m(:);
    accel_m = imu.accel_m(:);

    %% 2) Bias düzeltmesi
    omega_b = omega_m - state.b_g;   % body-frame angular rate
    f_b     = accel_m - state.b_a;   % body-frame specific force

    %% 3) Quaternion güncelle (RIGHT multiplication)
    dtheta_b = omega_b * dt;
    dq_b = small_angle_to_quat(dtheta_b);

    state.q_nb = quat_multiply(state.q_nb, dq_b);
    state.q_nb = quat_normalize(state.q_nb);

    %% 4) Rotation matrix
    R_nb = quat_to_rotmat(state.q_nb);

    %% 5) Navigation-frame ivme
    % a_n = R_nb * f_b + g_n
    a_n = R_nb * f_b + params.g_n;

    %% 6) Position + velocity update
    v_old = state.v_n;

    state.p_n = state.p_n + v_old * dt + 0.5 * a_n * dt^2;
    state.v_n = v_old + a_n * dt;
end