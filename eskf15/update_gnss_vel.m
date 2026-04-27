function [state, P, residual, K, dx_hat] = update_gnss_vel(state, P, z_gps_vel, params)
%UPDATE_GNSS_VEL
% GNSS hız ölçümü ile error-state EKF update yapar.
%
% Girdi:
%   state     : nominal state struct
%   P         : 15x15 error-state covariance
%   z_gps_vel : 3x1 GNSS NED velocity measurement [m/s]
%   params    : filtre parametreleri (R_gps_vel ve joseph form ayarı içerir)
%
% Çıktı:
%   state     : update edilmiş nominal state
%   P         : update edilmiş covariance
%   residual  : 3x1 innovation / residual
%   K         : 15x3 Kalman gain
%   dx_hat    : 15x1 error-state correction
%
% Ölçüm modeli:
%   z = v_n + noise
%
% Error-state lineer hali:
%   r = z - v_nom ≈ H * delta_x + noise
%
% Burada:
%   H = [0 I 0 0 0]

    %% =========================================================
    % 1) Boyut kontrolleri
    % ==========================================================
    z_gps_vel = z_gps_vel(:);

    if numel(z_gps_vel) ~= 3
        error('z_gps_vel, 3x1 GNSS velocity measurement olmalıdır.');
    end

    if ~isequal(size(P), [15 15])
        error('P matrisi 15x15 olmalıdır.');
    end

    %% =========================================================
    % 2) Beklenen ölçüm ve residual
    % ==========================================================
    % GNSS velocity doğrudan nominal hızı ölçüyor kabul ediliyor
    z_hat = state.v_n;                 % h(x_nom)
    residual = z_gps_vel - z_hat;      % innovation

    %% =========================================================
    % 3) H matrisi
    % ==========================================================
    % Error-state sırası:
    % [dp dv dtheta dbg dba]
    %
    % GNSS velocity yalnızca delta v bloğunu doğrudan görür
    H = zeros(3,15);
    H(:,4:6) = eye(3);

    %% =========================================================
    % 4) Innovation covariance
    % ==========================================================
    S = H * P * H.' + params.R_gps_vel;

    %% =========================================================
    % 5) Kalman gain
    % ==========================================================
    K = P * H.' / S;

    %% =========================================================
    % 6) Error-state düzeltmesi
    % ==========================================================
    dx_hat = K * residual;

    %% =========================================================
    % 7) Nominal state'e uygula
    % ==========================================================
    state = inject_error_state(state, dx_hat);

    %% =========================================================
    % 8) Covariance update
    % ==========================================================
    I15 = eye(15);

    if isfield(params, 'use_joseph_form') && params.use_joseph_form
        A = I15 - K * H;
        P = A * P * A.' + K * params.R_gps_vel * K.';
    else
        P = (I15 - K * H) * P;
    end

    % Sayısal simetri koruması
    P = 0.5 * (P + P.');

end