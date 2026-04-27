function [state, P, residual, K, dx_hat] = update_gnss_pos(state, P, z_gps_pos, params)
%UPDATE_GNSS_POS
% GNSS konum ölçümü ile error-state EKF update yapar.
%
% Girdi:
%   state     : nominal state struct
%   P         : 15x15 error-state covariance
%   z_gps_pos : 3x1 GNSS NED position measurement [m]
%   params    : filtre parametreleri (R_gps_pos ve joseph form ayarı içerir)
%
% Çıktı:
%   state     : update edilmiş nominal state
%   P         : update edilmiş covariance
%   residual  : 3x1 innovation / residual
%   K         : 15x3 Kalman gain
%   dx_hat    : 15x1 error-state correction
%
% Ölçüm modeli:
%   z = p_n + v
%
% Error-state lineer hali:
%   r = z - p_nom ≈ H * delta_x + v
%
% Burada:
%   H = [I 0 0 0 0]

    %% =========================================================
    % 1) Boyut kontrolleri
    % ==========================================================
    z_gps_pos = z_gps_pos(:);

    if numel(z_gps_pos) ~= 3
        error('z_gps_pos, 3x1 GNSS position measurement olmalıdır.');
    end

    if ~isequal(size(P), [15 15])
        error('P matrisi 15x15 olmalıdır.');
    end

    %% =========================================================
    % 2) Beklenen ölçüm ve residual
    % ==========================================================
    % GNSS position doğrudan nominal konumu ölçüyor kabul ediliyor
    z_hat = state.p_n;                 % h(x_nom)
    residual = z_gps_pos - z_hat;      % innovation

    %% =========================================================
    % 3) H matrisi
    % ==========================================================
    % Error-state sırası:
    % [dp dv dtheta dbg dba]
    %
    % GNSS position yalnızca delta p bloğunu doğrudan görür
    H = zeros(3,15);
    H(:,1:3) = eye(3);

    %% =========================================================
    % 4) Innovation covariance
    % ==========================================================
    S = H * P * H.' + params.R_gps_pos;

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
        % Joseph form daha numerik kararlı
        A = I15 - K * H;
        P = A * P * A.' + K * params.R_gps_pos * K.';
    else
        % Klasik form
        P = (I15 - K * H) * P;
    end

    % Sayısal simetri koruması
    P = 0.5 * (P + P.');

end