function [state, P, residual, K, dx_hat] = update_baro(state, P, z_baro, params)
%UPDATE_BARO
% Barometre yükseklik/irtifa ölçümü ile error-state EKF update yapar.
%
% Girdi:
%   state  : nominal state struct
%   P      : 15x15 error-state covariance
%   z_baro : 1x1 veya scalar baro ölçümü
%   params : filtre parametreleri (R_baro ve joseph form ayarı içerir)
%
% Çıktı:
%   state    : update edilmiş nominal state
%   P        : update edilmiş covariance
%   residual : scalar innovation / residual
%   K        : 15x1 Kalman gain
%   dx_hat   : 15x1 error-state correction
%
% Ölçüm modeli:
%   z = p_D + noise
%
% Error-state lineer hali:
%   r = z - p_D_nom ≈ H * delta_x + noise
%
% Burada:
%   H = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
%
% Not:
% Bu model, baro ölçümünün NED "Down" ekseni ile uyumlu hale getirildiğini varsayar.
% Eğer elindeki baro "yükseklik yukarı pozitif" veriyorsa, işaret dönüşümünü
% bu fonksiyona gelmeden önce yapman gerekir.

    %% =========================================================
    % 1) Giriş düzenleme
    % ==========================================================
    z_baro = z_baro(1);

    if ~isequal(size(P), [15 15])
        error('P matrisi 15x15 olmalıdır.');
    end

    %% =========================================================
    % 2) Beklenen ölçüm ve residual
    % ==========================================================
    % Baro sadece nominal position'ın down bileşenini görüyor
    z_hat = state.p_n(3);        % p_D_nom
    residual = z_baro - z_hat;   % scalar residual

    %% =========================================================
    % 3) H matrisi
    % ==========================================================
    % Error-state sırası:
    % [dp_N dp_E dp_D dv_N dv_E dv_D dthx dthy dthz dbgx dbgy dbgz dbax dbay dbaz]
    %
    % Baro yalnızca delta p_D bloğunu doğrudan görür.
    H = zeros(1,15);
    H(3) = 1;

    %% =========================================================
    % 4) Innovation covariance
    % ==========================================================
    S = H * P * H.' + params.R_baro;   % scalar

    %% =========================================================
    % 5) Kalman gain
    % ==========================================================
    K = P * H.' / S;                   % 15x1

    %% =========================================================
    % 6) Error-state düzeltmesi
    % ==========================================================
    dx_hat = K * residual;             % 15x1

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
        P = A * P * A.' + K * params.R_baro * K.';
    else
        P = (I15 - K * H) * P;
    end

    % Sayısal simetri koruması
    P = 0.5 * (P + P.');

end