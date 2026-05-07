function [state, P, residual, S, K] = update_baro(state, P, z_baro, params)
%UPDATE_BARO
% Barometer update for 16-state ESKF with estimated barometer offset.
%
% Measurement model:
%   z_baro = p_D + b_baro + noise
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro]
%
% H:
%   H(3)  = 1   -> Down position
%   H(16) = 1   -> barometer offset
%
% Not:
%   Bu model NED fixed-wing için p_n(3)=Down varsayımıyla kullanılır.
%   ENU modunda baro modelini ayrı düşünmek gerekir, yoksa eksenler yine
%   insanlığın üstüne kapanır.

    %% Boyut kontrolleri
    if ~isequal(size(P), [16 16])
        error('P 16x16 olmalıdır.');
    end

    if ~isscalar(z_baro) || ~isfinite(z_baro)
        error('z_baro sonlu scalar olmalıdır.');
    end

    if ~isfield(state, 'b_baro') || isempty(state.b_baro)
        state.b_baro = 0;
    end

    %% Measurement prediction
    z_hat = state.p_n(3) + state.b_baro;

    residual = z_baro - z_hat;

    %% Measurement Jacobian
    H = zeros(1,16);
    H(3)  = 1;
    H(16) = 1;

    %% Measurement covariance
    if isfield(params, 'R_baro')
        R = params.R_baro;
    elseif isfield(params, 'sigma_baro')
        R = params.sigma_baro^2;
    else
        R = 1.8^2;
    end

    %% Kalman update
    S = H * P * H.' + R;
    K = P * H.' / S;

    dx_hat = K * residual;

    state = inject_error_state(state, dx_hat);

    %% Covariance update
    I = eye(16);

    if isfield(params, 'use_joseph_form') && params.use_joseph_form
        P = (I - K*H) * P * (I - K*H).' + K * R * K.';
    else
        P = (I - K*H) * P;
    end

    P = 0.5 * (P + P.');
end