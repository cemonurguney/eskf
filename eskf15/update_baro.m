function [state, P, residual, S, K] = update_baro(state, P, z_baro, params)
%UPDATE_BARO
% Barometer update for 18-state ESKF with estimated barometer offset.
%
% Measurement model:
%   z_baro = p_D + b_baro + noise
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%
% H:
%   H(3)  = 1   -> Down position
%   H(16) = 1   -> barometer offset
%
% Wind state bu ölçümde direkt gözlenmez.

    %% Boyut kontrolleri
    if ~isequal(size(P), [18 18])
        error('P 18x18 olmalıdır.');
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
    H = zeros(1,18);
    H(3)  = 1;
    H(16) = 1;

    %% Measurement covariance
    if isfield(params, 'R_baro')
        R = params.R_baro;
    elseif isfield(params, 'sigma_baro')
        R = params.sigma_baro^2;
    else
        R = 2.5^2;
    end

    %% Kalman update
    S = H * P * H.' + R;
    K = P * H.' / S;

    dx_hat = K * residual;

    state = inject_error_state(state, dx_hat);

    %% Covariance update
    I = eye(18);

    if isfield(params, 'use_joseph_form') && params.use_joseph_form
        P = (I - K*H) * P * (I - K*H).' + K * R * K.';
    else
        P = (I - K*H) * P;
    end

    P = 0.5 * (P + P.');
end