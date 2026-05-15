function [state, P, residual, S, K, accepted] = update_baro(state, P, z_baro, params)
%UPDATE_BARO
% Barometer update for 18-state ESKF with optional barometer offset estimation.
%
% Measurement model:
%   z_baro = p_D + b_baro + noise
%
% If params.estimate_baro_bias = true:
%   H(3)  = 1
%   H(16) = 1
%
% If params.estimate_baro_bias = false:
%   H(3)  = 1
%   H(16) = 0
%
% In the second case, b_baro is used as fixed compensation in z_hat,
% but the baro measurement does not directly estimate db_baro.

    accepted = true;

    if ~isequal(size(P), [18 18])
        error('P 18x18 olmalıdır.');
    end

    if ~isscalar(z_baro) || ~isfinite(z_baro)
        error('z_baro sonlu scalar olmalıdır.');
    end

    if ~isfield(state, 'b_baro') || isempty(state.b_baro)
        state.b_baro = 0;
    end

    z_hat = state.p_n(3) + state.b_baro;
    residual = z_baro - z_hat;

    H = zeros(1,18);
    H(3) = 1;

    if isfield(params, 'estimate_baro_bias') && params.estimate_baro_bias
        H(16) = 1;
    else
        H(16) = 0;
    end

    if isfield(params, 'R_baro')
        R = params.R_baro;
    elseif isfield(params, 'sigma_baro')
        R = params.sigma_baro^2;
    else
        R = 3.0^2;
    end

    S = H * P * H.' + R;
    K = P * H.' / S;

    d2 = residual^2 / S;

    if isfield(params, 'baro_gate_chi2') && isfinite(params.baro_gate_chi2)
        if d2 > params.baro_gate_chi2
            accepted = false;
            K = zeros(size(P,1), 1);
            return;
        end
    end

    dx_hat = K * residual;

    state = inject_error_state(state, dx_hat);

    I = eye(18);

    if isfield(params, 'use_joseph_form') && params.use_joseph_form
        P = (I - K*H) * P * (I - K*H).' + K * R * K.';
    else
        P = (I - K*H) * P;
    end

    P = 0.5 * (P + P.');
end