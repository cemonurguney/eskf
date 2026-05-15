function [state, P, residual, S, K, accepted] = update_airspeed(state, P, z_tas, params)
%UPDATE_AIRSPEED
% TAS-aided horizontal wind update for 18-state ESKF.
%
% Measurement model:
%   z_TAS = || v_n - w_n || + noise
%
% Wind state:
%   state.wind_ne = [w_N; w_E]
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%
% Jacobian:
%   h = norm(v_n - [w_N; w_E; 0])
%   dh/dv = unit(v_air)
%   dh/dw_N = -unit_N
%   dh/dw_E = -unit_E

    accepted = true;

    if ~isequal(size(P), [18 18])
        error('P 18x18 olmalıdır.');
    end

    if ~isscalar(z_tas) || ~isfinite(z_tas)
        error('z_tas sonlu scalar olmalıdır.');
    end

    if ~isfield(state, 'wind_ne') || isempty(state.wind_ne)
        state.wind_ne = [0;0];
    end

    v_n = state.v_n(:);
    wind_ne = state.wind_ne(:);

    wind_n = [wind_ne(1); wind_ne(2); 0];

    v_air_n = v_n - wind_n;

    tas_hat = norm(v_air_n);

    if tas_hat < 1e-6
        residual = z_tas;
        S = nan;
        K = zeros(18,1);
        accepted = false;
        return;
    end

    residual = z_tas - tas_hat;

    e_air = v_air_n / tas_hat;

    H = zeros(1,18);

    % velocity sensitivity
    H(4:6) = e_air.';

    % wind sensitivity
    H(17) = -e_air(1);
    H(18) = -e_air(2);

    if isfield(params, 'R_tas')
        R = params.R_tas;
    elseif isfield(params, 'sigma_tas')
        R = params.sigma_tas^2;
    else
        R = 2.5^2;
    end

    S = H * P * H.' + R;
    K = P * H.' / S;

    d2 = residual^2 / S;

    if isfield(params, 'tas_gate_chi2') && isfinite(params.tas_gate_chi2)
        if d2 > params.tas_gate_chi2
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