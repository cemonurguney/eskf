function [state, P, residual, S, K] = update_airspeed(state, P, z_tas, params)
%UPDATE_AIRSPEED
% Airspeed/TAS update for 18-state ESKF with horizontal wind states.
%
% Measurement model:
%   z_tas = || v_n - wind_n || + noise
%
% where:
%   wind_n = [wind_N; wind_E; 0]
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%
% H:
%   d TAS / d v_n     = v_air_n' / ||v_air_n||
%   d TAS / d wind_NE = -[v_air_N, v_air_E] / ||v_air_n||
%
% Not:
%   Bu update TAS magnitude kullanır. Yani yalnız başına rüzgar yönünü
%   mucize gibi çözmez; gözlenebilirlik uçuş manevrası + GPS velocity +
%   airspeed değişimiyle gelir. Matematik de sihirbaz değil sonuçta.

    %% Boyut kontrolleri
    if ~isequal(size(P), [18 18])
        error('P 18x18 olmalıdır.');
    end

    if ~isscalar(z_tas) || ~isfinite(z_tas)
        error('z_tas sonlu scalar olmalıdır.');
    end

    if ~isfield(state, 'wind_ne') || isempty(state.wind_ne)
        state.wind_ne = [0; 0];
    end

    state.wind_ne = state.wind_ne(:);

    if numel(state.wind_ne) ~= 2
        error('state.wind_ne 2x1 olmalıdır: [wind_N; wind_E].');
    end

    if ~isfield(state, 'v_n') || numel(state.v_n) ~= 3
        error('state.v_n 3x1 olmalıdır.');
    end

    %% Measurement prediction
    wind_n = [state.wind_ne; 0];

    v_air_n = state.v_n(:) - wind_n;
    tas_hat = norm(v_air_n);

    % Çok düşük airspeed'te Jacobian güvenilmez.
    if tas_hat < 0.5
        residual = nan;
        S = nan;
        K = nan(18,1);
        return;
    end

    residual = z_tas - tas_hat;

    %% Measurement Jacobian
    H = zeros(1,18);

    d_tas_d_v = v_air_n(:).' / tas_hat;

    % Velocity states: dv_N, dv_E, dv_D
    H(4:6) = d_tas_d_v;

    % Wind states: dw_N, dw_E
    H(17) = -v_air_n(1) / tas_hat;
    H(18) = -v_air_n(2) / tas_hat;

    %% Measurement covariance
    if isfield(params, 'R_tas')
        R = params.R_tas;
    elseif isfield(params, 'sigma_tas')
        R = params.sigma_tas^2;
    else
        R = 1.5^2;
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