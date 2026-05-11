function [state, P, residual, S, K] = update_gnss_pos(state, P, z_gps_pos, params)
%UPDATE_GNSS_POS
% GNSS position update for 18-state ESKF.
%
% Measurement model:
%   z_gps_pos = p_n + noise
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%
% H:
%   H(:,1:3) = I3
%
% Baro offset ve wind state bu ölçümde direkt gözlenmez.

    %% Boyut kontrolleri
    if ~isequal(size(P), [18 18])
        error('P 18x18 olmalıdır.');
    end

    if ~isvector(z_gps_pos) || numel(z_gps_pos) ~= 3
        error('z_gps_pos 3 elemanlı olmalıdır.');
    end

    z_gps_pos = z_gps_pos(:);

    if any(~isfinite(z_gps_pos))
        error('z_gps_pos sonlu değerlerden oluşmalıdır.');
    end

    %% Measurement prediction
    z_hat = state.p_n(:);

    residual = z_gps_pos - z_hat;

    %% Measurement Jacobian
    H = zeros(3,18);
    H(:,1:3) = eye(3);

    %% Measurement covariance
    if isfield(params, 'R_gps_pos')
        R = params.R_gps_pos;
    elseif isfield(params, 'sigma_gps_pos')
        R = diag(params.sigma_gps_pos(:).^2);
    else
        R = diag([2.0; 2.0; 3.0].^2);
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