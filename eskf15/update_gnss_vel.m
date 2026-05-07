function [state, P, residual, S, K] = update_gnss_vel(state, P, z_gps_vel, params)
%UPDATE_GNSS_VEL
% GNSS velocity update for 16-state ESKF.
%
% Measurement model:
%   z_gps_vel = v_n + noise
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro]
%
% H:
%   H(:,4:6) = I3
%
% Baro offset state bu ölçümde direkt gözlenmez, bu yüzden H(:,16)=0.

    %% Boyut kontrolleri
    if ~isequal(size(P), [16 16])
        error('P 16x16 olmalıdır.');
    end

    if ~isvector(z_gps_vel) || numel(z_gps_vel) ~= 3
        error('z_gps_vel 3 elemanlı olmalıdır.');
    end

    z_gps_vel = z_gps_vel(:);

    if any(~isfinite(z_gps_vel))
        error('z_gps_vel sonlu değerlerden oluşmalıdır.');
    end

    %% Measurement prediction
    z_hat = state.v_n;

    residual = z_gps_vel - z_hat;

    %% Measurement Jacobian
    H = zeros(3,16);
    H(:,4:6) = eye(3);

    %% Measurement covariance
    if isfield(params, 'R_gps_vel')
        R = params.R_gps_vel;
    elseif isfield(params, 'sigma_gps_vel')
        R = diag(params.sigma_gps_vel(:).^2);
    else
        R = diag([0.5; 0.5; 0.8].^2);
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