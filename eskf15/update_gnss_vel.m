function [state, P, residual, S, K, accepted] = update_gnss_vel(state, P, z_gps_vel, params)
%UPDATE_GNSS_VEL
% 18-state ESKF GNSS velocity update with innovation gate.
%
% Measurement:
%   z = v_n + noise

    accepted = true;

    if ~isequal(size(P), [18 18])
        error('P 18x18 olmalıdır.');
    end

    z_gps_vel = z_gps_vel(:);

    if numel(z_gps_vel) ~= 3 || any(~isfinite(z_gps_vel))
        error('z_gps_vel 3x1 sonlu vektör olmalıdır.');
    end

    z_hat = state.v_n(:);
    residual = z_gps_vel - z_hat;

    H = zeros(3,18);
    H(:,4:6) = eye(3);

    if isfield(params, 'R_gps_vel')
        R = params.R_gps_vel;
    elseif isfield(params, 'sigma_gps_vel')
        R = diag(params.sigma_gps_vel(:).^2);
    else
        R = diag([0.5;0.5;0.8].^2);
    end

    S = H * P * H.' + R;
    S = 0.5 * (S + S.');

    K = P * H.' / S;

    d2 = residual.' / S * residual;

    if isfield(params, 'gps_vel_gate_chi2') && isfinite(params.gps_vel_gate_chi2)
        if d2 > params.gps_vel_gate_chi2
            accepted = false;
            K = zeros(size(P,1), size(H,1));
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