function [state, P, residual, S, K, accepted] = update_gnss_pos(state, P, z_gps_pos, params)
%UPDATE_GNSS_POS
% 18-state ESKF GNSS position update with innovation gate.
%
% Measurement:
%   z = p_n + noise
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]

    accepted = true;

    if ~isequal(size(P), [18 18])
        error('P 18x18 olmalıdır.');
    end

    z_gps_pos = z_gps_pos(:);

    if numel(z_gps_pos) ~= 3 || any(~isfinite(z_gps_pos))
        error('z_gps_pos 3x1 sonlu vektör olmalıdır.');
    end

    z_hat = state.p_n(:);
    residual = z_gps_pos - z_hat;

    H = zeros(3,18);
    H(:,1:3) = eye(3);

    if isfield(params, 'R_gps_pos')
        R = params.R_gps_pos;
    elseif isfield(params, 'sigma_gps_pos')
        R = diag(params.sigma_gps_pos(:).^2);
    else
        R = diag([3;3;5].^2);
    end

    S = H * P * H.' + R;
    S = 0.5 * (S + S.');

    K = P * H.' / S;

    d2 = residual.' / S * residual;

    if isfield(params, 'gps_pos_gate_chi2') && isfinite(params.gps_pos_gate_chi2)
        if d2 > params.gps_pos_gate_chi2
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