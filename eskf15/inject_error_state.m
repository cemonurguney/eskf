function state = inject_error_state(state, dx_hat)
%INJECT_ERROR_STATE
% Error-state EKF update sonrası bulunan hata düzeltmesini nominal state'e uygular.
%
% 16-state error-state sırası:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E;]
%
% Nominal state:
%   state.p_n     : 3x1 position, NED/ENU kullanılan main'e göre
%   state.v_n     : 3x1 velocity
%   state.q_nb    : 4x1 quaternion, body -> nav
%   state.b_g     : 3x1 gyro bias
%   state.b_a     : 3x1 accel bias
%   state.b_baro  : scalar barometer offset
%
% Quaternion convention:
% - q_nb : body -> nav
% - küçük attitude düzeltmesi body-frame küçük dönüş gibi uygulanır
% - bu yüzden quaternion düzeltmesi sağdan uygulanır:
%       q_new = q_old ⊗ dq

    %% Boyut kontrolü
    if ~isvector(dx_hat) || numel(dx_hat) ~= 18
        error('dx_hat 18 elemanlı olmalıdır: [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E;].');
    end

    dx_hat = dx_hat(:);

    %% Error-state blokları
    d_p      = dx_hat(1:3);
    d_v      = dx_hat(4:6);
    d_theta  = dx_hat(7:9);
    d_bg     = dx_hat(10:12);
    d_ba     = dx_hat(13:15);
    d_bbaro  = dx_hat(16);
    d_wind   = dx_hat(17:18);
    %% Additive nominal düzeltmeler
    state.p_n = state.p_n + d_p;
    state.v_n = state.v_n + d_v;
    state.b_g = state.b_g + d_bg;
    state.b_a = state.b_a + d_ba;

    if ~isfield(state, 'b_baro') || isempty(state.b_baro)
        state.b_baro = 0;
    end

    state.b_baro = state.b_baro + d_bbaro;
    if ~isfield(state, 'wind_ne') || isempty(state.wind_ne)
        state.wind_ne = [0;0];
    end
    state.wind_ne = state.wind_ne(:) + d_wind;
    
    %% Attitude düzeltmesi
    dq = small_angle_to_quat(d_theta);

    state.q_nb = quat_multiply(state.q_nb, dq);
    state.q_nb = quat_normalize(state.q_nb);
end