function state = inject_error_state(state, dx_hat)
%INJECT_ERROR_STATE
% Error-state EKF update sonrası bulunan hata düzeltmesini nominal state'e uygular.
%
% Bu sürümde kullanılan konvansiyon:
% - q_nb : body -> nav
% - küçük attitude düzeltmesi body-frame küçük dönüş gibi uygulanır
% - bu yüzden quaternion düzeltmesi sağdan uygulanır:
%       q_new = q_old ⊗ dq

    %% 1) Boyut kontrolü
    if ~isvector(dx_hat) || numel(dx_hat) ~= 15
        error('dx_hat, 15 elemanlı error-state düzeltmesi olmalıdır.');
    end

    dx_hat = dx_hat(:);

    %% 2) Error-state blokları
    d_p     = dx_hat(1:3);
    d_v     = dx_hat(4:6);
    d_theta = dx_hat(7:9);
    d_bg    = dx_hat(10:12);
    d_ba    = dx_hat(13:15);

    %% 3) Additive state düzeltmeleri
    state.p_n = state.p_n + d_p;
    state.v_n = state.v_n + d_v;
    state.b_g = state.b_g + d_bg;
    state.b_a = state.b_a + d_ba;

    %% 4) Attitude düzeltmesi (RIGHT multiplication)
    dq = small_angle_to_quat(d_theta);

    state.q_nb = quat_multiply(state.q_nb, dq);
    state.q_nb = quat_normalize(state.q_nb);
end