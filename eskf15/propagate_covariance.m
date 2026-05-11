function P = propagate_covariance(P, F, ~, Qd, dt)
%PROPAGATE_COVARIANCE
% 18-state ESKF covariance propagation.
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%
% Propagation:
%   Phi ≈ I + F dt + 0.5 F^2 dt^2
%   P_k+1 = Phi P Phi' + Qd

    %% Boyut kontrolleri
    if ~isequal(size(P), [18 18])
        error('P matrisi 18x18 olmalıdır.');
    end

    if ~isequal(size(F), [18 18])
        error('F matrisi 18x18 olmalıdır.');
    end

    if ~isequal(size(Qd), [18 18])
        error('Qd matrisi 18x18 olmalıdır.');
    end

    if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error('dt pozitif ve sonlu scalar olmalıdır.');
    end

    %% State transition
    I18 = eye(18);
    Phi = I18 + F*dt + 0.5*(F*F)*dt^2;

    %% Covariance propagation
    P = Phi * P * Phi.' + Qd;

    %% Symmetrize
    P = 0.5 * (P + P.');
end