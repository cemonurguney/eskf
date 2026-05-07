function P = propagate_covariance(P, F, ~, Qd, dt)
%PROPAGATE_COVARIANCE
% 16-state ESKF covariance propagation.
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro]
%
% Propagation:
%   Phi ≈ I + F dt + 0.5 F^2 dt^2
%   P_k+1 = Phi P Phi' + Qd

    %% Boyut kontrolleri
    if ~isequal(size(P), [16 16])
        error('P matrisi 16x16 olmalıdır.');
    end

    if ~isequal(size(F), [16 16])
        error('F matrisi 16x16 olmalıdır.');
    end

    if ~isequal(size(Qd), [16 16])
        error('Qd matrisi 16x16 olmalıdır.');
    end

    if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error('dt pozitif ve sonlu scalar olmalıdır.');
    end

    %% State transition
    I16 = eye(16);
    Phi = I16 + F*dt + 0.5*(F*F)*dt^2;

    %% Covariance propagation
    P = Phi * P * Phi.' + Qd;

    %% Symmetrize
    P = 0.5 * (P + P.');
end