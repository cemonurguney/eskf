function P = propagate_covariance(P, F, ~, Qd, dt)
%PROPAGATE_COVARIANCE
% 15-state error-state EKF için covariance propagation.

    if ~isequal(size(P), [15 15])
        error('P matrisi 15x15 olmalıdır.');
    end

    if ~isequal(size(F), [15 15])
        error('F matrisi 15x15 olmalıdır.');
    end

    if ~isequal(size(Qd), [15 15])
        error('Qd matrisi 15x15 olmalıdır.');
    end

    I15 = eye(15);
    Phi = I15 + F*dt + 0.5*(F*F)*dt^2;

    P = Phi * P * Phi.' + Qd;
    P = 0.5 * (P + P.');
end