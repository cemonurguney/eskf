function obs = analyze_observability(sim, params)
%ANALYZE_OBSERVABILITY
% 15-state lineerleştirilmiş sistem için yerel ve pencere bazlı
% gözlemlenebilirlik analizi yapar.

    t = sim.t;
    dt = sim.dt;
    N = numel(t);

    if isfield(params.obs, 'sample_step')
        sample_step = max(1, round(params.obs.sample_step));
    else
        sample_step = 1;
    end

    if isfield(params.obs, 'window_sec')
        window_len = max(1, round(params.obs.window_sec / dt));
    else
        window_len = max(1, round(5.0 / dt));
    end

    if isfield(params.obs, 'rank_tol')
        rank_tol = params.obs.rank_tol;
    else
        rank_tol = 1e-8;
    end

    idx_eval = 1:sample_step:N;
    M = numel(idx_eval);

    local_rank = nan(1,M);
    local_sigma_min = nan(1,M);
    window_rank = nan(1,M);
    window_sigma_min = nan(1,M);

    for ii = 1:M
        k = idx_eval(ii);

        state_k.p_n = sim.p_true(:,k);
        state_k.v_n = sim.v_true(:,k);
        state_k.q_nb = sim.q_true(:,k);
        state_k.b_g = sim.bg_true(:,k);
        state_k.b_a = sim.ba_true(:,k);

        imu_k.gyro_m = sim.imu_gyro(:,k);
        imu_k.accel_m = sim.imu_accel(:,k);

        [F, ~, ~] = compute_F_G_Qd(state_k, imu_k, params, dt);
        A = eye(15) + F * dt;

        Hk = measurement_matrix_local(sim, k);
        Oloc = build_obsv_local(A, Hk, 15);
        if isempty(Oloc)
            local_rank(ii) = 0;
            local_sigma_min(ii) = 0;
        else
            s = svd(Oloc);
            local_rank(ii) = sum(s > rank_tol);
            local_sigma_min(ii) = s(end);
        end

        Owin = [];
        Phi = eye(15);
        k_end = min(N, k + window_len - 1);
        for j = k:k_end
            state_j.p_n = sim.p_true(:,j);
            state_j.v_n = sim.v_true(:,j);
            state_j.q_nb = sim.q_true(:,j);
            state_j.b_g = sim.bg_true(:,j);
            state_j.b_a = sim.ba_true(:,j);

            imu_j.gyro_m = sim.imu_gyro(:,j);
            imu_j.accel_m = sim.imu_accel(:,j);
            [Fj, ~, ~] = compute_F_G_Qd(state_j, imu_j, params, dt);
            Aj = eye(15) + Fj * dt;

            Hj = measurement_matrix_local(sim, j);
            if ~isempty(Hj)
                Owin = [Owin; Hj * Phi]; %#ok<AGROW>
            end
            Phi = Aj * Phi;
        end

        if isempty(Owin)
            window_rank(ii) = 0;
            window_sigma_min(ii) = 0;
        else
            sw = svd(Owin);
            window_rank(ii) = sum(sw > rank_tol);
            window_sigma_min(ii) = sw(end);
        end
    end

    obs.t_eval = t(idx_eval);
    obs.local_rank = local_rank;
    obs.local_sigma_min = local_sigma_min;
    obs.window_rank = window_rank;
    obs.window_sigma_min = window_sigma_min;

    figure('Name','Observability Analysis');
    subplot(2,1,1);
    plot(obs.t_eval, obs.local_rank, 'LineWidth', 1.2); hold on;
    plot(obs.t_eval, obs.window_rank, 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Rank');
    legend('Local rank','Windowed rank','Location','best');
    title('Observability rank evolution');

    subplot(2,1,2);
    semilogy(obs.t_eval, max(obs.local_sigma_min, eps), 'LineWidth', 1.2); hold on;
    semilogy(obs.t_eval, max(obs.window_sigma_min, eps), 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Smallest singular value');
    legend('Local sigma min','Windowed sigma min','Location','best');
    title('Observability conditioning');
end

function H = measurement_matrix_local(sim, k)
    H = [];

    if sim.gps_pos_available(k)
        Hpos = zeros(3,15);
        Hpos(:,1:3) = eye(3);
        H = [H; Hpos]; %#ok<AGROW>
    end

    if sim.gps_vel_available(k)
        Hvel = zeros(3,15);
        Hvel(:,4:6) = eye(3);
        H = [H; Hvel]; %#ok<AGROW>
    end

    if sim.baro_available(k)
        Hb = zeros(1,15);
        Hb(1,3) = 1;
        H = [H; Hb]; %#ok<AGROW>
    end
end

function O = build_obsv_local(A, H, n)
    if isempty(H)
        O = [];
        return;
    end

    O = H;
    Ak = eye(n);
    for i = 2:n
        Ak = A * Ak;
        O = [O; H * Ak]; %#ok<AGROW>
    end
end
