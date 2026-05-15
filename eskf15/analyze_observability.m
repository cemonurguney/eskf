function obs = analyze_observability(sim, params)
%ANALYZE_OBSERVABILITY
% 18-state lineerleştirilmiş sistem için yerel ve pencere bazlı
% gözlemlenebilirlik analizi yapar.
%
% Error-state:
%   dx = [dp; dv; dtheta; dbg; dba; db_baro; dw_N; dw_E]
%      = 18x1
%
% Ölçüm modelleri:
%   GPS position : p
%   GPS velocity : v
%   Baro         : p_D + b_baro
%   TAS          : ||v - wind||
%
% Not:
% - Synthetic için sim.p_true / sim.v_true / sim.q_true gerçek truth olabilir.
% - Real veri için bu alanlara estimated logs konabilir:
%       p_true = log_p
%       v_true = log_v
%       q_true = log_q
%   Yani analiz estimated nominal trajectory etrafında yapılır.
%
% Plot:
% - Local rank/sigma hesaplanır ama varsayılan olarak çizdirilmez.
% - Sunum için windowed rank + windowed sigma_min daha anlamlıdır.
% - GPS OFF bölgeleri gri gölgelenebilir.

    %% ---------------- Basic setup ----------------
    t = sim.t(:).';

    if isfield(sim, 'dt')
        dt = sim.dt;
    else
        dt = median(diff(t), 'omitnan');
    end

    N = numel(t);
    n_err = 18;

    if isfield(params, 'obs') && isfield(params.obs, 'sample_step')
        sample_step = max(1, round(params.obs.sample_step));
    else
        sample_step = 1000;
    end

    if isfield(params, 'obs') && isfield(params.obs, 'window_sec')
        window_len = max(1, round(params.obs.window_sec / dt));
    else
        window_len = max(1, round(3.0 / dt));
    end

    if isfield(params, 'obs') && isfield(params.obs, 'rank_tol')
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

    %% ---------------- Observability loop ----------------
    for ii = 1:M
        k = idx_eval(ii);

        state_k = state_from_sim_local(sim, k);

        imu_k.gyro_m = sim.imu_gyro(:,k);
        imu_k.accel_m = sim.imu_accel(:,k);

        [F, ~, ~] = compute_F_G_Qd(state_k, imu_k, params, dt);

        A = eye(n_err) + F * dt;
        Hk = measurement_matrix_local(sim, k, n_err);

        %% Local observability
        Oloc = build_obsv_local(A, Hk, n_err);

        if isempty(Oloc)
            local_rank(ii) = 0;
            local_sigma_min(ii) = 0;
        else
            s = svd(Oloc);
            local_rank(ii) = sum(s > rank_tol);
            local_sigma_min(ii) = s(end);
        end

        %% Windowed observability
        Owin = [];
        Phi = eye(n_err);

        k_end = min(N, k + window_len - 1);

        for j = k:k_end
            state_j = state_from_sim_local(sim, j);

            imu_j.gyro_m = sim.imu_gyro(:,j);
            imu_j.accel_m = sim.imu_accel(:,j);

            [Fj, ~, ~] = compute_F_G_Qd(state_j, imu_j, params, dt);

            Aj = eye(n_err) + Fj * dt;
            Hj = measurement_matrix_local(sim, j, n_err);

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

    %% ---------------- Output struct ----------------
    obs.t_eval = t(idx_eval);

    obs.local_rank = local_rank;
    obs.local_sigma_min = local_sigma_min;

    obs.window_rank = window_rank;
    obs.window_sigma_min = window_sigma_min;

    obs.n_err = n_err;
    obs.sample_step = sample_step;
    obs.window_sec = window_len * dt;
    obs.rank_tol = rank_tol;

    %% ---------------- Plot options ----------------
    show_local = false;
    show_gps_schedule = true;

    if isfield(params, 'obs') && isfield(params.obs, 'show_local')
        show_local = params.obs.show_local;
    end

    if isfield(params, 'obs') && isfield(params.obs, 'show_gps_schedule')
        show_gps_schedule = params.obs.show_gps_schedule;
    end

    %% ---------------- Plot ----------------
    figure('Name','Observability Analysis');

    %% Rank plot
    subplot(2,1,1);
    hold on;
    grid on;

    ylim_rank = [0, n_err + 1];
    ylim(ylim_rank);

    if show_gps_schedule && isfield(sim, 'gps_meas_scheduled')
        shade_gps_off_regions_local(t, sim.gps_meas_scheduled);
    end

    if show_local
        plot(obs.t_eval, obs.local_rank, ...
            'Color', [0.55 0.55 0.55], ...
            'LineWidth', 0.9, ...
            'DisplayName', 'Local rank');
    end

    plot(obs.t_eval, obs.window_rank, ...
        'LineWidth', 1.7, ...
        'DisplayName', 'Windowed rank');

    yline(n_err, '--k', ...
        'LineWidth', 1.0, ...
        'DisplayName', 'n_{err}');

    xlim([0 t(end)]);
    ylim(ylim_rank);

    xlabel('Time [s]');
    ylabel('Rank');
    legend('Location','best');
    title('Windowed Observability Rank');

    %% Conditioning plot
    subplot(2,1,2);
    hold on;
    grid on;

    if show_gps_schedule && isfield(sim, 'gps_meas_scheduled')
        shade_gps_off_regions_local(t, sim.gps_meas_scheduled);
    end

    if show_local
        semilogy(obs.t_eval, max(obs.local_sigma_min, eps), ...
            'Color', [0.55 0.55 0.55], ...
            'LineWidth', 0.9, ...
            'DisplayName', 'Local \sigma_{min}');
    end

    semilogy(obs.t_eval, max(obs.window_sigma_min, eps), ...
        'LineWidth', 1.7, ...
        'DisplayName', 'Windowed \sigma_{min}');

    xlim([0 t(end)]);

    xlabel('Time [s]');
    ylabel('Smallest singular value');
    legend('Location','best');
    title('Windowed Observability Conditioning');
end

%% ============================================================
% Local helper: estimated/truth state from sim struct
% ============================================================
function state = state_from_sim_local(sim, k)
    state.p_n = sim.p_true(:,k);
    state.v_n = sim.v_true(:,k);
    state.q_nb = sim.q_true(:,k);

    if isfield(sim, 'bg_true')
        state.b_g = sim.bg_true(:,k);
    else
        state.b_g = [0;0;0];
    end

    if isfield(sim, 'ba_true')
        state.b_a = sim.ba_true(:,k);
    else
        state.b_a = [0;0;0];
    end

    if isfield(sim, 'b_baro_true')
        state.b_baro = sim.b_baro_true(k);
    else
        state.b_baro = 0;
    end

    if isfield(sim, 'wind_ref_available') && sim.wind_ref_available(k)
        state.wind_ne = sim.wind_ref(:,k);
    else
        state.wind_ne = [0;0];
    end
end

%% ============================================================
% Local helper: measurement matrix at time index k
% ============================================================
function H = measurement_matrix_local(sim, k, n_err)
    H = [];

    %% GPS schedule
    gps_scheduled = true;

    if isfield(sim, 'gps_meas_scheduled')
        gps_scheduled = sim.gps_meas_scheduled(k);
    end

    %% GPS position
    if gps_scheduled && sim.gps_pos_available(k)
        Hpos = zeros(3,n_err);
        Hpos(:,1:3) = eye(3);
        H = [H; Hpos]; %#ok<AGROW>
    end

    %% GPS velocity
    if gps_scheduled && sim.gps_vel_available(k)
        Hvel = zeros(3,n_err);
        Hvel(:,4:6) = eye(3);
        H = [H; Hvel]; %#ok<AGROW>
    end

    %% Barometer
    if sim.baro_available(k)
        Hb = zeros(1,n_err);

        % z_baro = p_D + b_baro
        Hb(1,3) = 1;

        % Bizim implementasyonda b_baro sadece GPS height recent iken
        % estimate ediliyor. Eğer flag yoksa theoretical/default olarak
        % H(16)=1 kabul edilir.
        estimate_baro_bias_here = true;

        if isfield(sim, 'estimate_baro_bias_scheduled')
            estimate_baro_bias_here = sim.estimate_baro_bias_scheduled(k);
        end

        if estimate_baro_bias_here
            Hb(1,16) = 1;
        else
            Hb(1,16) = 0;
        end

        H = [H; Hb]; %#ok<AGROW>
    end

    %% TAS / airspeed
    if isfield(sim, 'airspeed_available') && sim.airspeed_available(k)

        if isfield(sim, 'wind_ref_available') && sim.wind_ref_available(k)

            v_n = sim.v_true(:,k);
            w_ne = sim.wind_ref(:,k);

            w_n = [w_ne(1); w_ne(2); 0];
            v_air = v_n - w_n;
            tas = norm(v_air);

            if tas > 1e-6
                e_air = v_air / tas;

                Htas = zeros(1,n_err);

                % h = ||v - w||
                Htas(1,4:6) = e_air.';

                % wind_N, wind_E sensitivity
                Htas(1,17) = -e_air(1);
                Htas(1,18) = -e_air(2);

                H = [H; Htas]; %#ok<AGROW>
            end
        end
    end
end

%% ============================================================
% Local helper: observability matrix construction
% ============================================================
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

%% ============================================================
% Local helper: shade GPS OFF regions
% ============================================================
function shade_gps_off_regions_local(t, gps_flag)
%SHADE_GPS_OFF_REGIONS_LOCAL
% GPS OFF bölgelerini mevcut eksende gri arka plan olarak gösterir.

    gps_flag = logical(gps_flag(:).');
    t = t(:).';

    if numel(gps_flag) ~= numel(t)
        return;
    end

    off_flag = ~gps_flag;

    edges = diff([false, off_flag, false]);
    starts = find(edges == 1);
    stops  = find(edges == -1) - 1;

    yl = ylim;

    for i = 1:numel(starts)
        xs = t(starts(i));
        xe = t(stops(i));

        patch([xs xe xe xs], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.90 0.90 0.90], ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.45, ...
            'HandleVisibility', 'off');
    end

    % Patch arkada kalsın, çizgiler üstte görünsün.
    h_lines = findobj(gca, 'Type', 'line');

    if ~isempty(h_lines)
        uistack(h_lines, 'top');
    end
end