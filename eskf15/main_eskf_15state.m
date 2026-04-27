%function main_eskf_15state()
%MAIN_ESKF_15STATE
% 15-state error-state EKF ana calistirma dosyasi.
%
% VERI KAYNAGI SECIMI:
%   use_real_data = false  -> simulate_truth_and_sensors (sentetik)
%   use_real_data = true   -> build_sim_from_mat("hover_clean.mat") (gercek)
%
% REAL-TIME / ROKET NOTU:
% Bu kod offline calisir ama icindeki "pre-launch initialization" mantigi
% gercek bir roket/drone uygulamasi icin de dogrudur:
%   - Lift-off oncesi filtre stasyoner kabul edilip state'i GPS+truth'tan
%     aliyor, biaslar da gravity+durgunluk ile tahmin edilebilir durumda.
%   - Lift-off sonrasi (ivme veya harici bayrak) filtre full-mode, bu
%     pre-launch mantigi bir daha tetiklenmez.
% Real-time deployment'ta yapilacak degisiklikler en altta listelendi.

    clear; clc; close all;

    %% 0) Veri kaynagi anahtari
    use_real_data  = true;
    real_data_file = "rectangle_clean.mat";

    %% 1) Filtreyi baslat
    [state, P, params] = init_filter();

    %% 2) Simulasyon verisini uret / gercek veriyi yukle
    if use_real_data
        fprintf('[main] Gercek veri modu: %s yukleniyor...\n', real_data_file);
        sim = build_sim_from_mat(real_data_file);

        % ENU convention: data up-positive (z up).
        % Kodun default g_n NED; gercek veri modunda ENU olarak override.
        params.g_n = [0; 0; -9.81];

        % -----------------------------------------------------------------
        % PRE-LAUNCH INITIALIZATION DETECTION
        %
        % Amac: filtreyi guvenilir bir nokta ile baslatmak.
        % Kosul: drone/roket stasyoner (kucuk hiz, |a| ~ g, kucuk gyro).
        % Bu kosul saglananlar arasinda ilk "sustained" bolumun son noktasini
        % baslangic olarak sec. "Sustained" = win boyunca ustuste saglanan.
        %
        % ROKET NOTU: Bu detection kodu sadece lift-off oncesi bir kez
        % kullanilir. Lift-off sonrasi (ivme esigini astiktan itibaren)
        % filtre full-mode, bu blok bir daha calisrilmaz.
        % -----------------------------------------------------------------
        v_norm_true = vecnorm(sim.v_true,   2, 1);
        a_norm_imu  = vecnorm(sim.imu_accel,2, 1);
        w_norm_imu  = vecnorm(sim.imu_gyro, 2, 1);

        % Bu esikler esnek; her sensor setine gore tune edilebilir.
        v_thresh = 5.0;     % m/s  (glitch ~8 m/s'den kucuk olmali)
        g_band   = 3.0;     % m/s^2
        w_thresh = 0.5;     % rad/s
        win      = 50;      % 50 ornek @100Hz = 0.5 s sustained

        ok = (v_norm_true < v_thresh) ...
           & (abs(a_norm_imu - 9.81) < g_band) ...
           & (w_norm_imu < w_thresh);

        ok_sustained = movmin(double(ok), win) == 1;
        k_candidates = find(ok_sustained);

        if isempty(k_candidates)
            warning('[main] Pre-launch stasyoner bolge bulunamadi, k=1 kullaniliyor');
            k_start = 1;
        else
            % Ilk sustained bolumun orta noktasi — baslangic transient'inin
            % sonu ile takeoff'tan oncesini kapsar
            k_start = k_candidates(1);
            fprintf('[main] Pre-launch bolge: k=%d..%d arasinda uygun noktalar\n', ...
                k_candidates(1), k_candidates(end));
        end
        fprintf('[main] Baslangic indeksi: k_start = %d  (t_offset = %.2f s, %d ornek atlandi)\n', ...
            k_start, sim.t(k_start)-sim.t(1), k_start-1);

        % sim struct kirp
        t_offset              = sim.t(k_start);
        sim.t                 = sim.t(k_start:end) - t_offset;
        sim.imu_gyro          = sim.imu_gyro(:,k_start:end);
        sim.imu_accel         = sim.imu_accel(:,k_start:end);
        sim.p_true            = sim.p_true(:,k_start:end);
        sim.v_true            = sim.v_true(:,k_start:end);
        sim.q_true            = sim.q_true(:,k_start:end);
        sim.gps_pos           = sim.gps_pos(:,k_start:end);
        sim.gps_vel           = sim.gps_vel(:,k_start:end);
        sim.gps_pos_available = sim.gps_pos_available(k_start:end);
        sim.gps_vel_available = sim.gps_vel_available(k_start:end);
        sim.baro              = sim.baro(k_start:end);
        sim.baro_available    = sim.baro_available(k_start:end);
        sim.baro_gecerli =false;
        if isfield(sim,'range')
            sim.range           = sim.range(k_start:end);
            sim.range_available = sim.range_available(k_start:end);
        end

        % Guvenli baslangic state'i
        state.p_n  = sim.p_true(:,1);
        state.v_n  = sim.v_true(:,1);
        state.q_nb = sim.q_true(:,1) / norm(sim.q_true(:,1));

    else
        fprintf('[main] Sentetik simulasyon modu.\n');
        sim = simulate_truth_and_sensors(params);
    end

    t  = sim.t;
    dt = sim.dt;
    N  = numel(t);

    %% 3) Log alanlari
    log_p   = zeros(3, N);
    log_v   = zeros(3, N);
    log_q   = zeros(4, N);
    log_bg  = zeros(3, N);
    log_ba  = zeros(3, N);

    log_Pdiag = zeros(15, N);

    log_res_gps_pos = nan(3, N);
    log_res_gps_vel = nan(3, N);
    log_res_baro    = nan(1, N);

    %% 4) Ana dongu
    for k = 1:N

        imu.gyro_m  = sim.imu_gyro(:,k);
        imu.accel_m = sim.imu_accel(:,k);

        % Nominal propagation
        state = propagate_nominal(state, imu, params, dt);

        % Covariance propagation
        [F, G, Qd] = compute_F_G_Qd(state, imu, params, dt);
        P = propagate_covariance(P, F, G, Qd, dt);

        % GNSS position update
        if sim.gps_pos_available(k)
            z_gps_pos = sim.gps_pos(:,k);
            [state, P, residual, ~, ~] = update_gnss_pos(state, P, z_gps_pos, params);
            log_res_gps_pos(:,k) = residual;
        end

        % GNSS velocity update
        if sim.gps_vel_available(k)
            z_gps_vel = sim.gps_vel(:,k);
            [state, P, residual, ~, ~] = update_gnss_vel(state, P, z_gps_vel, params);
            log_res_gps_vel(:,k) = residual;
        end

        % Baro update
        if sim.baro_available(k) & sim.baro_gecerli
            z_baro = sim.baro(k);
            [state, P, residual, ~, ~] = update_baro(state, P, z_baro, params);
            log_res_baro(k) = residual;
        end

        % Log
        log_p(:,k)  = state.p_n;
        log_v(:,k)  = state.v_n;
        log_q(:,k)  = state.q_nb;
        log_bg(:,k) = state.b_g;
        log_ba(:,k) = state.b_a;

        log_Pdiag(:,k) = diag(P);
    end

    %% 5) Sonuclar
    plot_results(sim, t, log_p, log_v, log_q, log_bg, log_ba, log_Pdiag, ...
                 log_res_gps_pos, log_res_gps_vel, log_res_baro, params);

    %% 6) Observability (sadece sentetik)
    if ~use_real_data
        obs = analyze_observability(sim, params);
        disp('--- Observability summary ---');
        fprintf('Mean local rank        : %.2f\n', mean(obs.local_rank, 'omitnan'));
        fprintf('Mean windowed rank     : %.2f\n', mean(obs.window_rank, 'omitnan'));
        fprintf('Min local sigma_min    : %.3e\n', min(obs.local_sigma_min));
        fprintf('Min windowed sigma_min : %.3e\n', min(obs.window_sigma_min));
    else
        disp('[main] Gercek veri modunda observability analizi atlandi.');
    end

    %% 7) Ozet hata metrikleri
    err_p = log_p - sim.p_true;
    err_v = log_v - sim.v_true;
    fprintf('\n--- RUN FINISHED ---\n');
    fprintf('Mean pos error norm  [m]   : %.4f\n', mean(vecnorm(err_p,2,1)));
    fprintf('Final pos error norm [m]   : %.4f\n', norm(err_p(:,end)));
    fprintf('Mean vel error norm  [m/s] : %.4f\n', mean(vecnorm(err_v,2,1)));
    fprintf('Final vel error norm [m/s] : %.4f\n', norm(err_v(:,end)));
    %end

% =====================================================================
% ROKET REAL-TIME DEPLOYMENT NOTU
% =====================================================================
% Bu dosyadaki pre-launch init mantigi rokete aynen tasinir fakat
% calisma bicimi sudur:
%
% 1) Roket yerde, filtre ayakta, ama henuz "armed" degil.
%    Her IMU ornegi icin:
%       - propagate (bias tahminine yarar)
%       - GPS var ise update (pozisyon sabitlenir)
%       - "stasyoner mi" kosulu (|a|~g, gyro kucuk, GPS fix)
%         100 ornek ustuste saglandiginda "initialization OK" bayragini set et
%
% 2) Initialization OK olana kadar arm'i iptal et / go-no-go bayragi.
%
% 3) Lift-off: |a| > 2g esigini ilk gectigi an, veya external arming sinyali.
%    Pre-launch blok bir daha calisrilmaz. Filtre full-mode.
%
% 4) Ucus boyunca normal ESKF.
%
% Bu offline script ile real-time arasindaki tek yapisal fark: offline'da
% "k_start'dan sonrasi uc" mantigi var, real-time'da "init-ok'tan sonrasi
% yasal ucus" mantigi var. Ikisi de ayni kritere bakar.
% =====================================================================