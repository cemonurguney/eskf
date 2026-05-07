function [state, P, params] = init_filter()
    
%INIT_FILTER
% 15-state error-state EKF için başlangıç kurulumu.

    %% =========================================================
    % 1) NOMINAL STATE BAŞLANGICI
    % ==========================================================
    state.p_n = [0; 0; 0];
    state.v_n = [0; 0; 0];
    state.q_nb = [1; 0; 0; 0];
    state.b_g = [0; 0; 0];
    state.b_a = [0; 0; 0];
    state.b_baro = 0;   % [m] barometer Down-axis measurement offset

    %% =========================================================
    % 2) FİZİKSEL SABİTLER
    % ==========================================================
    params.g = 9.80665;
    params.g_n = [0; 0; params.g];

    %% =========================================================
    % 3) IMU GÜRÜLTÜ KABULLERİ
    % ==========================================================
    params.sigma_g = deg2rad(0.08);   % tactical-grade gyro white noise std [rad/s]
    params.sigma_a = 0.020;            % tactical-grade accel white noise std [m/s^2]

    params.sigma_bg_rw = deg2rad(0.003); % gyro bias RW std [rad/s/sqrt(s)]
    params.sigma_ba_rw = 0.0015;         % accel bias RW std [m/s^2/sqrt(s)]
    %% =========================================================
    % 3) BARO BİAS GÜRÜLTÜ KABULLERİ
    % ==========================================================
    params.sigma_baro_bias0 = 12.0;       % [m] initial uncertainty of baro offset
    params.sigma_baro_bias_rw = 0.01;     % [m/sqrt(s)] slow random walk

    %% =========================================================
    % 4) ÖLÇÜM GÜRÜLTÜ KABULLERİ
    % ==========================================================
    params.sigma_gps_pos = [0.2; 0.2; 0.4];
    params.sigma_gps_vel = [0.08; 0.08; 0.12];
    params.sigma_baro    = 0.6;

    params.R_gps_pos = diag(params.sigma_gps_pos.^2);
    params.R_gps_vel = diag(params.sigma_gps_vel.^2);
    params.R_baro    = params.sigma_baro^2;

    %% =========================================================
    % 5) BAŞLANGIÇ COVARIANCE P0
    % ==========================================================
    sigma_p0  = [1.0; 1.0; 1.5];
    sigma_v0  = [0.08; 0.08; 0.10];
    sigma_th0 = deg2rad([0.4; 0.4; 2.0]);
    sigma_bg0 = deg2rad([0.03; 0.03; 0.03]);
    sigma_ba0 = [0.020; 0.020; 0.030];
    
    P = diag([
        sigma_p0.^2;
        sigma_v0.^2;
        sigma_th0.^2;
        sigma_bg0.^2;
        sigma_ba0.^2
    ]);

    %% =========================================================
    % 6) DİĞER AYARLAR
    % ==========================================================
    params.use_joseph_form = true;

    %% =========================================================
    % 7) ZAMAN / ANİMASYON AYARLARI
    % ==========================================================
    params.sim.dt = 0.005;
    params.sim.T_end = 100;

    params.anim.playback_rate = 1.0;   % 1.0 = gerçek zamanlı
    params.anim.step = 5;
    params.anim.trail_len = 120;
    params.anim.axis_len_truth = 14;
    params.anim.axis_len_est   = 12;

    %% =========================================================
    % 7.1) OBSERVABILITY ANALİZ AYARLARI
    % ==========================================================
    params.obs.enable = true;
    params.obs.window_sec = 5.0;
    params.obs.rank_tol = 1e-8;
    params.obs.sample_step = 10;

    %% =========================================================
    % 8) TRUTH / SENARYO AYARLARI
    % ==========================================================
    params.truth.bg0 = deg2rad([0.03; -0.02; 0.025]);
    params.truth.ba0 = [0.015; -0.010; 0.020];

    % -------------------------
    % Attitude senaryosu
    % -------------------------
    params.truth.attitude_mode = 'sinusoid';
    % seçenekler:
    % 'fixed'
    % 'sinusoid'
    % 'piecewise'

    params.truth.attitude.fixed_euler_deg = [0; 0; 0];

    params.truth.attitude.sin.roll_amp_deg       = 8;
    params.truth.attitude.sin.pitch_amp_deg      = 5;
    params.truth.attitude.sin.yaw_bias_rate_deg  = 0.6;  % deg/s
    params.truth.attitude.sin.yaw_amp_deg        = 10;

    params.truth.attitude.sin.roll_freq  = 0.18;
    params.truth.attitude.sin.pitch_freq = 0.11;
    params.truth.attitude.sin.yaw_freq   = 0.05;

    % interval-based örnek:
    % 0-20 s   : [0; 0; 0]
    % 20-50 s  : [5; 2; 15]
    % 50-100 s : [-3; 4; 35]
    params.truth.attitude.piecewise.t_breaks = [0 20 50 100];
    params.truth.attitude.piecewise.euler_deg = [ ...
         0   5   -3;
         0   2    4;
         0  15   35];
    % boyut: 3 x (numel(t_breaks)-1)

    % -------------------------
    % Translasyon senaryosu
    % -------------------------
    params.truth.translation_mode = 'random_smooth';
    % seçenekler:
    % 'constant_velocity'
    % 'constant_accel'
    % 'random_smooth'
    % 'piecewise_accel'

    params.truth.translation.v0_body = [5; 2; 1.5];

    params.truth.translation.const_vel_body = [5; 2; 1.5];
    params.truth.translation.const_acc_body = [0.03; -0.01; 0.00];

    params.truth.translation.random.nom_v_body = [5.2; 2.4; 1.7];
    params.truth.translation.random.tau_a = [1.8; 2.3; 2.0];
    params.truth.translation.random.sigma_a_rand = [0.55; 0.35; 0.22];
    params.truth.translation.random.tau_v = [4.0; 4.5; 4.0];
    params.truth.translation.random.seed = 7;

    params.truth.translation.piecewise.t_breaks = [0 20 40 70 100];
    params.truth.translation.piecewise.acc_body = [ ...
        0.00   0.08  -0.04   0.00;
        0.00   0.00   0.03  -0.02;
        0.00   0.00   0.00   0.01];
    % boyut: 3 x (numel(t_breaks)-1)
end