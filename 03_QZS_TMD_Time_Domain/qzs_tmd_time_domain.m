function result = qzs_tmd_time_domain()
    % Author: Yu Gao
    % Purpose: Compute increasing/decreasing harmonic-frequency sweeps of the nonlinear QZS-TMD.
    % QZS_TMD_TIME_DOMAIN
    % Source implementation: Time_domain_comparison_FrequencyRatio.m.
    % States: [x1; r; dx1; dr], with r = x2 - x1.
    % The nonlinear restoring force and time-integration/FFT algorithm are retained.
    % The source parameter case is mu = 0.05, alpha = 0.25, xi_2 = 0.01, f = 1 N.
    % Output frequencies are in Hz; amplitudes are normalised by f/k_1.
    % primary_* uses x1; absorber_* uses x1+r, as in the source calculation.
    % Input data: none. Outputs are returned in memory; no figures/files are written.
    % Requires MATLAB. A full sweep is computationally intensive.

    %% Physical parameters and linear matrices
    m_1   = 100;
    k_1   = 100;
    L_0   = 2.5;
    m_tmd = 5;
    mu    = m_tmd / m_1;
    f_opt   = 1 / (1 + mu);
    xi_2    = 0.01;
    Omega_1 = sqrt(k_1 / m_1);
    Omega_2 = Omega_1 * f_opt;
    k_tmd = Omega_2^2 * m_tmd;
    c_tmd = 2 * xi_2 * Omega_2 * m_tmd;
    c_1   = 0.01 * (2 * Omega_1 * m_1);

    M_rel = [m_1 + m_tmd,   m_tmd
             m_tmd,         m_tmd];
    C_rel = [c_1,           0
             0,             c_tmd];
    K_rel = [k_1,           0
             0,             k_tmd];

    %% Frequency sweep and integration settings
    Fs   = 20.48 / 5;
    nfft = 4096 * 4;
    tmax = 4100;
    ts   = (0:1 / Fs:tmax).';
    options = odeset('AbsTol', [1e-10 1e-10 1e-8 1e-8]);
    f = 1.0;
    alpha_list = [0.25];

    ff_up   = 0.12:0.0005:0.22;
    ff_down = fliplr(ff_up);
    X0 = f / k_1;

    %% Forward and reverse sweeps
    result = struct([]);
    for j = 1:length(alpha_list)
        alpha = alpha_list(j);
        gamma = 2 * alpha / (2 * alpha + 1);
        a = L_0 * gamma;
        h = sqrt(L_0^2 - a^2);

        r_bias = 0.0;
        Fo_exact = @(r) -2 * k_tmd * alpha * (h + (r_bias + r)) .* ...
                        (1 - sqrt(h^2 + a^2) ./ sqrt((h + (r_bias + r)).^2 + a^2));
        fprintf('Time-domain forward sweep...\n');
        [X1_td_up_norm, X2_td_up_norm] = run_td_rel_sweep_keep_method( ...
            ff_up, ts, M_rel, C_rel, K_rel, f, X0, Fo_exact, options, Fs, nfft);

        fprintf('Time-domain down sweep...\n');
        [X1_td_down_tmp, X2_td_down_tmp] = run_td_rel_sweep_keep_method( ...
            ff_down, ts, M_rel, C_rel, K_rel, f, X0, Fo_exact, options, Fs, nfft);

        X1_td_down_norm = fliplr(X1_td_down_tmp);
        X2_td_down_norm = fliplr(X2_td_down_tmp);
        result(j).alpha = alpha;
        result(j).frequency_hz = ff_up;
        result(j).primary_up = X1_td_up_norm;
        result(j).primary_down = X1_td_down_norm;
        result(j).absorber_up = X2_td_up_norm;
        result(j).absorber_down = X2_td_down_norm;
    end
end

% Continuation and FFT amplitude extraction.
function [X1_td_norm, X2_td_norm] = run_td_rel_sweep_keep_method( ...
    ff_sweep, ts, M_rel, C_rel, K_rel, f, X0, Fo_exact, options, Fs, nfft)
    IC = [0; 0; 0; 0];
    w0 = 2 * pi * ff_sweep(1);
    [~, Z0] = ode45( ...
        @(t, z) DOF2_relCoord_Stiffness(t, z, M_rel, C_rel, K_rel, f, w0, Fo_exact), ...
        ts, IC, options);
    Tw = ceil(Fs / ff_sweep(1));
    [~, idx_peak] = max(Z0(end - Tw:end, 1));
    IC = Z0(end - Tw + idx_peak - 1, :).';
    X1_td_norm = zeros(size(ff_sweep));
    X2_td_norm = zeros(size(ff_sweep));
    for i = 1:length(ff_sweep)
        w = 2 * pi * ff_sweep(i);
        [~, Z] = ode45( ...
            @(t, z) DOF2_relCoord_Stiffness(t, z, M_rel, C_rel, K_rel, f, w, Fo_exact), ...
            ts, IC, options);
        seg_x1 = Z(end - nfft + 1:end, 1);
        seg_r  = Z(end - nfft + 1:end, 2);
        seg_x2 = seg_x1 + seg_r;
        tseg   = ts(end - nfft + 1:end);
        seg_f  = f * sin(w * tseg);
        Zfft_x1 = fft(seg_x1, nfft);
        Zfft_x2 = fft(seg_x2, nfft);
        Ffft = fft(seg_f, nfft);
        [~, idx_z] = max(abs(Zfft_x1));
        Tw = ceil(Fs / ff_sweep(i));
        [~, idx_peak] = max(Z(end - Tw:end, 1));
        IC = Z(end - Tw + idx_peak - 1, :).';
        X1_td_amp = abs(Zfft_x1(idx_z) / Ffft(idx_z)) * f;
        X2_td_amp = abs(Zfft_x2(idx_z) / Ffft(idx_z)) * f;
        X1_td_norm(i) = X1_td_amp / X0;
        X2_td_norm(i) = X2_td_amp / X0;
        if mod(i, 20) == 0
            fprintf('  TD sweep: %d/%d  (%.4f Hz)\n', i, numel(ff_sweep), ff_sweep(i));
        end
    end
end

% Nonlinear time-domain state equation.
function dzdt = DOF2_relCoord_Stiffness(t, z, M_rel, C_rel, K_rel, f, w, Fo_exact)
    q  = z(1:2);
    qd = z(3:4);
    P   = [f * sin(w * t); Fo_exact(q(2))];
    qdd = M_rel \ (P - C_rel * qd - K_rel * q);
    dzdt = [qd; qdd];
end
