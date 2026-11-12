function Run_Active_Damping_Optimal_FullModel()
    % Author: Yu Gao
    % Purpose: Compare uncontrolled, passive-TMD and AVD-controlled platform responses in Simulink.
    % RUN_ACTIVE_DAMPING_OPTIMAL_FULLMODEL
    % Run from this folder with the supplied SLX file on the MATLAB path.
    % Source configuration: alpha = 0.5:0.1:1.5; target frequency = 0.036 Hz.
    % The alpha sweep includes the thesis reference value 0.8.
    % The controller applies acceleration, velocity and displacement feedback.
    % The SLX filename retains DA for compatibility; the driver also sets Ca1/Ca2.
    % The SLX generates its excitation; no archived response data are required.
    % Outputs: in-memory simulation results, comparison plots and RMS values.
    % Requires MATLAB, Simulink, Symbolic Math Toolbox and Signal Processing Toolbox.
    % Input data: none. The retained configuration is specific to this ROM study.

    clc;
    close all;

    fprintf('========================================================\n');
    fprintf(' ATMD Analysis: Active vs Passive vs No-TMD \n');
    fprintf(' - Mode: Single Simulation (Parallel Outputs) \n');
    fprintf(' - Passive: out.yout{2} | Active: out.yout{1} \n');
    fprintf('========================================================\n');

    %% Geometry, physical parameters and tuning settings
    % k_TMD* denotes the retuned physical spring; k_pa* is the passive target.
    P = struct( ...
               'g', 9.81, 'h_R', 147.082, 'h_p', 14.94, 'h_T', 15, 'z', 14.94, 'L', 51.75, ...
               'm_N', 9.49e5, 'm_T', 1.483e5, 'm_p', 1.7834e7, 'xi_TA', 0.01, 'I_p', 1.2507e10, ...
               'm_as', 9.64e6, 'm_asp', -1.01e8, 'I_a', 1.16e10, ...
               'c_s', 7.94e4, 'c_sp', -2.2e7, 'c_p', 5.80e8, ...
               'k_p', 2.60e9, 'k_h', 4.470e6, 'k_hp', 4.239e3, ...
               'k_mooringP', -0.76e8, 'k_mooringS', 7.965e4, 'k_mooringSP', 1.162e6, ...
               'coes', 9.225e5, 'coesp', 8.92e6, 'coep', 1.68e10, 'coeh', 2.296e6, ...
               'm_TMD1', 2.0e5, 'm_TMD2', 4.0e5, ...
               'k_TMD1', 177523, 'k_TMD2', 177523 * 2, ...
               'c_TMD1', 8717, 'c_TMD2', 8717 * 2, ...
               'k_pa', 10232, 'k_pa1', 20464, ... % Passive Parameters
               'f_des', 0.036 ...
              );

    % P.L and P.L1 correspond to the thesis lever arms L1 and L2.
    P.L1 = P.L / 2;

    wd = 2 * pi * P.f_des;

    static_deflection1 = P.m_TMD1 * P.g / P.k_TMD1;
    static_deflection2 = P.m_TMD2 * P.g / P.k_TMD2;

    [m_TA, m_1, m_2, m_3, m_4, k_TA] = tower_equiv(P.g, P.m_N);

    [A_o, B_o, E_o, D_o, C_V_o] = base3dof(P, m_TA, m_1, m_2, m_3, m_4, k_TA);

    [M_1, C_1, K_1] = full6dof_MCK(P, m_TA, m_1, m_2, m_3, m_4, k_TA);
    [A, B] = so2ss(M_1, C_1, K_1);
    M_Structure_Pitch = M_1(4, 4);

    K_pa = K_1;
    K_pa(5, 5) = P.k_pa;
    K_pa(6, 6) = P.k_pa1;
    [A_Pa, ~] = so2ss(M_1, C_1, K_pa);

    n = size(M_1, 1);
    O = zeros(n);
    I = eye(n);
    E = [I O];
    E1 = [O I];
    D = zeros(size(E1, 1), size(B, 2));

    C_V = [0, 0, 0, 0, 0, 0
           0, -P.coes, 0, -P.coesp, 0, 0
           0, 0, -P.coeh, 0, 0, 0
           0, -P.coesp, 0, -P.coep, 0, 0
           0, 0, 0, 0, 0, 0
           0, 0, 0, 0, 0, 0];

    F_G = [0; 0; 0; 0; -P.m_TMD1 * P.g; -P.m_TMD2 * P.g];

    x_eq1 = -P.m_TMD1 * P.g / P.k_TMD1;
    x_eq2 = -P.m_TMD2 * P.g / P.k_TMD2;

    wc_hz = 1.0;
    wc = 2 * pi * wc_hz;

    X2 = [0; 0; 0; 0; -static_deflection1; -static_deflection2; zeros(6, 1)];
    X4 = [0; 0; 0; 0; -P.m_TMD1 * P.g / P.k_TMD1; -P.m_TMD2 * P.g / P.k_TMD2];

    X2_pa = [0; 0; 0; 0; -P.m_TMD1 * P.g / P.k_pa; -P.m_TMD2 * P.g / P.k_pa1; zeros(6, 1)];
    X4_pa = [0; 0; 0; 0; -P.m_TMD1 * P.g / P.k_pa; -P.m_TMD2 * P.g / P.k_pa1];

    Xo_o = zeros(6, 1);
    Xo_o1 = zeros(3, 1);

    fnames = fieldnames(P);
    for i = 1:numel(fnames)
        assignin('base', fnames{i}, P.(fnames{i}));
    end

    assignin('base', 'A', A);
    assignin('base', 'B', B);
    assignin('base', 'E', E);
    assignin('base', 'E1', E1);
    assignin('base', 'D', D);
    assignin('base', 'F_G', F_G);
    assignin('base', 'C_V', C_V);

    assignin('base', 'A_Pa', A_Pa);

    assignin('base', 'A_o', A_o);
    assignin('base', 'B_o', B_o);
    assignin('base', 'E_o', E_o);
    assignin('base', 'D_o', D_o);
    assignin('base', 'C_V_o', C_V_o);

    assignin('base', 'x_eq1', x_eq1);
    assignin('base', 'x_eq2', x_eq2);
    assignin('base', 'wc', wc);
    assignin('base', 'wc_hz', wc_hz);
    assignin('base', 'X2', X2);
    assignin('base', 'X4', X4);
    assignin('base', 'X2_pa', X2_pa);
    assignin('base', 'X4_pa', X4_pa);
    assignin('base', 'Xo_o', Xo_o);
    assignin('base', 'Xo_o1', Xo_o1);

    %% Load the model and initialise the control subsystem
    mdl = 'PlatformDamper_Simulink_DA_Passive_MassRatio_Damping';
    load_system(mdl);

    subsys = [mdl '/Control subsystem'];
    pblk = @(s) [subsys '/' s];

    Data_Active = struct();
    Data_Passive = struct('t', [], 'Pitch', [], 'f', [], 'A', [], 'rms', 0);

    fprintf(' Running Passive TMD Simulation...\n');

    set_param([mdl, '/Passive SS'], 'InitialCondition', 'X2_pa');
    set_param([mdl, '/Results Output_p/Integrator'], 'InitialCondition', 'X4_pa');

    %% Sweep the equivalent-mass coefficient
    alpha_list = [0.5:0.1:1.5];
    rms_values = zeros(size(alpha_list));

    set_param(pblk('Constant_xeq'),  'Value', sprintf('%.15g', x_eq1));
    set_param(pblk('Constant_xeq1'), 'Value', sprintf('%.15g', x_eq2));

    set_param([mdl, '/Active SS'], 'InitialCondition', 'X2');
    set_param([mdl, '/Results Output/ActiveIntegrator'], 'InitialCondition', 'X4');

    fprintf('\n Alpha | Mass Ratio (mu) | Opt Xi | Ca1 (Active Damp) | RMS (deg)\n');
    fprintf(' ------------------------------------------------------------------------------\n');

    for k = 1:numel(alpha_list)

        alpha_mass = alpha_list(k);

        assignin('base', 'A', A); % Reset to Active Matrix A

        % Target mass, stiffness, damping and AVD feedback gains.
        m_eq1 = alpha_mass * P.m_TMD1;
        m_eq2 = alpha_mass * P.m_TMD2;
        k_eq1 = m_eq1 * wd^2;
        k_eq2 = k_eq1 * 2;

        mu1 = (m_eq1 * P.L^2)  / M_Structure_Pitch;
        mu2 = (m_eq2 * P.L1^2) / M_Structure_Pitch;

        xi_opt1 = sqrt(3 * mu1 / (8 * (1 + mu1)^3));
        xi_opt2 = sqrt(3 * mu2 / (8 * (1 + mu2)^3));

        c_target1 = 2 * xi_opt1 * sqrt(m_eq1 * k_eq1);
        c_target2 = 2 * xi_opt2 * sqrt(m_eq2 * k_eq2);

        Ka1 = m_eq1 - P.m_TMD1;
        Ka2 = m_eq2 - P.m_TMD2;
        Kp1 = k_eq1 - P.k_TMD1;
        Kp2 = k_eq2 - P.k_TMD2;

        Ca1 = c_target1 - P.c_TMD1;
        Ca2 = c_target2 - P.c_TMD2;

        assignin('base', 'alpha_mass', alpha_mass);
        assignin('base', 'Ka1', Ka1);
        assignin('base', 'Ka2', Ka2);
        assignin('base', 'Kp1', Kp1);
        assignin('base', 'Kp2', Kp2);
        assignin('base', 'Ca1', Ca1);
        assignin('base', 'Ca2', Ca2);

        set_param(pblk('Gain KP1'), 'Gain', sprintf('%.15g', Kp1));
        set_param(pblk('Gain KP2'), 'Gain', sprintf('%.15g', Kp2));
        set_param(pblk('Gain KA1'), 'Gain', sprintf('%.15g', Ka1));
        set_param(pblk('Gain KA2'), 'Gain', sprintf('%.15g', Ka2));
        set_param(pblk('Gain CA1'), 'Gain', sprintf('%.15g', Ca1));
        set_param(pblk('Gain CA2'), 'Gain', sprintf('%.15g', Ca2));

        out = sim([mdl '.slx']);
        assignin('base', 'out', out);

        t = out.tout(:);
        Pitch = out.yout{1}.Values.Data(:, 2); % Active uses yout{1}

        t_uni = t(1):(1 / 50):t(end); % Fixed Fs=50Hz
        Pitch_uni = interp1(t, Pitch, t_uni, 'linear');

        rms_val = rms(Pitch_uni);
        rms_values(k) = rms_val;

        fprintf(' %.2f | mu1=%.4f (TMD1) | %.4f | %12.2e | %.4f\n', ...
                alpha_mass, mu1, xi_opt1, Ca1, rms_val);

        [f, Amp] = fft_amp_optimized(Pitch_uni, 50);

        Data_Active(k).t = t;
        Data_Active(k).Pitch = Pitch;
        Data_Active(k).f = f;
        Data_Active(k).A = Amp;
        Data_Active(k).Legend = sprintf('\\alpha=%.1f (Active)', alpha_mass);

        if isempty(Data_Passive.t) && numel(out.yout) >= 2
            try
                Pitch_pa_raw = out.yout{2}.Values.Data(:, 2);
                Data_Passive.t = t;
                Data_Passive.Pitch = Pitch_pa_raw;

                Pitch_pa_uni = interp1(t, Pitch_pa_raw, t_uni, 'linear');
                Data_Passive.rms = rms(Pitch_pa_uni);

                [f_pa, A_pa] = fft_amp_optimized(Pitch_pa_uni, 50);
                Data_Passive.f = f_pa;
                Data_Passive.A = A_pa;

                fprintf(' [Passive Data] Extracted from yout{2}. RMS=%.4f\n', Data_Passive.rms);
            catch ME
                fprintf(' [Warning] Passive extraction failed: %s\n', ME.message);
            end
        end
    end

    styles = {'-', '--', ':', '-.'}; % keep original variable (harmless)
    colors = lines(length(alpha_list));

    freq_vec = Data_Active(1).f;

    baseStyles = {'-', '--', ':', '-.'};
    markers    = {'none', 'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x', '+', '*'};

    figure('Name', 'Time Domain Analysis', 'Color', 'w', 'Position', [50, 400, 700, 500]);
    hold on;

    for k = 1:numel(alpha_list)

        nLS = numel(baseStyles);
        nMK = numel(markers);

        style_idx = mod(k - 1, nLS) + 1;

        if k <= nLS
            mk = 'none'; % line style only
        else
            marker_idx = floor((k - 1) / nLS);
            marker_idx = mod(marker_idx, nMK - 1) + 2; % skip 'none'
            mk = markers{marker_idx};
        end

        h = plot(Data_Active(k).t, Data_Active(k).Pitch, ...
                 'Color', colors(k, :), ...
                 'LineStyle', baseStyles{style_idx}, ...
                 'Marker', mk, ...
                 'LineWidth', 1.5, ...
                 'DisplayName', Data_Active(k).Legend);

        if ~strcmp(mk, 'none')
            try
                set(h, 'MarkerIndices', round(linspace(1, numel(Data_Active(k).t), 12)));
            catch
                idx = round(linspace(1, numel(Data_Active(k).t), 12));
                plot(Data_Active(k).t(idx), Data_Active(k).Pitch(idx), mk, ...
                     'Color', colors(k, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');
            end
        end
    end

    if ~isempty(Data_Passive.t)
        plot(Data_Passive.t, Data_Passive.Pitch, 'k--', 'LineWidth', 2, 'DisplayName', 'Passive TMD');
    end

    grid on;
    box on;
    xlabel('Time (s)');
    ylabel('Pitch (deg)');
    title('Time Domain: Active vs Passive', 'FontSize', 14);
    legend('Location', 'bestoutside');

    figure('Name', 'Frequency Domain Analysis', 'Color', 'w', 'Position', [800, 400, 700, 500]);
    hold on;

    for k = 1:numel(alpha_list)

        nLS = numel(baseStyles);
        nMK = numel(markers);

        style_idx = mod(k - 1, nLS) + 1;

        if k <= nLS
            mk = 'none';
        else
            marker_idx = floor((k - 1) / nLS);
            marker_idx = mod(marker_idx, nMK - 1) + 2;  % 2..nMK, skip 'none'
            mk = markers{marker_idx};
        end

        x = Data_Active(k).f;
        y = 20 * log10(max(Data_Active(k).A, eps));

        semilogx(x, y, ...
                 'Color', colors(k, :), ...
                 'LineStyle', baseStyles{style_idx}, ...
                 'LineWidth', 1.5, ...
                 'DisplayName', Data_Active(k).Legend);

        if ~strcmp(mk, 'none')
            N = numel(x);

            nPts = 12;
            idx = unique(round(logspace(log10(2), log10(N), nPts)));
            idx(idx < 2) = [];           % 确保不取到 1
            idx(idx > N) = [];           % 防越界

            semilogx(x(idx), y(idx), ...
                     'LineStyle', 'none', ...
                     'Marker', mk, ...
                     'MarkerSize', 6, ...
                     'Color', colors(k, :), ...
                     'HandleVisibility', 'off');  % 不重复进 legend
        end
    end

    if ~isempty(Data_Passive.f)
        plot(Data_Passive.f, 20 * log10(max(Data_Passive.A, eps)), 'k--', ...
             'LineWidth', 2.5, 'DisplayName', 'Passive TMD');
    end

    set(gca, 'XScale', 'log'); % <-- keep log axis
    grid on;
    box on;
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (dB)');
    title('Frequency Response: Active vs Passive', 'FontSize', 14);
    xlim([0.01 0.2]);
    legend('Location', 'bestoutside');

    figure('Name', 'RMS Comparison', 'Color', 'w', 'Position', [400, 50, 600, 300]);
    hold on;

    b = bar(alpha_list, rms_values);
    b.FaceColor = 'flat';
    b.CData = colors;

    if Data_Passive.rms > 0
        yline(Data_Passive.rms, 'k--', 'Passive RMS', 'LineWidth', 2, ...
              'LabelHorizontalAlignment', 'left', 'FontSize', 10);
    end

    xlabel('Alpha (\alpha)');
    ylabel('Pitch RMS (deg)');
    title('RMS Performance Comparison', 'FontSize', 14);
    grid on;
    box on;

end

% Spectral post-processing.
function [f, A] = fft_amp_optimized(x, fs)

    x = detrend(x(:), 'constant');
    N = numel(x);
    w = hann(N);
    x_win = x .* w;
    nfft = 2^nextpow2(N);
    X = fft(x_win, nfft);
    correction_factor = 2;
    A = 2 * abs(X(1:nfft / 2 + 1)) / N * correction_factor;
    f = (fs * (0:nfft / 2) / nfft).';
    A = smoothdata(A, 'movmean', 5);
    f = f(:);
    A = A(:);

end

% Equivalent tower quantities.
function [m_TA, m_1, m_2, m_3, m_4, k_TA] = tower_equiv(g, m_N)

    H_T = 129.582;
    syms z_T;
    mu = 0.0084 * z_T^3 - 1.077 * z_T^2 - 171.5 * z_T + 2.92e04;
    EI = 1.905e06 * z_T^3 - 2.47e08 * z_T^2 - 5.208e10 * z_T + 6.851e12;
    Phi = 0.9414 * (z_T / H_T)^2 + 0.3468 * (z_T / H_T)^3 - 1.073 * (z_T / H_T)^4 + ...
        1.3139 * (z_T / H_T)^5 - 0.5289 * (z_T / H_T)^6;

    m_TA = double(int(mu * Phi^2, z_T, 0, H_T));
    m_1  = double(int(mu * Phi, z_T, 0, H_T));
    m_2  = double(int(mu * z_T * Phi, z_T, 0, H_T));
    m_3  = double(int(mu * z_T, z_T, 0, H_T));
    m_4  = double(int(mu * z_T^2, z_T, 0, H_T));

    D2 = diff(Phi, z_T, 2);
    D1 = diff(Phi, z_T, 1);

    f1 = int(EI * D2^2, 0, H_T);
    f2 = int(mu, z_T, H_T);
    f3 = int(g * (m_N + f2) * D1^2, 0, H_T);
    k_TA = double(f1 - f3);

end

function [A_o, B_o, E_o, D_o, C_V_o] = base3dof(P, m_TA, m_1, m_2, m_3, m_4, k_TA)

    M_o = [P.m_N + m_TA, P.m_N + m_1, P.m_N * P.h_R + m_2
           P.m_N + m_1, P.m_N + P.m_p + P.m_T, (P.m_N * P.h_R - P.m_p * P.z + m_3)
           P.m_N * P.h_R + m_2, (P.m_N * P.h_R - P.m_p * P.z + m_3), ...
           P.m_N * P.h_R^2 + P.I_p + m_4 + P.m_p * P.z^2];

    C_o = [2 * P.xi_TA * sqrt(m_TA * k_TA) 0 0; 0 0 0; 0 0 0];

    K_o = [k_TA, 0, -(P.m_N + P.m_T) * P.g
           0, 0, 0
           -(P.m_N + P.m_T) * P.g, 0, -(P.m_N * P.h_R + m_3 - P.m_p * P.z) * P.g];

    M_add_o = [0 0 0; 0 P.m_as P.m_asp; 0 P.m_asp P.I_a];

    C_add_o = [0 0 0; 0 P.c_s P.c_sp; 0 P.c_sp P.c_p];

    K_add_o = [0 0 0; 0 0 0; 0 0 P.k_p];

    K_moor_o = [0 0 0; 0 P.k_mooringS P.k_mooringSP; 0 P.k_mooringSP P.k_mooringP];

    M_2 = M_o + M_add_o;
    C_2 = C_o + C_add_o;
    K_2 = K_o + K_add_o + K_moor_o;

    O3 = zeros(3);
    I3 = eye(3);

    A_o = [O3 I3; -M_2 \ K_2 -M_2 \ C_2];
    B_o = [O3; M_2 \ eye(3)];
    E_o = [O3 I3];
    D_o = zeros(3);

    C_V_o = [0 0 0; 0 -P.coes -P.coesp; 0 -P.coesp -P.coep];

end

% Six-DOF matrices in the Appendix C coordinate convention.
function [M_1, C_1, K_1] = full6dof_MCK(P, m_TA, m_1, m_2, m_3, m_4, k_TA)

    % Lever-arm names in Appendix C: L1 for the front TMD, L2 for the pair.
    L1 = P.L;
    L2 = P.L1;
    mSigma = P.m_N + P.m_p + P.m_T + P.m_TMD1 + P.m_TMD2;
    rSigma = P.m_TMD1 * L1 - P.m_TMD2 * L2;
    ISigma = P.m_N * P.h_R^2 + P.I_p + m_4 + P.m_p * P.z^2 ...
           + P.m_TMD1 * L1^2 + P.m_TMD2 * L2^2;
    mTp = P.m_N * P.h_R + m_2;
    mSp = P.m_N * P.h_R - P.m_p * P.z + m_3;

    M_pp = [P.m_N + m_TA, P.m_N + m_1, 0, mTp
            P.m_N + m_1, mSigma, 0, mSp
            0, 0, mSigma, rSigma
            mTp, mSp, rSigma, ISigma];
    M_pa = [0, 0
            0, 0
            P.m_TMD1, P.m_TMD2
            P.m_TMD1 * L1, -P.m_TMD2 * L2];
    M_aa = diag([P.m_TMD1, P.m_TMD2]);
    M = [M_pp, M_pa; M_pa.', M_aa];

    C = diag([2 * P.xi_TA * sqrt(m_TA * k_TA), 0, 0, 0, P.c_TMD1, P.c_TMD2]);

    K = [k_TA, 0, 0, -(P.m_N + P.m_T) * P.g, 0, 0
         0, 0, 0, 0, 0, 0
         0, 0, 0, 0, 0, 0
         -(P.m_N + P.m_T) * P.g, 0, 0, -(P.m_N * P.h_R + m_3 - P.m_p * P.z) * P.g, 0, 0
         0, 0, 0, 0, P.k_TMD1, 0
         0, 0, 0, 0, 0, P.k_TMD2];

    M_add = zeros(6);
    M_add(2, 2) = P.m_as;
    M_add(2, 4) = P.m_asp;
    M_add(3, 3) = 2.48e7;
    M_add(4, 2) = P.m_asp;
    M_add(4, 4) = P.I_a;

    C_add = zeros(6);
    C_add(2, 2) = P.c_s;
    C_add(2, 4) = P.c_sp;
    C_add(3, 3) = 3e6;
    C_add(4, 2) = P.c_sp;
    C_add(4, 4) = P.c_p;

    K_add = zeros(6);
    K_add(3, 3) = P.k_h;
    K_add(4, 4) = P.k_p;
    K_add(3, 4) = P.k_hp;
    K_add(4, 3) = P.k_hp;

    K_mooring = zeros(6);
    K_mooring(2, 2) = P.k_mooringS;
    K_mooring(2, 4) = P.k_mooringSP;
    K_mooring(3, 3) = 5.9e4;
    K_mooring(4, 2) = P.k_mooringSP;
    K_mooring(4, 4) = P.k_mooringP;

    M_1 = M + M_add;
    C_1 = C + C_add;
    K_1 = K + K_add + K_mooring;

end

function [A, B] = so2ss(M, C, K)

    n = size(M, 1);
    O = zeros(n);
    I = eye(n);
    A = [O I; -M \ K -M \ C];
    B = [O; M \ eye(n)];

end
