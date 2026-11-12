function model = platform_tmd_model()
    % Author: Yu Gao
    % Purpose: Assemble the six-DOF platform-TMD ROM used by the accompanying Simulink driver.
    % PLATFORM_TMD_MODEL
    % Structural matrices follow Appendix C of the thesis.
    % Coordinates: [q_T; q_s; q_h; q_p; x_TMD1; x_TMD23].
    % TMD23 is the equivalent lateral pair; masses are 200,000 and 400,000 kg.
    % Outputs: model.M, C, K, A, B and the parameter structure model.parameters.
    % This function assembles the model; excitation and control are in folder 05.
    % Input data: none. Requires MATLAB and Symbolic Math Toolbox.

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
    [m_TA, m_1, m_2, m_3, m_4, k_TA] = tower_equiv(P.g, P.m_N);
    [M_1, C_1, K_1] = full6dof_MCK(P, m_TA, m_1, m_2, m_3, m_4, k_TA);
    [A, B] = so2ss(M_1, C_1, K_1);
    model = struct('M', M_1, 'C', C_1, 'K', K_1, 'A', A, 'B', B, 'parameters', P);
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
