function model = three_dof_model()
    % Author: Yu Gao
    % Purpose: Assemble the tower fore-aft, platform surge and platform pitch ROM.
    % THREE_DOF_MODEL
    % Model matrices follow Appendix A of the thesis.
    % State: X = [q_T; q_s; q_p; dq_T; dq_s; dq_p]. Pitch is in radians.
    % Outputs: model.M, C, K, A, B, rhs and integration options.
    % M, C and K retain the common 1/1000 scaling used by the source.
    % Use ode45(model.rhs, time_vector, initial_state, model.options).
    % Input data: none. The caller supplies time and initial state.
    % Requires MATLAB and Symbolic Math Toolbox.

    %% Geometry and model parameters
    syms z_T;
    h_R = 147.082; % RNA CoM height above the platform reference point [m].
    H_T = 129.582; % Tower height [m].
    z = 14.94; % Platform CoM depth below the platform reference point [m].
    g = 9.81;
    m_T = 1.483e5; % kg / tower
    m_N = 9.91e5; % RNA mass [kg].
    m_p = 1.7854e7;  % kg /platform mass
    xi_TA = 0.01; % Tower damping ratio specified in Appendix A.
    I_p = 1.2507 * 10^10; % kg m2 /mass moment of inertia of platform
    m_as = 9.64 * 10^6;  % kg / Added mass for platform surge
    I_a = 1.16 * 10^10;  % kg m2 / Added mass for platform pitch
    m_asp = -1.01 * 10^8; % kg m / Coupling effects of added mass between surge and pitch
    c_s = 9.94e5;  % Linear surge damping coefficient
    c_sp = -2.25e7; % coupled damping value between surge and platform pitch
    c_p = 5.60e8; % damping in pitch motion
    k_p = 2.49e9; % rotational stiffness of platform
    k_mooringS = 7.965e4;
    k_mooringSP = 1.162e6;
    k_mooringP = 2.65e8;

    %% Equivalent tower properties
    mu = 0.0084 * z_T^3 - 1.077 * z_T^2 - 171.5 * z_T + 2.72e04; % mass per length
    EI = 1.905e06 * z_T^3 - 2.47e08 * z_T^2 - 5.208e10 * z_T + 6.851e12; % tower bending stiffness
    % First tower fore-aft mode shape.
    Phi_TA = 0.9414 * (z_T / H_T)^2 + 0.3468 * (z_T / H_T)^3 ...
             - 1.073 * (z_T / H_T)^4 + 1.3139 * (z_T / H_T)^5 ...
             - 0.5289 * (z_T / H_T)^6;
    fun1 = mu * Phi_TA^2;
    m_TA = double(int(fun1, z_T, 0, H_T));
    fun2 = mu * Phi_TA;
    m_1 = double(int(fun2, z_T, 0, H_T));
    fun3 = mu * (z_T) * Phi_TA;
    m_2 = double(int(fun3, z_T, 0, H_T));
    fun4 = mu * (z_T);
    m_3 = double(int(fun4, z_T, 0, H_T));
    fun5 = mu * (z_T)^2;
    m_4 = double(int(fun5, z_T, 0, H_T));
    D2y = diff(Phi_TA, z_T, 2);
    Dy = diff(Phi_TA, z_T, 1);
    fun6 = EI * D2y^2;
    f1 = int(fun6, 0, H_T);
    fun7 = mu;
    f2 = int(fun7, z_T, H_T);
    fun8 = g * (m_N + f2) * Dy^2;
    f3 = int(fun8, 0, H_T);
    k_TA = double(f1 - f3);

    %% Structural, added-mass, damping and stiffness matrices
    M = [m_N + m_TA,      m_N + m_1,                   m_N * h_R + m_2
         m_N + m_1,       m_N + m_p + m_T,             (m_N * h_R - m_p * z + m_3)
         m_N * h_R + m_2, (m_N * h_R - m_p * z + m_3), m_N * h_R^2 + I_p + m_4 + m_p * z^2];
    C = [2 * xi_TA * sqrt(m_TA * k_TA) 0 0
         0 0 0
         0 0 0];
    K = [k_TA 0 -(m_N + m_T) * g
         0 0 0
         -(m_N + m_T) * g 0 -(m_N * h_R + m_3 - m_p * z) * g];
    M_add = [0 0 0
             0 m_as m_asp
             0 m_asp I_a];
    C_add = [0 0 0
             0 c_s c_sp
             0 c_sp c_p];
    K_add = [0 0 0
             0 0 0
             0 0 k_p];
    K_mooring = [0 0 0
                 0 k_mooringS k_mooringSP
                 0 k_mooringSP k_mooringP];
    M_1 = (M + M_add) / 1000;
    K_1 = (K + K_mooring + K_add) / 1000;
    C_1 = (C + C_add) / 1000;

    %% State-space representation
    O = zeros(3, 3);
    I = eye(3);
    A = [O I; -inv(M_1) * K_1 -inv(M_1) * C_1];
    B = [O; inv(M_1)];
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    model = struct('M', M_1, 'C', C_1, 'K', K_1, 'A', A, 'B', B, 'options', options);
    model.rhs = @(t, X) reducedmodel(t, A, B, X);
end

function dXdt = reducedmodel(~, A, B, X)
    % Generalised-load function retained from the source implementation.
    m_N = 1.017e3;
    g = 9.80655;
    F_g = m_N * g * 6.6673;

    coes = 9.225e5;
    coesp = 8.92e6;
    coep = 1.68e10;

    F1 = -(coes * X(5) * abs(X(5)) - coesp * X(6) * abs(X(6))) / 1000;
    F2 = -(coep * X(6) * abs(X(6)) - coesp * X(5) * abs(X(5))) / 1000;

    F_offset = [m_N * 12.1673 * g / 129.386; F1; F2 + F_g];
    dXdt = A * X + B * F_offset;
end
