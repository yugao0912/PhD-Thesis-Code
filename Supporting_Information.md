# Supporting Information: Computational Code

Author: Yu Gao

I provide the source code for the ROMs, TMD models and control implementations listed below. Each file identifies its author, function, inputs, outputs and software dependencies. Raw simulation histories, environmental input files, archived results and unpublished study datasets are not included.

| Folder | Function | Entry point and inputs |
|---|---|---|
| 01_Three_DOF_ROM | Assemble the tower fore-aft, platform surge and platform pitch model | `model = three_dof_model()`. Time and initial state are supplied by the caller for integration. |
| 02_Nacelle_TMD | Assemble the four-DOF nacelle-TMD model | `model = nacelle_tmd_model(k_TMD, c_TMD)`, with stiffness in N/m and damping in N s/m. |
| 03_QZS_TMD_Time_Domain | Compute nonlinear harmonic-frequency sweeps | `result = qzs_tmd_time_domain()`. Parameters and sweep settings are defined in the function. |
| 04_Platform_TMD_ROM | Assemble the six-DOF platform-TMD model | `model = platform_tmd_model()`. The equivalent lateral pair is represented by TMD23. |
| 05_Platform_TMD_Simulink | Compare uncontrolled, passive and AVD-controlled responses | Run `Run_Active_Damping_Optimal_FullModel` from this folder with its accompanying `.slx` file. |
| 06_OrcaFlex_AVD | Compute the three TMD control forces and their platform reactions | Load the six external functions through an OrcaFlex model configured with the specified objects and tags. |

The MATLAB model-assembly functions expose the matrices and, where applicable, the state-equation function for integration. The baseline, nacelle-TMD and platform-TMD matrix expressions follow Appendices A, B and C, respectively. The nacelle-TMD stiffness and damping are explicit inputs: the thesis pitch-tuned values are 2,570 N/m and 3,685 N s/m. No response data are supplied.

The QZS-TMD function retains the source time-domain integration and FFT calculations. It returns the primary-mass response and the absorber-response quantity calculated from `x1+r`, each normalised by `f/k_1`. It contains the requested time-domain analysis; the separate GA and Pareto studies are not included in this package.

The Simulink driver contains all three feedback terms and sweeps the equivalent-mass coefficient from 0.5 to 1.5. The retained `.slx` filename contains `DA`; its control subsystem also includes the velocity-feedback gains. The OrcaFlex external functions use the reference coefficient 0.8 and compute opposing control and platform-reaction forces.

I use MATLAB numerical and symbolic routines, Simulink simulation blocks, and the OrcFxAPI interface as software dependencies. Their built-in solvers and libraries are not part of my source-code contribution. MATLAB, Simulink, OrcaFlex and their libraries must be obtained separately. The `.slx` model was saved in MATLAB R2025a.
