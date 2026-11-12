# Author: Yu Gao
# Purpose: Apply the equal-and-opposite AVD reaction to the platform.
# PIDcontrollerforDis_Cons4_structure.py
# Application: Chapter 7, vertical platform-mounted TMD configuration.
# Feedback object: Constraint4; channels: z, z velocity, z acceleration.
# State access also expects the object named 6D buoy2 in the OrcaFlex model.
# Used tags: Mass, k_tmd, c_tmd, g value, coe.
# k_target and c_target are read for source compatibility; gains are computed below.
# Use the consistent mass, force, length and time units of the OrcaFlex model.
# Parameters: alpha = 0.8, mu = 0.01, target frequency = 0.036 Hz.
# Kp, Ca and Ka are displacement, velocity and acceleration feedback gains.
# The PID filename is retained for the model link; no integral feedback is used.
# Pair this file with the same-number control/reaction file.
# Requires OrcaFlex, OrcFxAPI and NumPy; the OrcaFlex case is not supplied.
# Input data: live model states and tags; no recorded time series are embedded.

import OrcFxAPI
import numpy as np
import math


class InertialForce(object):

    def Initialise(self, info):

        def GetParam(paramName, default=None):
            param = info.ModelObject.tags.get(paramName, None)
            if param is None:
                if default is None:
                    raise Exception(
                        "Parameter {} is required but is not included in "
                        "the object tags.".format(paramName)
                    )
                return default
            return param

        def GetFloatParam(paramName, default=None):
            return float(GetParam(paramName, default))

        self.Mass = GetFloatParam("Mass", 0.0)
        self.k_tmd = GetFloatParam("k_tmd", 0.0)
        self.k_target = GetFloatParam("k_target", 0.0)
        self.c_tmd = GetFloatParam("c_tmd", 0.0)
        self.c_target = GetFloatParam("c_target", 0.0)
        self.g = GetFloatParam("g value", 0.0)
        self.coe = GetFloatParam("coe", 0.0)
        self.alpha = 0.8
        self.mu = 0.01

    def Calculate(self, info):
        self.period = OrcFxAPI.Period(OrcFxAPI.pnInstantaneousValue)
        Damper = info.Model["6D buoy2"]
        Constraint4 = info.Model["Constraint4"]
        dis = Constraint4.TimeHistory("z", self.period)
        vel = Constraint4.TimeHistory("z velocity", self.period)
        acc = Constraint4.TimeHistory("z acceleration", self.period)

        # Equivalent mass and target stiffness.
        m_eq = self.alpha * self.Mass
        k_eq = m_eq * math.pow((2 * math.pi * 0.036), 2)

        xi = self.coe * math.sqrt(3 * self.mu / (8 * (1 + self.mu)))

        c_target = 2 * xi * math.sqrt(m_eq * k_eq)

        # Gains relative to the physical spring, damper and mass.
        Ca = c_target - self.c_tmd
        Kp = k_eq - self.k_tmd
        Ka = m_eq - self.Mass

        # Static gravity equilibrium along the vertical guide.
        xeq = -self.Mass * self.g / self.k_tmd
        info.Value = Kp * (dis - xeq) + Ca * vel + Ka * acc
