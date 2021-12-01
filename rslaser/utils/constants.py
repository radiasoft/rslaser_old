# -*- coding: utf-8 -*-
u"""Physical and numerical constants
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""
import math
import scipy.constants as const

TWO_PI = 2 * math.pi
RT_TWO_PI = math.sqrt(2*math.pi)
RT_2_OVER_PI = math.sqrt(2/math.pi)

C_SQ = const.c**2
C_INV  = 1./const.c
MKS_FACTOR = 1./(4.*math.pi*const.epsilon_0)
M_E_EV = const.m_e * C_SQ / (-const.e)
