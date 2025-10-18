# -*- coding: utf-8 -*-
"""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
This file will contain - for student's reference - expressions from Roloff & Matek. 
These expressions are only for internal use and e.g. for studying/answering the exam. 
These expressions may never and under no exception be distributed or shared to anyone else 
and serve only for the purpose of students that follow the KULeuven courses on dimensioning of machine parts 
By using these files you commit to adhering to this condition ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Formula's of chapter C2: Tolerances


"""

import sympy as sp
from . import MySymbol
from ..Units.Units import mm_


class Tolerance:

    MySymbolDict = {'G_G'   : '[mm] Upper limit value',
                    'G_K'   : '[mm] Lower limit value',
                    'G_GB'  : '[mm] Upper limit value Hole',
                    'G_GA'  : '[mm] Upper limit value Shaft',
                    'G_KB'  : '[mm] Lower limit value Hole',
                    'G_KA'  : '[mm] Lower limit value Shaft',
                    'T'     : '[mm] tolerance field',
                    'T_A'   : '[mm] tolerance field Shaft',
                    'T_B'   : '[mm] tolerance field Hole',
                    'ES'    : '[mm] Upper limit deviation Hole',
                    'EI'    : '[mm] Lower limit deviation Hole',
                    'es'    : '[mm] Upper limit deviation Shaft',
                    'ei'    : '[mm] Lower limit deviation Shaft',
                    'N'     : '[mm] Nominal value',
                    'D_1'   : '[mm] Lower limit of nominal tolerance',
                    'D_2'   : '[mm] Upper limit of nominal tolerance',
                    'D'     : '[mm] gemometrical average lower and upper limit nominal value=MySymbol( D = sqrt(D_1*D_2))',
                    'i'     : '[E-6m] Tolerance-Unit Shaft',
                    'I'     : '[E-6m] Tolerance-Unit Hole',
                    'I_B'   : '[mm] Efective diamter Hole ',
                    'I_A'   : '[mm] Efective diameter Shaft',
                    'P'     : '[mm] fitting',
                    'P_max' : '[mm] maximum clearance',
                    'P_min' : '[mm] mimimum clearance',
                    'P_T'   : '[mm] fitting tolerance',
                    }

    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))
            
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

    def E2_1A_HoleUpperLimVal(self):
        """ G_GB = N + ES """
        G_GB = self.N + self.ES
        return G_GB

    def E2_1A_ShaftUpperLimVal(self):
        """ G_GA = N + es """
        G_GA = self.N + self.es
        return G_GA

    def E2_2A_HoleLowerLimVal(self):
        """ G_KB = N + EI """
        G_KB = self.N + self.EI
        return G_KB

    def E2_2A_ShaftLowerLimVal(self):
        """ G_KA = N + ei """
        G_KA = self.N + self.ei
        return G_KA

    def E2_3A_Tolerance(self):
        """ T = G_G - G_K """
        T = self.G_G - self.G_K
        return T

    def E2_3B1_HoleTolerance(self):
        """ T_B = G_GB - G_KB """
        T_B = self.G_GB - self.G_KB
        return T_B

    def E2_3B2_HoleTolerance(self):
        """ T_B = ES - EI """
        T_B = self.ES - self.EI
        return T_B

    def E2_3C1_ShaftTolerance(self):
        """ T_A = G_GA - G_KA """
        T_A = self.G_GA - self.G_KA
        return T_A

    def E2_3C2_ShaftTolerance(self):
        """ T_A = es - ei """
        T_A = self.es - self.ei
        return T_A

    def E2_4A_ToleranceUnit(self):
        """ i = (0.45*D**(1/3) + 0.001*D )*1e-3
            to be used for nominal values <= 500
        """
        i = (0.45*(self.D*mm_**2)**(1/3) + 0.001*self.D)*1e-3
        return i

    def E2_4B_ToleranceUnit(self):
        """ I = (0.004*D+2.1)*1e-3 
            to be used for numinal values 500 < N <= 3150 
        """
        I = (0.0004*self.D + 2.1*mm_)*1e-3
        return I

    def E2_4_h_geometricMiddleField(self):
        """ D = sqrt(D_1*D_2)"""
        D = sp.sqrt(self.D_1*self.D_2)
        return D

    def E2_5A_FittingGeneral(self):
        """ P = I_B - I_A """
        P = self.I_B - self.I_A
        return P

    def E2_5B1_MaximumClearance(self):
        """ P_max = G_GB - G_KA """
        P_max = self.G_GB - self.G_KA
        return P_max

    def E2_5B2_MaximumClearance(self):
        """ P_max = ES - ei """
        P_max = self.ES - self.ei
        return P_max

    def E2_5C1_MinimumClearance(self):
        """ P_max = G_KB - G_GA """
        P_min = self.G_KB - self.G_GA
        return P_min

    def E2_5C2_MinimumClearance(self):
        """ P_max = EI - es """
        P_min = self.EI - self.es
        return P_min

    def E2_6A1_FittingTollerance(self):
        """ P_T = P_max - P_min """
        P_T = self.P_max - self.P_min
        return P_T

    def E2_6A2_FittingTollerance(self):
        """ P_T = (G_GB - G_GA) - (G_GK-G_GA) """
        P_T = (self.G_GB - self.G_GA) - (self.G_KB-self.G_GA)
        return P_T

    def E2_6B1_FittingTollerance(self):
        """ P_T = T_B + T_A """
        P_T = self.T_B - self.T_A
        return P_T

    def E2_6B2_FittingTollerance(self):
        """ P_T = (ES - EI) + (es - ei) """
        P_T = (self.ES - self.EI) + (self.es - self.ei)
        return P_T
