# -*- coding: utf-8 -*-
"""
This folder will contain - for student's reference - expressions from Roloff & Matek. 
These expressions are only for internal use and e.g. for studying/answering the exam. 
These expressions may never and under no exception be distributed or shared to anyone else 
and serve only for the purpose of students that follow the KULeuven courses on dimensioning of machine parts 
By using these files you commit to adhering to this condition ! 
"""


import sympy as sp
from . import MySymbol
#from ..Units.Units import mm_


class GearDesign:

    
    MySymbolDict = {'omega_1'  : '[rad/s]' ,\
                    'i'        : '[]' ,\
                    'omega_2'  : '[rad/s]' ,\
                    'r_n1'     : '[mm]' ,\
                    'r_n2'     : '[mm]' ,\
                    'r_2'      : '[mm]' ,\
                    'r_1'      : '[mm]' ,\
                    'k_s'      : '[MPa]' ,\
                    'v'        : '[m/s] circumferential velocity 𝑑_1 ⋅ 𝜋 ⋅ 𝑛_1 ∕ 60 in m/s' ,\
                    'F_t'      : '[N] Tangential (circumferential) force' ,\
                    'b'        : '[mm] Tooth width' ,\
                    'd_1'      : '[mm] Pitch circle diameter' ,\
                    'u'        : '[] Tooth ration Z_1/Z_2' ,\
                    'T_2'      : '[Nm] Outgoing torque' ,\
                    'a'        : '[m] Shaft distance ' ,\
                    'n_s'      : '[1/min]' ,\
                    'P_out'    : '[W]' ,\
                    'P_in'     : '[W]' ,\
                    'P_1'      : '[W]' ,\
                    'P_2'      : '[W]' ,\
                    'P_v'      : '[W] Loss of power',\
                    'T_1'      : '[Nm] Incoming Torque' ,\
                    'i_tot'    : '[] Total gear ratio' ,\
                    'nu_t'     : '[] Total gear efficiency' ,\
                    'n_tot'    : '[] Total power loss = nu_ZTot + nu_DTot + nu_LTot' ,\
                    'nu_ZTot'  : '[] Power losses due to roller sliding of the tooth flanks' ,\
                    'nu_LTot'  : '[] Power losses due to bearing friction' ,\
                    'nu_DTot'  : '[] Power losses due to shaft seals' ,\
                    'nu_z'     : '[] Gear efficiency' ,\
                    'beta_1'   : '[__o] Tooth angle' ,\
                    'beta_2'   : '[__O] Tooth angle' ,\
                    'rho__1'   : '[__O] Wedge friction angle: for 𝜇 ≈ 0,05 … 0,01 and 𝛼𝑛 = 20° 𝜚__1 𝑏𝑒𝑐𝑜𝑚𝑒𝑠 1 ≈ 3°. . .6°' ,\
                    'gama_av'  : '[__o] Mean pitch angle of the worm (chapter 23)',\
                    }


    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))
            
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------


    def E20_1A_TransmissionGearRation(self):
        """ i = omega_1 / omega_2"""
        i = self.omega_1 / self.omega_2
        return i

    def E20_1B_TransmissionGearRation(self):
        """ i = r_n2 / r_n1"""
        i = self.r_n2 / self.r_n1
        return i

    def E20_1C_TransmissionGearRation(self):
        """ i = r_2 / r_1"""
        i = self.r_2 / self.r_1
        return i

    def E20_2_ForceSpeedFactorForRollerTransmissions(self):
        """ k_s/v = (3*(F_t/(b*d_1))*((u+1)/u))*(1/v) """
        p_gem = (3*(self.F_t/(self.b * self.d_1))*((self.u + 1)/self.u))*(1/self.v)
        return p_gem

    def E20_3_ForceSpeedFactorForHelicalGears(self):
        """ k_s/v = T_2/(a**3 * n_s """
        p_gem = self.T_2/(self.a**3 * self.n_s)
        return p_gem

    def E20_4A_TotalTransmissionEfficiency(self):
        """ nu_t = P_out / P_in"""
        p_gem = self.P_out / self.P_in
        return p_gem

    def E20_4B_TotalTransmissionEfficiency(self):
        """ nu_t = P_2 / P_1"""
        p_gem = self.P_2 / self.P_1
        return p_gem

    def E20_4C_TotalTransmissionEfficiency(self):
        """ nu_t = P_2 / (P_2 + P_v)"""
        p_gem = self.P_2 / (self.P_2 + self.P_v)
        return p_gem

    def E20_4D_TotalTransmissionEfficiency(self):
        """ nu_gt = (T_2 * omega_2) / (T_1 * omega_1)"""
        p_gem = (self.T_2 * self.omega_2) / (self.T_1 * self.omega_1)
        return p_gem

    def E20_4E_TotalTransmissionEfficiency(self):
        """ nu_t = T_2 / (T_1 * i_tot)"""
        p_gem = self.T_2 / (self.T_1 * self.i_tot)
        return p_gem

    def E20_5_TotalPowerLoss(self):
        """ n_tot = nu_ZTot + nu_DTot + nu_LTot"""
        p_gem = self.nu_DTot + self.nu_ZTot + self.nu_LTot
        return p_gem

    def E20_6A_BevelGearTransmissionsEfficiency(self):
        """ nu_z = (cos(beta_2 + rho__1) * cos(beta_1)) / (cos(beta_1 + rho__1) * cos(beta_2))"""
        p_gem = (sp.cos(self.beta_2 + self.rho__1) * sp.cos(self.beta_1)) / (sp.cos(self.beta_1 + self.rho__1) * sp.cos(self.beta_2))
        return p_gem

    def E20_6B_BevelGearTransmissionsEfficiency(self):
        """ nu_Z = tan(beta_1 - rho__1) / tan(beta_1)"""
        p_gem = sp.tan(self.beta_1 - self.rho__1) / sp.tan(self.beta_1)
        return p_gem

    def E20_7A_DrivingWormGearEfficiency(self):
        """ nu_z = tan(gama_av) / tan(gama_av + rho__1)"""
        p_gem = sp.tan(self.gama_av) / sp.tan(self.gama_av + self.rho__1)
        return p_gem

    def E20_7A_DrivenWormWheelGearEfficiency(self):
        """ nu_z__1 = tan(gama_av - rho__1) / tan(gama_av)"""
        p_gem = sp.tan(self.gama_av - self.rho__1) / sp.tan(self.gama_av)
        return p_gem