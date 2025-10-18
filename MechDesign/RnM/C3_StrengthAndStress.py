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
# from ..Units.Units import mm_


class StrengthAndStress:
    
    MySymbolDict = {'sigma_z'           :'[N/mm²]' ,\
                    'sigma_zn'          :'[N/mm²]' ,\
                    'sigma_n'           :'[N/mm²]' ,\
                    'sigma_max'         :'[N/mm²] Maximum stress' ,\
                    'sigma_nom'         :'[N/mm²] Nominal stress' ,\
                    'sigma_sigmaMax'    :'[N/mm²]' ,\
                    'sigma_nomSigma'    :'[N/mm²]' ,\
                    'alpha_k'           :'[] Form factor to get from tables' ,\
                    'alpha_kSigma'      :'[]' ,\
                    'alpha_kTao'        :'[]' ,\
                    'tao_taoMax'        :'[N/mm²]' ,\
                    'tao_nomTao'        :'[N/mm²]' ,\
                    'sigma_E'           :'[N/mm²]' ,\
                    'beta_k'            :'[] Dynamic notch factor, see table 3-9 or calculate from the notch shape, factor alpha_k',\
                    'sigma_w'           :'[N/mm²] Long-term acceptable alternating strength of uncut polish bar' ,\
                    'sigma_wk'          :'[N/mm²] Long-term acceptable alternating strength of notched bar',\
                    'n'                 :'[] Support factor for the unnotched or notched construction part' ,\
                    'n_0'               :'[] Support factor for the unnotched or notched construction part' ,\
                    'beta_kTest'        :'[] Experimentally determined notch factor, valid for the test rod diameter' ,\
                    'K_alpha'           :'[] Shape-dependent size factor of the part and test rod respectively,\ (see 3.5.1-3.), values according to table 3-11d' ,\
                    'K_alphaTest'       :'[] Shape-dependent size factor of the part and test rod respectively,\ (see 3.5.1-3.), values according to table 3-11d' ,\
                    'beta_k1'           :'[]' ,\
                    'beta_k2'           :'[]' ,\
                    'K_Dtd'             :'[] Construction factor for tension/compression' ,\
                    'beta_ktD'          :'[]' ,\
                    'K_g'               :'[] Geometric size factor, see table 3-11c' ,\
                    'K_0Sigma'          :'[] Surface coefficient, see table 3-10' ,\
                    'K_v'               :'[] Surface reinforcement factor, see table 3-12' ,\
                    'K_Ds'              :'[] Construction factor for shear' ,\
                    'beta_ks'           :'[]' ,\
                    'K_0Tao'            :'[] Surface coefficient, see table 3-10' ,\
                    'K_Db'              :'[] Construction factor for bending' ,\
                    'beta_kb'           :'[]' ,\
                    'k_Dt'              :'[] Construction factor for torsion' ,\
                    'beta_kt'           :'[]',\
                     }

    
    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))
            
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
 

    def E3_14A_StrengthReducingFormFactor(self):
        """ alpha_k =  sigma_max / sigma_nom"""
        alpha_k = self.sigma_max / self.sigma_nom
        return alpha_k

    def E3_14B_StrengthReducingFormFactor(self):
        """ alpha_kSigma = sigma_sigmaMax / sigma_nomSigma """
        alpha_kSigma = self.sigma_sigmaMax / self.sigma_nomSigma
        return alpha_kSigma

    def E3_14C_StrengthReducingFormFactor(self):
        """ alpha_kTao = tao_taoMax / tao_nomTao"""
        alpha_kTao = self.tao_taoMax / self.tao_nomTao
        return alpha_kTao

    def E3_15A_DynamicNotchFactor(self):
        """ beta_k = sigma_w / sigma_wk"""
        beta_k = self.sigma_w / self.sigma_wk
        return beta_k

    def E3_15B_DynamicNotchFactor(self):
        """ p_gem = alpha_k / (n_0 * n)"""
        p_gem = self.alpha_k / (self.n_0 * self.n)
        return p_gem

    def E3_15C_DynamicNotchFactor(self):
        """ beta_k = beta_kTest * (K_alphaTest / K_alpha)"""
        beta_k = self.beta_kTest * (self.K_alphaTest / self.K_alpha)
        return beta_k

    def E3_15D_DynamicNotchFactor(self):
        """ beta_k <= 1 + (beta_k1 -1) + (beta_k2 -1)"""
        beta_k = 1 + (self.beta_k1 - 1) + (self.beta_k2 -1)
        return beta_k

    def E3_16A_ConstructionFactorTensionCompression(self):
        """ K_Dtd = ((beta_ktd / K_g) + (1 / K_0Sigma) -1) * (1 / K_v) """
        K_Dtd = ((self.beta_ktD / self.K_g) + (1 / self.K_0Sigma) -1) * (1 / self.K_v)
        return K_Dtd

    def E3_16B_ConstructionFactorShear(self):
        """ K_Ds = ((beta_ks / K_g) + (1 / K_0Tao -1) * (1 / K_v) """
        K_Ds = ((self.beta_ks / self.K_g) + (1 / self.K_0Tao) -1) * (1 / self.K_v)
        return K_Ds

    def E3_16C_ConstructionFactorBending(self):
        """ K_Db = ((beta_kb / K_g) + (1 / K_0Sigma -1) * (1 / K_v) """
        K_Db = ((self.beta_kb / self.K_g) + (1 / self.K_0Sigma) - 1) * (1 / self.K_v)
        return K_Db

    def E3_16D_ConstructionFactorTorsion(self):
        """ K_Dt = ((beta_kt / K_g) + (1 / K_0Tao -1) * (1 / K_v) """
        K_Dt = ((self.beta_kt / self.K_g) + (1 / self.K_0Tao) - 1) * (1 / self.K_v)
        return K_Dt

