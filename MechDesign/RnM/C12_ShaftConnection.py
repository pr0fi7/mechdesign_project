# -*- coding: utf-8 -*-
"""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
This file will contain - for student's reference - expressions from Roloff & Matek.
These expressions are only for internal use and e.g. for studying/answering the exam.
These expressions may never and under no exception be distributed or shared to anyone else
and serve only for the purpose of students that follow the KULeuven courses on dimensioning of machine parts
By using these files you commit to adhering to this condition !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Formula's of chapter C12: Shaft Connections


"""

import sympy as sp
from . import MySymbol


class ShaftConnection:

    # create dictionary of symbols to be used
    MySymbolDict = {
      'p_gem'   : '[N/mm²] average pressure',
      'F_t'     : '[N] tangential load',
      'K_lambda': '[] load distribution factor method C: K_lambda = 1;  method B: {tabel 12-2c}',
      'h'       : '[mm] height of the key {tabel 12-2a}',
      'hprime'  : "[mm] h' is the effective height ",
      'l'       : '[mm] length of the key / length of the press fit /length of cone',
      'lprime'  : '[mm] effective length of the key',
      'b'       : '[mm] curvature width {tabel 12-2a}',
      'n'       : '[] number of keys or number of teeth of a spline',
      'phi'     : '[] correction factor based on number of keys',
      'd'       : '[mm] diameter of the shaft for splines: {tabel 12-3a} or addition for weakend section due to press fit screws (multiple ring press fit)',
      'T_eq'    : '[Nm] equivalent transmitted torque',
      'T_nom'   : '[Nm] nominal transmitted torque',
      'K_A'     : '[] load factor {tabel 3-5}',
      'K_t'     : '[] technology factor {equation 3-7; tabel 3-11}',
      'F_tNom'  : '[N] nominal tangential force',
      'F_tMax'  : '[N] maximum tangential force',
      'T'       : '[Nm] Torque',
      'T_max'   : '[Nm] maximum static torque to be transmitted',
      'R_e'     : '[N/mm²] yield strength {tabel 1-1}',
      'S_V'     : '[] safety factor',
      'R_m'     : '[N/mm²] tensile strength {tabel 1-1 of 1-2}',
      'S_B'     : '[]',
      'p_toel'  : '[N/mm²] permissible pressure',
      'f_S'     : '[] support factor {tabel 12-2d}',
      'f_H'     : '[] hardness factor {tabel 12-2d}',
      'd_gem'   : '[mm] average profile diameter (D+d)/2 {tabel 12-3a}',
      'D'       : '[mm] largest profile diameter {tabel 12-3a}',
      'e_r'     : '[mm] arithmetic eccentricity (d_1 -d_2)/4 {eq 12..3}',
      'e_1'     : '[mm] eccentricity of polygon {tabel 12-5}',
      'd_r'     : '[mm] arithmetic theoretical diameter (d_2 + 2e_r) {eq 12.3}',
      'd_1'     : '[mm] dimension according {tabel 12-5 tabel 12-4}',
      'd_2'     : '[mm] dimension according {tabel 12-5}',
      'r'       : '[mm] dimension according {tabel 12-5}',
      's'       : '[mm] minimum shaft wall thickness',
      'c'       : '[] profile factor',
      'sigma_toel': '[N/mm²] to be calculated helper eq 12.4',
      'L'       : '[mm] shaft length',
      'd_3'     : '[mm] diamter 3 {tabel 12-4b}',
      'd_a1'    : '[mm] diamter a1 {tabel 12-4a}',
      'd_a2'    : '[mm] diamter a2 {tabel 12-4a}',
      'm'       : '[] modulus {tabel 12-4b}',
      'Q_U'     : '[] diamter ratio outer component',
      'D_F'     : '[mm] diamter: shared diamter or inner diameter Press Ring',
      'D_Uu'    : '[mm] outer diamter outer component',
      'Q_I'     : '[] dimater ratio inner component',
      'D_Ii'    : '[mm] inner diamter inner component',
      'sigma_tUi': '[N/mm²] tengential tension outer edge outer component',
      'sigma_rUi': '[N/mm²] radial tension inner edge outer component',
      'sigma_tIi': '[N/mm²] tengential tension inner edge inner component',
      'sigma_tIu': '[N/mm²] tengential tension outer edge inner component',
      'sigma_rIu': '[N/mm²] radial tension uiter edge inner component',
      'p_F'     : '[N/mm²] contact pressure',
      'sigma_vUi': '[N/mm²] tension  outer component',
      'sigma_vIi': '[N/mm²] tension  inner component',
      'sigma_vUiLim': '[N/mm²] tension limit outer component',
      'sigma_vIiLim': '[N/mm²] tension limit inner component',
      'R_eU'    : '[N/mm²] strain limit outer component',
      'R_eI'    : '[N/mm²] strain limit inner component',
      'S_pU'    : '[] safety factor plastic deformation outer component',
      'S_pI'    : '[] safety factor plastic deformation inner component',
      'F_S'     : '[N] resulting "tangential force"',
      'S_S'     : '[] slip safety factor {eq 12.8} ',
      'F_l'     : '[N] "tangential force" longitudinal component',
      'F_res'   : '"tangential force" combination of torque and longitudinal load',
      'p_Fmin'  : '[N/mm²] minimum required contact pressure',
      'A_F'     : '[mm²] contact area of connection',
      'mu'      : '[] friction coefficient {tabel 12.6a}',
      'E_U'     : '[N/mm²] elastic modulus outer component',
      'E_I'     : '[N/mm²] elastic modulus inner component',
      'nu_U'    : '[] poisson ratio outer component',
      'nu_I'    : '[] possion ratio inner component',
      'Z'       : '[mm] press fit',
      'Z_min'   : '[mm] minimum required press fit',
      'K'       : '[] helper variable {eq 12.12}',
      'epsilon_Ui': '[] relative strain outer component inner side',
      'epsilon_Iu': '[] relative strain inner component outer side',
      'G'       : '[mm] roughness decrease',
      'Rz_UI'   : '[µm] roughness outer component inner side',
      'Rz_Iu'   : '[µm] roughness inner component outer side',
      'S_nmin'  : '[mm] minimum permissible over size',
      'S_nmax'  : '[mm] maximum permissible over size',
      'Sprime_nmax': '[mm] effectively used maximum oversize',
      'p_Fmax'  : '[N/mm²] maximum permissible contact pressure',
      'Z_max'   : '[mm] maximum permissible press fit',
      'P_T'     : '[µm] fitting tolerance',
      'T_B'     : '[µm] tolerance field outer component',
      'T_A'     : '[µm] tolerance field inner component',
      'F_i'     : '[N] maximum pressing force',
      'pprime_Fmax': '[N/mm²] effective maximal contact pressure',
      'mu_i'    : '[] friction coeficient for press fitting {tabel 12-6a}',
      'T_U'     : '[°C] mounting temperature outer component',
      'T_I'     : '[°C] mounting temperature inner component',
      'T_room'  : '[°C] room temperature',
      'S_m'     : '[mm] minimum required assembling clearance',
      'alpha_U' : '[1/K] outer component expansion coefficient {tabel 12-6}',
      'alpha_I' : '[1/K] inner component expansion coefficient {tabel 12-6}',
      'n_g'     : '[rpm] critical speed at which the minimum required pressure is reached',
      'pprime_Fmin': '[N/mm²] effective minimum pressure @ N=0 rpm',
      'rho_U'   : '[kg/m³] outer component denisity',
      'T_n'     : '[Nm] transmissable torque as function of rotational speed',
      'N'       : '[rpm] rotationl speed of shaft',
      'C'       : '[] cone shape ratio or form factor {eq 12.36}',
      'alpha'   : '[°] angle of conical press fit',
      'D_1'     : '[mm] minimum diameter of conical press fit',
      'D_2'     : '[mm] maximum diameter of conical press fit',
      'a_min'   : '[mm] minimum displacement distance',
      'a_max'   : '[mm] maximum displacement distance',
      'D_Fgem'  : '[mm] average contact diameter',
      'S_s'     : '[] slip safety factor {1.2<->1.5}',
      'rho_i'   : '[°] friction angle press fitting',
      'T_Tab'   : '[Nm] Transmittable torque of a single pressure ring',
      'T_tot'   : '[Nm] Transmittable torque of f_n pressure rings',
      'T_res'   : '[Nm] resulting Torque for pressure rings',
      'F_a'     : '[N] Axial pressure pressure ring',
      'F_aTab'  : '[N] Axial pressure on a single pressure ring',
      'F_atot'  : '[N] Total axial pressure of f_n pressure rings',
      'f_n'     : '[] numer of pressure rings',
      'p_N'     : '[N/mm²] local surface pressure on shaft',
      'p_A'     : '[N/mm²] local surface pressure on hub',
      'C_A'     : '[] form factor for hollow shafts {0.8<->1.0}',
      'l_F'     : '[mm] length of clamp connection ',
      'F_N'     : '[N] required clamping force',
      'F_Kl'    : '[N] required clamping force / bolt',
      'l_1'     : '[mm] length pinching point to center',
      'l_2'     : '[mm] length pinching point to bolt center',
      'F_VM'    : '[N] Tensioning force of bolt',
          }
    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------    
                               
    def E12_1A_KeyAveragePressure(self):
        """ p_gem = F_t * K_Lambda / (hprime * lprime * n * phi) """
        p_gem = self.F_t * self.K_lambda / (self.hprime * self.lprime * self.n * self.phi)
        return p_gem
    
    def E12_1B_KeyAveragePressure(self):
        """ p_gem = 2 * T_eq  * K_Lambda / (d * hprime * lprime * n * phi) """
        p_gem = 2 * self.T_eq  * self.K_lambda / (self.d * self.hprime * self.lprime * self.n * self.phi)
        return p_gem
    
    def E12_1_hA_DynamicForce(self):
        """ F_t = K_A * F_tNom """
        F_t = self.K_A * self.F_tNom
        return F_t
    
    def E12_1_hB_StaticForce(self):
        """ F_t = F_tMax """
        F_t = self.F_tMax
        return F_t
    
    def E12_1_hC_DynamicLoadTorque(self):
        """ T_eq = K_A * T_nom """
        T_eq = self.K_A * self.T_nom
        return T_eq
    
    def E12_1_hD_StaticLoadTorque(self):
        """ T_eq = T_max """
        T_eq =  self.T_max
        return T_eq
    
    def E12_1_hE_PermissiblePressure(self):
        """ p_toel = R_e / S_V """
        p_toel = self.R_e / self.S_V
        return p_toel
    
    def E12_1_hF_PermissiblePressure(self):
        """ p_toel = R_m / S_B """
        p_toel = self.R_m / self.S_B
        return p_toel

    def E12_1_hG_PermissiblePressure(self):
        """ p_toel = f_S * f_H * R_e / S_V """
        p_toel = self.f_S * self.f_H * self.R_e / self.S_V
        return p_toel

    def E12_1_hH_PermissiblePressure(self):
        """ p_toel = f_S * R_m / S_B """
        p_toel = self.f_S * self.R_e / self.S_V
        return p_toel       
    
    def E12_1_hI_KeyEffectiveLength(self):
        """ lprime = l - b """
        lprime = self.l - self.b
        return lprime
    
    def E12_1_hJ_KeyEffectiveHeight(self):
        """ hprime = 0.45 * h """
        hprime = 0.45 * self.h
        return hprime    
    
    def E12_2_SplineAveragePressure(self):
        """ p_gem = 2 * T *1E3 / (d_gem * L * hprime * 0.75 * n) """
        p_gem = 2 * self.T*1E3  / (self.d_gem * self.L * self.hprime * 0.75 * self.n)
        return p_gem
    
    def E12_2_hA_d_gem(self):
        """ d_gem = (D +d)/2 """
        d_gem = (self.D + self.d) / 2
        return d_gem
    
    def E12_2_hB_h_1(self):
        """ hprime = 0.4 * (D -d) """
        hprime = 0.4 * (self.D - self.d)
        return hprime
    
    def E12_2_hC_h_1_DIN5480(self):
        """ hprime = 0.5 * (d_a1 -(d_a2 +0.16*m))"""
        hprime = 0.5 * (self.d_a1 -(self.d_a2 +0.16*self.m))
        return hprime
    
    def E12_2_hD_h_1_DIN5481(self):
        """ hprime = 0.5*(d_3-d_1) """
        hprime = 0.5*(self.d_3-self.d_1)
        return hprime
    
    def E12_3A_PolygonAveragePressure(self):
        """ p_gem = T *1E3 / (lprime*(0.75*pi*e_1*d_1+0.05*d_1**2)) """
        p_gem = self.T*1E3 / (self.lprime * (0.75 -sp.pi * self.e_1 * self.d_1 + 0.05*self.d_1**2))
        return p_gem
    
    def E12_3B_PolygonAveragePressure(self):
        """ p_gem = T*1E3 / (lprime * (pi *e_r *d_r + 0.05*d_r**2)) """
        p_gem = self.T*1E3 / (self.lprime * (sp.pi * self.e_r * self.d_r + 0.05* self.d_r**2))
        return p_gem
    
    def E12_3_hA_T(self):
        """ T_eq = K_A * T_nom """
        T_eq = self.K_A * self.T_nom
        return T_eq
    
    def E12_3_hB_e_r(self):
        """ e_r = (d_1 - d_2) / 4 """
        e_r = (self.d_1 - self.d_2)/4
        return e_r
    
    def E12_3_hC_d_r(self):
        """ d_r = d_2 + 2*e_r """
        d_r = self.d_2 + 2*self.e_r
        return d_r
    
    def E12_4_MinShaftThickness(self):
        """ s = c * sp.sqrt( T*1E3 / (sigma_toel * L)) """
        s = self.c * sp.sqrt(self.T*1E3/(self.sigma_toel*self.L))
        return s
    
    def E12_4_hA_T(self):
        """ T_eq = K_A * T_nom """
        T_eq = self.K_A * self.T_nom
        return T_eq
    
    def E12_4_hB_sigma_toel(self):
        """ sigma_toel = R_e / S_V """
        sigma_toel = self.R_e / self.S_V
        return sigma_toel
    
    def E12_4_hC_sigma_toel(self):
        """ sigma_toel = R_m / S_B """
        sigma_toel = self.R_m / self.S_B
        return sigma_toel
    
    def E12_5_hA_Q_U(self):
        """ Q_U = (D_F/D_Uu) """
        Q_U = self.D_F / self.D_Uu
        return Q_U

    def E12_5_hB_Q_I(self):
        """ Q_I = (D_Ii /D_F) """
        Q_I = self.D_Ii / self.D_F
        return Q_I
    
    def E12_5A_sigma_tUi(self):
        """ sigma_tUi = p_F * (1+Q_U**2)/(1-Q_U**2) """
        sigma_tUi = self.p_F * (1+ self.Q_U**2)/(1 - self.Q_U**2)
        return sigma_tUi
    
    def E12_5B_sigma_tUu(self):
        """ sigma_tUu = p_F * (1+Q_U**2)/(1-Q_U**2) - p_F"""
        sigma_tUu = self.p_F * (1+ self.Q_U**2)/(1 - self.Q_U**2) - self.p_F
        return sigma_tUu
    
    def E12_6_hA_Q_I(self):
        """ Q_I = D_Ii / D_F """
        Q_I = self.D_Ii / self.D_F
        return Q_I

    def E12_6A_sigma_tIi(self):
        """ sigma_tIi = -1 * (p_F * (1+Q_I**2)/(1-Q_I**2) + p_F) """
        sigma_tIi = -1 * (self.p_F * (1+self.Q_I**2)/(self.Q_I) + self.p_F)
        return sigma_tIi
    
    def E12_6B_sigma_tIu(self):
        """ sigma_tIu = -1 * (p_F * (1+Q_I**2)/(1-Q_I**2)) """
        sigma_tIu = -1 * (self.p_F *  (1+self.Q_I**2)/(1-self.Q_I**2))
        return sigma_tIu
    
    def E12_7A_sigma_vUi(self):
        """ sigma_vUi = 2 * p_F / (1-Q_u**2) 
            the limit is to be set using sigma_vUiLim """
        sigma_vUi = 2* self.p_F / (1- self.Q_U**2)
        return sigma_vUi
    
    def E12_7B_sigma_vUiLim(self):
        """ sigma_vUiLim = 2/sqrt(3) * ReU/SpU """
        sigma_vUiLim = (2/sp.sqrt(3)) * (self.R_eU/self.S_pU)
        return sigma_vUiLim 
    
    def E12_7C_sigma_vIi(self):
        """ sigma_vIi = -1*2*p_F / (1-Q_I**2) 
            the limit is to be set using sigma_vIiLim"""
        sigma_vIi = sp.Abs(-2*self.p_F / (1-self.Q_I**2))
        return sigma_vIi
    
    def E12_7D_sigma_vIiLim(self):
        """ sigma_vIiLim = 2/sqrt(3) * R__eI/S_pI """
        sigma_vIiLim = (2/sp.sqrt(3)) * self.R_eI/self.S_pI
        return sigma_vIiLim
    
    def E12_8_hA_F_res(self):
        """ F_res = sqrt(F_t**2+F_l**2) """
        F_res = sp.sqrt(self.F_t**2 + self.F_l**2)
        return F_res
    
    def E12_8_F_S(self):
        """ F_S = S_S * K_A * F_res """
        F_S = self.S_S * self.K_A * self.F_res 
        return F_S
    
    def E12_9_p_Fmin(self):
        """ p_Fmin =  F_S / (A_F * mu) """
        p_Fmin = self.F_S / (self.A_F * self.mu)
        return p_Fmin

    def E12_9_hA_AF(self):
        """D *  pi * L"""
        A_F = self.D_F * sp.pi *  self.l_F
        return A_F
    
    def E12_10_epsilon_Ui(self):
        """ epsilon_Ui = (p_F / E_U) * ( (1+Q_U**2)/(1-Q_U**2) + nu_U) """
        epsilon_Ui = (self.p_F/self.E_U) * ( (1+self.Q_U**2)/(1-self.Q_U**2) + self.nu_U)
        return epsilon_Ui
    
    def E12_11_epsilon_Iu(self):
        """ epsilon_Iu = (p_F / E_I) *  ( (1+Q_I**2) / (1-Q_I**2) -nu_I ) """
        epsilon_Iu = (self.p_F / self.E_I) * ( (1+self.Q_I**2 ) / (1-self.Q_I**2) - self.nu_I)
        return epsilon_Iu
    
    def E12_12_K(self):
        """ K = E_U/E_I * ( (1+Q_I**2)/(1-Q_I**2)-nu_I) + (1+Q_U**2)/(1-Q_U**2) + nu_U """
        K = self.E_U/self.E_I * ( (1+self.Q_I**2)/(1-self.Q_I**2) - self.nu_I) + (1+self.Q_U**2)/(1-self.Q_U**2) + self.nu_U
        return K
    
    def E12_12_hA_Z(self):
        """ Z = D_F * (epsilon_Ui + epsilon_Iu) """
        Z = self.D_F * (self.epsilon_Ui + self.epsilon_Iu)
        return Z
        
    def E12_13_Z_min(self):
        """ Z_min = p_Fmin * D_F / E_U * K """
        z_min = self.p_Fmin * self.D_F / self.E_U * self.K
        return z_min
    
    def E12_14_G(self):
        """ G = 0.8*(Rz_Ui+Rz_Iu)*1E-3 """
        G = 0.8 * (self.Rz_UI + self.Rz_Iu)
        return G
    
    def E12_15_S_nmin(self):
        """ S_nmin = Z_min + G
        S_nmin:[mm] """
        S_nmin = self.Z_min + self.G
        return S_nmin
    
    def E12_16A_p_Fmax_Outer(self):
        """ p_Fmax = (ReU/SpU)*(1-Q_u**2/(sqrt(3))) """
        p_Fmax = (self.R_eU/self.S_pU) * ((1-self.Q_U**2)/sp.sqrt(3))
        return p_Fmax
    
    def E12_16B_p_Fmax_InnerHollow(self):
        """ p_Fmax = (ReI/SpI)*(1-Q_I**2/(sqrt(3))) """
        p_Fmax = (self.R_eI/self.S_pI) * ((1-self.Q_I**2)/sp.sqrt(3))
        return p_Fmax
        
    def E12_16B_p_Fmax_InnerMassive(self):
        """ p_Fmax = (ReI/SpI)*(2/(sqrt(3))) """
        p_Fmax = (self.R_eI/self.S_pI) * ((2)/sp.sqrt(3))
        return p_Fmax    
    
    def E12_17_Z_max(self):
        """ Z_max = (p_Fmax * D_F / E_U ) * K """
        Z_max = (self.p_Fmax * self.D_F / self.E_U) * self.K
        return Z_max
    
    def E12_18_Snmax(self):
        """ S_nmax = Z_max + G """
        S_nmax = self.Z_max + self.G
        return S_nmax
    
    def E12_19_P_T(self):
        """ P_T = (S_nmax - S_nmin) *1E3 """
        P_T = (self.S_nmax - self.S_nmin)
        return P_T
    
    def E12_20_P_T(self):
        """ P_T = T_B + T_A """
        P_T = self.T_B + self.T_A
        return P_T

    def E12_20_hA_TB(self):
        """ T_B = 0.6*P_T 
        for groups 1-3 0.6 is adviced"""
        T_B = 0.6  * self.P_T
        return T_B
    
    def E12_21_F_i_PressFittingForce(self):
        """ F_i = A_F * pprime_Fmax * mu_i """
        F_i = self.A_F * self.pprime_Fmax * self.mu_i
        return F_i
        
    def E12_22_T_U_TemperatureOuterPart(self):
        """ T_U = T_room + (Sprime_nmax+Sm)/(alpha_U.D_F) + alpha_I/alpha_U * (T_I - T_room) """
        T_U = self.T_room + (self.Sprime_nmax+self.S_m)/(self.alpha_U*self.D_F) + self.alpha_I/self.alpha_U * (self.T_I -self.T_room)
        return T_U
    
    def E12_23a_n_g_BoundarySpeed(self):
        """ n_g = 60* 2 / (pi *D_Uu) * sqrt(2*pprime_Fmin / (((3 + nu_U)*(1-Q_U**2)*rho))) 
            [n_g] in rpm , formula 12_23a * 60 seconds
            [rho] in kg/mm³  
            [D_Uu] in mm
        """
        n_g = 60 *1E6 * 2 / (sp.pi *  self.D_Uu) * sp.sqrt( 2*self.pprime_Fmin / ((3+self.nu_U)*(1-self.Q_U)*self.rho_U)) 
        return n_g
    
    def E12_23b_n_g_BoundarySpeed(self):
        """ n_g = 29.7*1E6 * sqrt(2*pprime_Fmin / (((D_Uu**2)*(1-Q_U**2)*rho))) 
            [n_g] in rpm
            [rho] in kg/mm³  
            [D_Uu] in mm
            [pprime_Fmin] in N/mm²
        """
        n_g = 29.7*1E6 * sp.sqrt( 2*self.pprime_Fmin / ((self.D_Uu**2)*(1-self.Q_U)*self.rho_U)) 
        return n_g

    def E12_24_T_n_TransmissbleTorque(self):
        """ T_n = T * (1 - (N/n_g)**2) """
        T_n = self.T * (1 - (self.N / self.n_g)**2)
        return T_n

    def E12_25_C_ConeRatio(self):
        """ C = (D_1 - D_2)/l """
        C = (self.D_1-self.D_2)/self.l 
        return C
        
    def E12_26_alpha_ConicalAngle(self):
        """ tan(alpha/2) = (D_1-D_2)/(2*l) """
        alpha = sp.atan((self.D_1-self.D_2)/(2*self.l)) / sp.pi *180
        return alpha
    
    def E12_27A1_a_min_MinPressingDist(self):
        """ a_min = (S_nmin/2) / (tan(alpha/2)) """
        a_min = (self.S_nmin/2) / sp.tan(self.alpha/2 /180*sp.pi)
        return a_min
    
    def E12_27A2_a_min_MinPressingDist(self):
        """ a_min = (Z_min + G)/2 / (tan(alpha/2)) """
        a_min = ((self.Z_min + self.G) /2) / sp.tan(self.alpha/2 /180*sp.pi)
        return a_min
    
    def E12_27B1_a_max_MaxPressingDist(self):
        """ a_min = (S_nmax/2) / (tan(alpha/2)) """
        a_min = (self.S_nmax/2) / sp.tan(self.alpha/2 /180*sp.pi)
        return a_min
    
    def E12_27B2_a_max_MaxPressingDist(self):
        """ a_min = (Z_max + G)/2 / (tan(alpha/2)) """
        a_min = ((self.Z_max + self.G) /2) / sp.tan(self.alpha/2 /180*sp.pi)
        return a_min
    
    def E12_28A_Z_min_MinPressFit(self):
        """ Z_min = (p_Fmin * D_Fgem*K)/ (E_U*cos(alpha/2)) """
        Z_min = (self.p_Fmin * self.D_Fgem * self.K ) / (self.E_U - sp.cos(self.alpha/2 /180*sp.pi))
        return Z_min
    
    def E12_28B_Z_max_MaxPressFit(self):
        """ Z_max = (p_Fmax * D_Fgem*K)/ (E_U*cos(alpha/2)) """
        Z_max = (self.p_Fmax * self.D_Fgem * self.K ) / (self.E_U - sp.cos(self.alpha/2 /180*sp.pi))
        return Z_max
    
    def E12_28_A_D_Fgem_AvergeDiameter(self):
        """ D_Fgem = (D_1+D_2)/2 """
        D_Fgem = (self.D_1 + self.D_2)/2
        return D_Fgem
    
    def E12_29_F_i_PressFitForce(self):
        """ F_i = (2*S_S*T*sin(rho_i+alpha/2)) / (D_Fgem*mu*cos(rho_i)) """
        F_i = (2*self.S_S * self.T * sp.sin((self.rho_i+self.alpha/2)/180*sp.pi)) / (self.D_Fgem*self.mu*self.rho_i)
        return F_i
    
    def E12_29_A_rho_i_PressFitFrictionAngle(self):
        """ tan(rho_i) = mu_i """
        rho_i = sp.atan(self.mu_i) /sp.pi*180
        return rho_i
    
    def E12_30_p_F_SurfacePressure(self):
        """ p_F = (F_i * cos(rho_i)*cos(alpha/2)) / (D_Fgem * pi * l * sin (rho_i+alpha/2)) """
        p_F = (self.F_i * sp.cos(self.rho_i/180*sp.pi) * sp.cos(self.alpha/2/180*sp.pi)) / (self.D_Fgem * sp.pi * self.l * sp.sin((self.rho_i+self.alpha/2)/180*sp.pi)) 
        return p_F
    
    def E12_31_p_Fmin_MinSurfacePressure(self):
        """ p_Fmin = (2 * S_S * T * cos(alpha/2)) / (D_Fgem**2 * pi * mu * l) """
        p_Fmin = (2 * self.S_S * self.T * sp.cos(self.alpha/2/180*sp.pi)) / (sp.pi * self.mu * self.l * self.D_Fgem**2)
        return p_Fmin
    
    def E12_32_T_MaxTransmissbleTorque(self):
        """ T = (Z_min * E_U * D_Fgem * pi * mu * l) / (2 * S_S * K) """
        T = (self.Z_min * self.E_U *self.D_Fgem * sp.pi * self.mu * self.l) / (2 * self.S_S * self.K)
        return T
    
    def E12_33A_T_tot_TransmissbleTorque(self):
        """ T_tot = T_Tab*f_n """
        T_tot = self.T_Tab * self.f_n
        return T_tot
    
    def E12_33b_F_atot_TotalAxialForce(self):
        """ F_atot = F_aTab * f_n """
        F_atot = self.F_aTab * self.f_n
        return F_atot
    
    def E12_34_T_res_ResultungTorque(self):
        """ T_res = sqrt(T**2 + (F_a*D_F/2)**2) """
        T_res = sp.sqrt(self.T**2 + (self.F_a * self.D_F/2)**2)
        return T_res
    
    def E12_36A_D_Uu_HubOuterDiamter(self):
        """ D_Uu = D*sqrt((R_eU+p_N*C)/(R_eU-pN+C))+d """
        D_Uu = self.D * sp.sqrt((self.R_eU+self.p_N*self.C) /(self.R_eU-self.p_N*self.C) ) + self.d
        return D_Uu
    
    def E12_36B_DIi_HollowShaftInnerDiameter(self):
        """ D_Ii = D_F*sqrt((R_eI-2*P_A*C)/R_eI)-d """
        D_Ii = self.D_F * sp.sqrt((self.R_eI-2*self.p_A*self.C)/self.R_eI)
        return D_Ii
    
    def E12_37_p_Fmin_MinContactPressureClamp(self):
        """ p_Fmin = (2*K_A*T_nom*S_S) / (pi*l_F*mu*D_f**2) * K"""
        p_Fmin = (2 * self.K_A * self.T_nom * self.S_S) / (sp.pi * self.l_F * self.mu * self.D_F**2) * self.K
        return p_Fmin
    
    def E12_38_F_Kl_CalmpingForceBolt(self):
        """ F_Kl = (2*K_A*T_nom*S_S*K) / (n*pi*D_F*mu) """
        F_Kl = (2 * self.K_A * self.T_nom * self.S_S * self.K ) /  (self.n * sp.pi * self.D_F * self.mu)
        return F_Kl
    
    def E12_39_p_F_EffectiveClampingForce(self):
        """ p_F  = n * F_VM / (D_F * l_F) """
        p_F = (self.n * self.F_VM) / (self.D_F * self.l_F)
        return p_F
    
    def E12_40_F_N_RequiredClampingForce(self):
        """ F_N = (K_A * T_nom) / (D_F * mu) """
        F_N = (self.K_A * self.T_nom) / (self.D_F * self.mu)
        return F_N
    
    def E12_41_F_Kl_RequiredClampingForceBolt(self):
        """ F_Kl = (K_A * T_nom * S_S * l_1) / (n * D_F * mu * l2) """
        F_Kl = (self.K_A * self.T_nom * self.S_S * self.l_1) / (self.n * self.D_F * self.mu * self.l_2)
        return F_Kl
    
    def E12_42_p_F_EffectiveClampingForce(self):
        """ p_F = (n * F_VM * l_2) / (D_F * l_F *l_1) """
        p_F = (self.n * self.F_VM * self.l_2) / (self.D_F * self.l_F * self.l_1)
        return p_F
    
    
        
    