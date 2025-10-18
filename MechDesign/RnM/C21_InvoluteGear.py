# -*- coding: utf-8 -*-
"""
Some variables are written in the book as x_1,2 but hence cannot be denoted in that way in the code so are referred to as x_12
a variable which in the book is written as x__* is written in the code as x__star


This folder will contain - for student's reference - expressions from Roloff & Matek. 
These expressions are only for internal use and e.g. for studying/answering the exam. 
These expressions may never and under no exception be distributed or shared to anyone else 
and serve only for the purpose of students that follow the KULeuven courses on dimensioning of machine parts 
By using these files you commit to adhering to this condition ! 
"""

import sympy as sp
from . import MySymbol
# from ..Units.Units import mm_



class InvoluteGear:
    
 
    MySymbolDict = {'d'             : '[mm] pitch circle diameter' ,\
                    'F'             : '[N]' ,\
                    'epsilon_gama'  : '[] paasage ratio' ,\
                    'sigma_H'       : '[]' ,\
                    'y_alpha'       : '[]' ,\
                    'p'             : '[rad] pitch' ,\
                    'm'             : '[mm] modulus which is pitch dependent' ,\
                    'F_t'           : '[N] nominal circumferential force' ,\
                    'f_HBeta'       : '[]' ,\
                    'alpha'         : '[rad] pressure angle (profile angle) according to DIN 867' ,\
                    'alpha_sp'      : '[rad]' ,\
                    'alpha_n'       : '[rad]' ,\
                    'z'             : '[] number of teeth' ,\
                    'd_b'           : '[mm] circle base diameter' ,\
                    'p_b'           : '[mm] perimeter base pitch' ,\
                    'h_a'           : '[] roll of circle' ,\
                    'h_ap'          : '[]' ,\
                    'd_a'           : '[]' ,\
                    'h_v'           : '[]' ,\
                    'h_vp'          : '[]' ,\
                    'c'             : '[mm] top clearance for gear according to DIN867' ,\
                    'd_v'           : '[]' ,\
                    'd_12'          : '[mm] pitch circle diameters 𝑑1 and 𝑑2 (index 1 for driver 1, index 2 for follower 2)' ,\
                    'z_12'          : '[]' ,\
                    'd_1'           : '[mm] pitch circle diameter of small wheel' ,\
                    'd_2'           : '[mm] pitch circle diameter of big wheel' ,\
                    'z_1'           : '[] # of tooth of small wheel' ,\
                    'z_2'           : '[] # of tooth of big wheel' ,\
                    'omega_1'       : '[rad/s] angular velocity' ,\
                    'omega_2'       : '[rad/s] angular velocity' ,\
                    'n_1'           : '[]' ,\
                    'n_2'           : '[]' ,\
                    'z_bigWheel'    : '[] # of tooth of big wheel' ,\
                    'z_smallWheel'  : '[] # of tooth of small wheel' ,\
                    'd_a1'          : '[mm] top circle diameter according to equation (21.6)' ,\
                    'a_d'           : '[mm] zero-profile shift distance according to equation (21.8)' ,\
                    'u'             : '[] tooth ratio z_2/z_1' ,\
                    'd_a2'          : '[mm] top circle diameter according to equation (21.6)' ,\
                    'd_b1'          : '[mm] base circle diameter according to equation (21.2)' ,\
                    'd_b2'          : '[mm] base circle diameter according to equation (21.2)' ,\
                    'g_a'           : '[] engagement path 21.12' ,\
                    'p_i'           : '[] engagement pitch' ,\
                    'epsilon_alpha' : '[] engagement quotient 21.13' ,\
                    'x'             : '[] The profile shift factor 𝑥 (as part of the modulus)' ,\
                    'V'             : '[mm] distance shifted from pitch of circle 21.15' ,\
                    'z__1_g'        : '[] practical limit number of teeth' ,\
                    'z_g'           : '[] limit number or minimum number of teeth' ,\
                    'd_r1'          : '[mm] operating rolling circles diameter' ,\
                    'd_r2'          : '[mm] operating rolling circles diameter' ,\
                    'alpha_r'       : '[rad] roller pressure angle' ,\
                    'a'             : '[mm] specific center distance' ,\
                    'k'             : '[] addendum factor' ,\
                    'k__star'       : '[] 𝑦 − 𝛴𝑥 Addendum change factor, when y is the pitch circle distance factor from (21.20)' ,\
                    'x_1'           : '[] Profile shift factor of pinion 1' ,\
                    'x_2'           : '[] Profile shift factor of wheel 2' ,\
                    'V_1'           : '[mm] distance shifted from pitch of circle 21.15' ,\
                    'V_2'           : '[mm] distance shifted from pitch of circle 21.15' ,\
                    'd_y'           : '[rad] halve tooth thickness angle 𝜓 (see also figure 21-1)' ,\
                    'a_a'           : '[]' ,\
                    'a_sp'          : '[]' ,\
                    'e_y'           : '[mm] hole width' ,\
                    'e'             : '[]' ,\
                    'Nu'            : '[] Half tooth space angle' ,\
                    'SigmaA_sti'    : '[] minimum theoretical radial play follows from the tooth thickness deviations ( transverse cross-section)' ,\
                    'Deltaj_ae'     : '[]' ,\
                    'SigmaA_sni'    : '[] minimum theoretical radial play follows from the tooth thickness deviations ( normal cross-section)' ,\
                    'beta'          : '[rad]' ,\
                    'SigmaA_ste'    : '[] maximum theoretical radial play follows from the tooth thickness deviations ( transverse cross-section)' ,\
                    'Deltaj_ai'     : '[]' ,\
                    'SigmaA_sne'    : '[] maximum theoretical radial play follows from the tooth thickness deviations ( normal cross-section)' ,\
                    'j_tMax'        : '[] Maximal theoretical radial play' ,\
                    'j_tMin'        : '[] Minimum theoretical radial play' ,\
                    'A_a'           : '[]' ,\
                    'a_n'           : '[]' ,\
                    'z_n'           : '[]virtual number of teeth (number of teeth of the replacement cylindrical wheel, with straight teeth = 𝑧)' ,\
                    'W_k'           : '[] Nominal size for tooth width' ,\
                    'j_a'           : '[]' ,\
                    'd_sp'          : '[] teeth diameter when becoming pointed' ,\
                    's_a'           : '[] tooth thickness at top of circle' ,\
                    's_y'           : '[] halve tooth thickness angle 𝜓 (see also figure 21-1)' ,\
                    'epsilon_a'     : '[] engagment quotient' ,\
                    'x_limit'       : '[] boundary value for x' ,\
                    'i'             : '[] gear ratio' ,\
                    'd_v12'         : '[] pitch circl diameters for driver(1) and follower(2)' ,\
                    'i_1'           : '[] gear ratio of wheel 1' ,\
                    'i_2'           : '[] gear ratio of wheel 2' ,\
                    'i_n'           : '[] gear ration of wheel n' ,\
                    'u_1'           : '[] Ratio number of teeth of wheel 1' ,\
                    'u_2'           : '[] Ratio number of teeth of wheel 1' ,\
                    'u_n'           : '[] ratio number of teeth of wheel n' ,\
                    'm__1_n'        : '[] modulus for pinion on shaft or pinion-shaft' ,\
                    'd_shaft'       : '[mm] diameter of shaft' ,\
                    'm__2_n'        : '[] modulus for shaft distance a' ,\
                    'm__3_n'        : '[] modulus for shaft distance a' ,\
                    'psi_d'         : '[] width-diameter ratio according to table 21-14a' ,\
                    'sigma_FVLim1'  : '[] resistance to tooth fracture of the pinion tooth base according to tables 20-1 and 20-2' ,\
                    'T_1eq'         : '[Nm] torque to be transmitted through the gear pair' ,\
                    'T_1'           : '[NM]' ,\
                    'sigma_HLim'    : '[] resistance to flank damage of the softest material according to tables 20-1 and 20-2' ,\
                    'T_1Nom'        : '[Nm] nominal torque' ,\
                    'F_t1'          : '[N] nominal circumferential force' ,\
                    'T_2Nom'        : '[Nm] Nominal torque' ,\
                    'F_t2'          : '[N] nominal circumferential force' ,\
                    'F_t12'         : '[N] nominal circumferential force' ,\
                    'F_bn12'        : '[N] Tooth normal forces' ,\
                    'T_12'          : '[Nm] torque for circles 1 and 2' ,\
                    'd_r12'         : '[mm] diameters for roller circles 1 and 2' ,\
                    'F_r12'         : '[N] radial forces' ,\
                    'F_a12'         : '[N] axial forces' ,\
                    'K_v'           : '[] Dynamic factor' ,\
                    'K_1'           : '[] factors according to table 21-15' ,\
                    'K_A'           : '[] Operating factor table 3-5' ,\
                    'b'             : '[mm] tooth width, with unequal wheel width the smallest tooth width' ,\
                    'K_2'           : '[] factors according to table 21-15' ,\
                    'K_3'           : '[] = 0,01 ⋅ 𝑧1⋅ 𝑣𝑡⋅ √𝑢2 ∕ (1 +𝑢2) ≤ 10 m/s' ,\
                    'K_HBeta'       : '[] Load concentration factors table 21-18' ,\
                    'K_VBeta'       : '[] Load concentration factors table 21-18' ,\
                    'F_betaY'       : '[] effective flank line deviation tbable 21-17' ,\
                    'F_m'           : '[N] Average load' ,\
                    'N_f'           : '[]' ,\
                    'f_as'          : '[] flank line deviation due to deformation table 21-16a or equation 21.75' ,\
                    'K__1'          : '[] factor that takes into account the pinion position relative to the bearings, depending on 𝑠 and 𝑙; values according to table 21-16b' ,\
                    'l'             : '[]' ,\
                    's'             : '[]' ,\
                    'd_as'          : '[mm] pinion shaft diameter' ,\
                    'f_ma'          : '[micro meter] manufacturing-dependent flank deviation' ,\
                    'q_H'           : '[1]' ,\
                    'F_betaX'       : '[] Deviation from the flank line' ,\
                    'F_betaXMin'    : '[] Minimum deviation from the flank line' ,\
                    'y_beta'        : '[] run in value table 21-17' ,\
                    'N_v'           : '[] Factor for tooth base' ,\
                    'h'             : '[mm] tooth height' ,\
                    'K_HAlpha'      : '[] Perimeter factor tab;e 21-19' ,\
                    'K_VAlpha'      : '[] Perimeter factor tab;e 21-19' ,\
                    'c_gama'        : '[] Engagement rigidity' ,\
                    'f_pe'          : '[] Maximum value for engagement pitch deviation table 21-19b' ,\
                    'y_a'           : '[] inloopwaarde; value according to table 21-19c' ,\
                    'F_tH'          : '[] relevant circumferential force, 𝐹𝑡𝐻 = 𝐹𝑡⋅ 𝐾𝐴 ⋅ 𝐾𝐻⋅ 𝐾' ,\
                    'K_VTot'        : '[] resulting load influence according to equation (21.81)' ,\
                    'K_V'           : '[]' ,\
                    'K_VX'          : '[]' ,\
                    'M'             : '[]' ,\
                    'W'             : '[]' ,\
                    'sigma_b'       : '[] Nominal bending stress' ,\
                    'alpha_Fin'     : '[]' ,\
                    'H_Fin'         : '[]' ,\
                    's_Vn'          : '[]' ,\
                    'alpha_w'       : '[]' ,\
                    'Y_Va'          : '[] sahpe factor table 21-20a' ,\
                    'sigma_V0'      : '[] local tooth root stress according to equation (21.82)' ,\
                    'm_n'           : '[] resulting load influence according to equation (21.81)' ,\
                    'Y_Ka'          : '[] Notch factor table 21-20a' ,\
                    'Y_epsilon'     : '[] Engagment factor table 21-20b' ,\
                    'Y_beta'        : '[] Tooth angle factor table 21-20c' ,\
                    'sigma_V012'    : '[] local tooth root stress according to equation (21.82' ,\
                    'sigma_VG'      : '[] permissible tooth-base stress' ,\
                    'sigma_VLim'    : '[] nominal tooth base fatigue strength against bending of the test wheels according to table 20-1 and table 20-2' ,\
                    'Y_ST'          : '[] From notch factor for test wheels' ,\
                    'Y_NT'          : '[] Life span factor table 21-21a' ,\
                    'Y_deltaRelT'   : '[] Relative support factor from table 21-21b' ,\
                    'Y_RRelT'       : '[] Relative roughness factor from table 21-21d' ,\
                    'Y_X'           : '[] Size factor table 21-21d' ,\
                    'S_V12'         : '[] safety of the tooth-base strength for pinion and wheel' ,\
                    'sigma_VG12'    : '[]permissible tooth-base stress for pinion and wheel' ,\
                    'sigma_V12'     : '[] nominal tooth base fatigue strength for pinion and wheel' ,\
                    'S_VMin'        : '[] minimum safety factor for the foot load.' ,\
                    'sigma_HMax'    : '[] max stress according to Hertz result' ,\
                    'sigma_max'     : '[] max stress' ,\
                    'rho'           : '[] reduced radius of curveture' ,\
                    'E'             : '[] Young,\'s Modulus' ,\
                    'nu'            : '[] Poisson coefficient' ,\
                    'alpha_y2'      : '[] Pressure angle from 𝑐𝑜𝑠 𝛼𝑦 = 𝑑 ⋅ 𝑐𝑜𝑠 𝛼 ∕ 𝑑�' ,\
                    'alpha_y1'      : '[] Pressure angle from 𝑐𝑜𝑠 𝛼𝑦 = 𝑑 ⋅ 𝑐𝑜𝑠 𝛼 ∕ 𝑑�' ,\
                    'rho_y'         : '[] radius of curveture' ,\
                    'sigma_HY'      : '[] Egagement position Y' ,\
                    'sigma_HC'      : '[] Stress in pole C' ,\
                    'Z_H'           : '[] Geometry facotr 𝑍𝐻 = 2 ⋅ 𝑐𝑜𝑠 𝛽𝑏 ∕ 𝑐𝑜𝑠2 𝛼𝑡⋅ 𝑡𝑎𝑛 𝛼𝑟𝑡,' ,\
                    'Z_E'           : '[] Elasticity factor √0,175 ⋅ 𝐸., 𝐸 = (2 ⋅ 𝐸1 ⋅ 𝐸2) ∕ (𝐸1 + 𝐸2)' ,\
                    'sigma_H0'      : '[] Contact stress' ,\
                    'Z_epsilon'     : '[] Engagement factor values from table 21-22c' ,\
                    'Z_beta'        : '[] Tooth angle factor 𝑍𝛽 = √cos' ,\
                    'K_HTot'        : '[] Load influence factor equation 21.81' ,\
                    'Z_NT'          : '[] Life span factor table 21-23d' ,\
                    'Z_L'           : '[] Viscosity factor table 21-23a' ,\
                    'Z_vir'         : '[] Sliding factor table 21-23b' ,\
                    'Z_R'           : '[] Roughness factor table 21-23c' ,\
                    'Z_M'           : '[] Material Combination factor table 21-23e' ,\
                    'Z_X'           : '[] Size factor table 21-21d' ,\
                    'sigma_HG12'    : '[] Allowable flank stress for circle 1 and 2' ,\
                    'S_HMin'        : '[] required minimum safety for contact strength' ,\
                    'S_H12'         : '[] Safety flank strength' ,\
                     }

    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(kkk) for kkk in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))
            
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------





    def inv(self, angle):
        """ inv(alpha) = tan(alpha) - alpha"""
        p_gem = sp.tan(angle) - angle
        return p_gem

    def E21_1A_PitchCircleDiameter(self):
        """ d = 2p/pi """
        p_gem = 2*self.p/sp.pi
        return p_gem

    def E21_1B_PitchCircleDiameter(self):
        """ d = z* m"""
        p_gem = self.z * self.m
        return p_gem

    def E21_2A_BaseCircleDiameter(self):
        """ d_b = d * cos(alpha) """
        p_gem = self.d * sp.cos(self.alpha)
        return p_gem

    def E21_2B_BaseCircleDiameter(self):
        """ d_b = z * m * cos(alpha)"""
        p_gem = self.z * self.m * sp.cos(self.alpha)
        return p_gem

    def E21_3A_PerimeterBasePitch(self):
        """ p_b = d_b * pi / z"""
        p_gem = self.d_b * sp.pi / self.z
        return p_gem

    def E21_3B_PerimeterBasePitch(self):
        """ p_b = p * cos(alpha)"""
        p_gem = self.p * sp.cos(self.alpha)
        return p_gem

    def E21_4A_PitchOnLineOfContact(self):
        """ P_i = P_b = p * cos(alpha)"""
        p_gem = self.p * sp.cos(self.alpha)
        return p_gem

    def E21_4b_PitchOnLineOfContact(self):
        """ P_i = pi * m * cos(alpha)"""
        p_gem = sp.pi * self.m * sp.cos(self.alpha)
        return p_gem

    def E21_5A_ToothDimension(self):
        """ h_a = h_ap = m"""
        return self.m

    def E21_5B_ToothDimensiosn(self):
        """ h_a = 0.5(d_a - d)"""
        p_gem = 0.5*(self.d_a - self.d)
        return p_gem

    def E21_5CToothDimension(self):
        """ h_v = h_vp = m + c"""
        p_gem = self.m + self.c
        return p_gem
    def E21_5D_ToothDimension(self):
        """ h_v = 0.5(d - d_v_"""
        p_gem = 0.5*(self.d - self.d_v)
        return p_gem

    def E21_5EToothDimension(self):
        """ h = h_a + h_v = 2m+c"""
        p_gem = 2*self.m + self.c
        return p_gem

    def E21_5F_ToothDimension(self):
        """ h = 0.5(d_a - d_v_"""
        p_gem = 0.5*(self.d_a - self.d_v)
        return p_gem

    def E21_6A_PitchCircleDiameter(self):
        """ d_a12 = d_12 + 2*h_a"""
        p_gem = self.d_12 + 2 * self.h_a
        return p_gem

    def E21_6B_PitchCircleDiameter(self):
        """ d_a12 = m(z_12 + 2)"""
        p_gem = self.m*(self.z_12 + 2)
        return p_gem

    def E21_7A_DedendumCircleDiameters(self):
        """ d_v12 = d_12 + 2"""
        p_gem = self.d_12 + 2
        return p_gem

    def E21_7B_PitchCircleDiameter(self):
        """ d_v12 = m(z_12 - 2.5)"""
        p_gem = self.m*(self.z_12 - 2.5)
        return p_gem

    def E21_8A_ZeroProfileShiftCenterDistance(self):
        """ a_d = (d_1 + d_2)/2"""
        p_gem = (self.d_1 + self.d_2) / 2
        return p_gem

    def E21_8B_ZeroProfileShiftCenterDistance(self):
        """ a_d = m/2 * (z_1 + z_2)"""
        p_gem = (self.m / 2) * (self.z_1 + self.z_2)
        return p_gem

    def E21_9A_Transmission(self): #Not clear about the name
        """i = omega_1/omega_2"""
        p_gem = self.omega_1 / self.omega_2
        return p_gem

    def E21_9B_Transmission(self):
        """i = n_1/n_2"""
        p_gem = self.n_1 / self.n_2
        return p_gem

    def E21_9C_Transmission(self):
        """i = d_2/d_1"""
        p_gem = self.d_2 / self.d_1
        return p_gem

    def E21_9D_Transmission(self):
        """i = z_2/z_1"""
        p_gem = self.z_2 / self.z_1
        return p_gem
    
    def E21_10A_ToothRation(self):
        """ u = z_bigWheel / z_smallWheel"""
        p_gem = self.z_bigWheel / self.z_smallWheel
        return p_gem

    def E21_10B_ToothRatio(self):
        """ u = z_2 / z_1"""
        p_gem = self.z_2 / self.z_1
        return p_gem

    def E21_11A_PitchCircleDiameterPinion1(self):
        """ d_1 = m*z_1"""
        p_gem = self.m * self.z_1
        return p_gem

    def E21_11B_PitchCircleDiameterPinion1(self):
        """ d_1 = d_a1 - 2*m"""
        p_gem = self.d_a1 - 2* self.m
        return p_gem

    def E21_11C_PitchCircleDiameterPinion1(self):
        """ d_1 = (z_1 * d_a1) / (z_1 +2)"""
        p_gem = (self.z_1 * self.d_a1) / (self.z_1 + 2)
        return p_gem

    def E21_11D_PitchCircleDiameterPinion1(self):
        """ d_1 = 2a_d / (1+u)"""
        p_gem = 2 * self.a_d / (1 + self.u)
        return p_gem

    def E21_11E_PitchCircleDiameterWheel2(self):
        """ d_2 = m*z_2"""
        p_gem = self.m * self.z_2
        return p_gem

    def E21_11F_PitchCircleDiameterWheel2(self):
        """ d_2 = d_a2 - 2*m"""
        p_gem = self.d_a2 - 2* self.m
        return p_gem

    def E21_11G_PitchCircleDiameterWheel2(self):
        """ d_2 = (z_2 * d_a2) / (z_2 +2)"""
        p_gem = (self.z_2 * self.d_a2) / (self.z_2 + 2)
        return p_gem

    def E21_11H_PitchCircleDiameterWheel2(self):
        """ d_1 = 2* a_d * o / (1+u) with o = i"""
        p_gem = 2 * self.a_d * self.i / (1 + self.u)
        return p_gem

    def E21_12_EngagementPath(self):
        """ g_a = 0.5 * (sqrt(d_a1**2-d_b1**2) + z_2/|z_2| * sqrt(d_a2**2 - d_b2**2) - a_d*sin(alpha)"""
        p_gem = 0.5 * (sp.sqrt(self.d_a1**2 - self.d_b1**2) + (self.z_2 / abs(self.z_2)) * sp.sqrt(self.d_a2**2 - self.d_b2**2)) - self.a_d * sp.sin(self.alpha)
        return p_gem

    def E21_13A_EngagementQuotient(self):
        """ epsilon_a = g_a / p_i"""
        p_gem = self.g_a / self.p_i
        return p_gem

    def E_21_13B_EngagementQuotient(self):
        """ epsilon_a = (0.5 * (sqrt(d_a1**2-d_b1**2) + z_2/|z_2| * sqrt(d_a2**2 - d_b2**2) - a_d*sin(alpha)) / (pi * m * cos(alpha)"""
        p_gem = (0.5 * (sp.sqrt(self.d_a1 ** 2 - self.d_b1 ** 2) + (self.z_2 / abs(self.z_2)) * sp.sqrt(
            self.d_a2 ** 2 - self.d_b2 ** 2)) - self.a_d * sp.sin(self.alpha)) / (sp.pi * self.m * sp.cos(self.alpha))
        return p_gem

    def E21_14_TheoreticalNumberOfTeeth(self):
        """ z_g 2 / sin(alpha)**2"""
        p_gem = 2 / (sp.sin(self.alpha) * sp.sin(self.alpha))
        return p_gem

    def E21_15_ProfileShift(self):
        """ V = x*m"""
        p_gem = self.x * self.m
        return p_gem

    def E21_16A_ProfileShiftFactor(self):
        """ x_limit = (z__1_g - z)/z_g"""
        p_gem = (self.z__1_g - self.z) / self.z_g
        return p_gem

    def E21_16B_ProfileShiftFactor(self):
        """ x_limit = (14 - z)/17"""
        p_gem = (14 - self.z) / 17
        return p_gem

    def E21_17A_NominalSizeOfToothThickness(self):
        """ s = p/2 + 2*V*taN(alpha)"""
        p_gem = self.p/2 + 2 * self.V * sp.tan(self.alpha)
        return p_gem

    def E21_17B_NominalSizeOfToothThickness(self):
        """ s = m*(pi/2 + 2 * x * tan(alpha)"""
        p_gem = self.m*(sp.pi / 2 + 2 * self.x * sp.tan(self.alpha))
        return p_gem

    def E21_18A_ToothSpace(self):
        """ e = p/2 + 2*V*taN(alpha)"""
        p_gem = self.p / 2 + 2 * self.V * sp.tan(self.alpha)
        return p_gem

    def E21_18B_ToothSpace(self):
        """ e = m*(pi/2 + 2 * x * tan(alpha)"""
        p_gem = self.m*(sp.pi / 2 + 2 * self.x * sp.tan(self.alpha))
        return p_gem

    def E21_19A_RollerPressureAngle(self):
        """ a = ((d_r1 + d_r2) / 2)"""
        p_gem = (self.d_r1 + self.d_r2) /2
        return p_gem

    def E21_19B_RollerPressureAngle(self):
        """ a = ((d_1 + d_2) / 2) * cos(alpha)/cos(alpha_r)"""
        p_gem = ((self.d_1 + self.d_2) / 2) * sp.cos(self.alpha) / sp.cos(self.alpha_r)
        return p_gem

    def E21_19C_RollerPressureAngle(self):
        """ a = a_d * cos(alpha)/cos(alpha_r)"""
        p_gem = self.a_d * sp.cos(self.alpha) / sp.cos(self.alpha_r)
        return p_gem

    def E21_20_ShaftShiftFactor(self):
        """ y = 0.5 * (z_1 + z_2) * [cos(alpha)/cos(alpha_r) -1]"""
        p_gem = 0.5 * (self.z_1 + self.z_2) * (sp.cos(self.alpha) / sp.cos(self.alpha_r) -1)
        return p_gem

    def E21_21_RollPressureAngle(self):
        """ alpha_r = arccos(a_d * cos(alpha) / a"""
        p_gem = sp.acos(self.a_d * sp.cos(self.alpha) / self.a)
        return p_gem

    def E21_22AA_OperatingRollingCircleDiameter(self):
        """ d_r1 = d_1 * cos(alpha) / cos(alpha_r)"""
        p_gem = self.d_1 * sp.cos(self.alpha) / sp.cos(self.alpha_r)
        return p_gem

    def E21_22AB_OperatingRollingCircleDiameter(self):
        """ d_r1 = 2*a / (1+u)"""
        p_gem = 2 * self.a / (1 + self.u)
        return p_gem

    def E21_22AC_OperatingRollingCircleDiameter(self):
        """ d_r1 = 2* z_1 * a / (z_1 + z_2)"""
        p_gem = 2 * self.z_1 * self.a / (self.z_1 + self.z_2)
        return p_gem

    def E21_22BA_OperatingRollingCircleDiameter(self):
        """ d_r2 = d_2 * cos(alpha) / cos(alpha_r)"""
        p_gem = self.d_2 * sp.cos(self.alpha) / sp.cos(self.alpha_r)
        return p_gem

    def E21_22BB_OperatingRollingCircleDiameter(self):
        """ d_r2 = 2*a - d_r1"""
        p_gem = 2 * self.a - self.d_r1
        return p_gem

    def E21_22BC_OperatingRollingCircleDiameter(self):
        """ d_r2 = 2*a*u / (1+u)"""
        p_gem = 2 * self.a * self.u / (1 + self.u)
        return p_gem

    def E21_22BD_OperatingRollingCircleDiameter(self):
        """ d_r2 = 2* z_2 * a / (z_1 + z_2)"""
        p_gem = 2 * self.z_2 * self.a / (self.z_1 + self.z_2)
        return p_gem

    def E21_23A_AddendumFactor(self):
        """ k = k__star * m"""
        p_gem = self.k__star * self.m
        return p_gem

    def E21_23B_AddendumFactor(self):
        """ k = a - a_d - m*(x_1 + x_2)"""
        p_gem = self.a - self.a_d - self.m * (self.x_1 + self.x_2)
        return p_gem

    def E21_24A_TopCircleDiameter(self):
        """ d_a1 = d_1 + 2* m + 2*V_1 + 2*k"""
        p_gem = self.d_1 + 2 * self.m + 2 * self.V_1 + 2 * self.k
        return p_gem

    def E21_24B_TopCircleDiameter(self):
        """ d_a2 = d_2 + 2(m + V_2 = K)"""
        p_gem = self.d_2 + 2 * (self.m + self.V_2 + self.k)
        return p_gem

    def E21_25A_DedendumCicleDiameter(self):
        """ d_v = d - 2 * h_v + 2* nu"""
        p_gem = self.d - 2 * self.h_v + 2 * self.nu
        return p_gem

    def E21_25B_DedendumCicleDiameter(self):
        """ d_v = d - 2[m + c - V]"""
        p_gem = self.d - 2 * (self.m + self.c - self.V)
        return p_gem

    def E21_26_EngagementQuotient(self):
        """ epsilon_a = (0.5 * (sqrt(d_a1**2 - d_b1**2) + z_2 / abs(z_2) * sqrt(d_a2**2 - d_b2**2)) - a * sin(alpha_r) / (pi * m * cos(alpha)) """
        p_gem = (0.5 * (sp.sqrt(self.d_a1 ** 2 - self.d_b1 ** 2) + (self.z_2 / abs(self.z_2)) * sp.sqrt(self.d_a2 ** 2 - self.d_b2 ** 2)) - self.a * sp.sin(self.alpha_r)) / (sp.pi * self.m * sp.cos(self.alpha))
        return p_gem

    def E21_27A_NominalToothThickness(self):
        """ s_y = d_y * ((pi + 4 * x * tan(alpha) / (2 * z) + inv(alpha) - inv(alpha_y))"""
        p_gem = self.d_y * ((sp.pi + 4 * self.x * sp.tan(self.alpha) / (2 * self.z)) + self.inv(self.alpha) - self.inv(self.alpha_r))
        return p_gem

    def E21_27B_NominalToothThickness(self):
        """ s_y = d_y * (s/d + inv(alpha) - inv(alpha_r)"""
        p_gem = self.d_y * (self.s / self.d + self.inv(self.alpha) - self.inv(self.alpha_r))
        return p_gem

    def E21_28_ToothThicknessTopCircle(self):
        """ s_a = d_a * (s/d + inv(alpha) - inv(alpha_r)"""
        p_gem = self.d_a * (self.s / self.d + self.inv(self.alpha) - self.inv(self.alpha_r))
        return  p_gem

    def E21_29_DiameterWhereToothSharpens(self):
        """ d_sp = d * cos(alpha) / cos(alpha_sp)"""
        p_gem = self.d  * sp.cos(self.alpha) / sp.cos(self.alpha_sp)
        return p_gem

    def E21_30A_ToothSpace(self):
        """ e_y = d_y * ((pi - 4 * x * tan(alpha))/ (2*z) - inv(alpha) + inv(alpha_r)"""
        p_gem = self.d_y * ((sp.pi - 4 * self.x * sp.tan(self.alpha)) / (2 * self.z) - self.inv(self.alpha) + self.inv(self.alpha_r))
        return p_gem

    def E21_30B_ToothSpace(self):
        """ e_y = d_y - (e/d - inv(alpha) + inv(alpha_r))"""
        p_gem = self.d_y - (self.e / self.d - self.inv(self.alpha) + self.inv(self.alpha_r))
        return p_gem

    def E21_31_InvAlpha_r(self):
        """inv(alpha_r) = 2 * (x_1 + x_2) / (z_1 + z_2) * tan(alpha) + inv(alpha)"""
        p_gem = 2 * (self.x_1 + self.x_2) / (self.z_1 + self.z_2) * sp.tan(self.alpha) + self.inv(self.alpha)
        return p_gem

    def E21_32A_SumOfProfileShiftFactors(self):
        """ Sigma (x) = x_1 + x_2 """
        p_gem = self.x_1 + self.x_2
        return p_gem

    def E21_32B_SumOfProfileShiftFactors(self):
        """ Sigma(x) = (inv(alpha_r) - inv(alpha)) / (2 * tan(alpha) * (z_1 + z_2)"""
        p_gem = (self.inv((self.alpha_r) - self.inv(self.alpha)) / (2 * sp.tan(self.alpha) * (self.z_1 + self.z_2)))
        return p_gem

    def E21_33_ProfileShiftFactor(self):
        """ x_1 = (x_1 + x_2) / 2 + (0.5 - (x_1 - x_2) / 2) * log(u) / log(z_1 * z_2 / 100)"""
        p_gem = (self.x_1 + self.x_2) / 2 + (0.5 - (self.x_1 - self.x_2) / 2) * sp.log(self.u) / sp.log(self.z_1 * self.z_2 / 100)
        return p_gem

    def E21_58AA_MaxTheoreticalRadialPlay(self):
        """ j_tMax = -SigmaA_sti +_Deltaj_ae"""
        p_gem = -1 * self.SigmaA_sti + self.Deltaj_ae
        return p_gem

    def E21_58AB_MaxTheoreticalRadialPlay(self):
        """ j_tMax = -(SigmaA_sni / cos(beta)) + Deltaj_ae"""
        p_gem = -1* (self.SigmaA_sni / sp.cos(self.beta)) + self.Deltaj_ae
        return p_gem

    def E21_58BA_MinTheoreticalRadialPlay(self):
        """ j_tMin = -SigmaA_ste +_Deltaj_ai"""
        p_gem = -1 * self.SigmaA_ste + self.Deltaj_ai
        return p_gem

    def E21_58BB_MinTheoreticalRadialPlay(self):
        """ j_tMin = -(SigmaA_sne / cos(beta)) + Deltaj_ai"""
        p_gem = -1 * (self.SigmaA_sne / sp.cos(self.beta)) + self.Deltaj_ai
        return p_gem

    def E21_59_ChangeInPlay(self):
        """ Delataj_a = 2 * A_a * tan(alpha_n) / cos(beta)"""
        p_gem = 2 * self.A_a * sp.tan(self.alpha_n) / sp.cos(self.beta)
        return p_gem

    def E21_60_ToothWidth(self):
        """ W_k = m_n * cos(alpha_n) * [(k-0.5) + 2 * x * m_n * sin(alpha_n)"""
        p_gem = self.m_n * sp.cos(self.alpha_n) * ((self.k - 0.5) + 2 * self.x * self.m_n * sp.sin(self.alpha_n))
        return p_gem

    def E21_61_NumberOfTeethForMeasurement(self): # not clear on name for function
        """ k = z_n * alpha_n / 180 + 0.5"""
        p_gem = self.z_n * self.alpha_n / 180 + 0.5
        return p_gem

    def E21_62A_TotalGearRation(self): # The function is a product of 'n' terms. not sure how to represent in code form
        """ i = i_1 * i_2 * ... * i_n"""
        p_gem = self.i_1 * self.i_2 * self.i_n
        return p_gem

    def E21_62B_TotalGearRation(self):  # The function is a product of 'n' terms. not sure how to represent in code form
        """ u = u_1 * u_2 * ... * u_n"""
        p_gem = self.u_1 * self.u_2 * self.u_n
        return p_gem

    def E21_63A_ModulusForPinionOnShaft(self):
        """ m__1_n = (1.8 * d_shaft * cos(beta))/(z_1 - 2.5)"""
        p_gem = (1.8 * self.d_shaft * sp.cos(self.beta)) / (self.z_1 - 2.5)
        return p_gem

    def E21_63B_ModulusForPinionShaft(self):
        """ m__1_n = (1.1 * d_shaft * cos(beta)) / (z_1 - 2.5)"""
        p_gem = (1.1 * self.d_shaft * sp.cos(self.beta)) / (self.z_1 - 2.5)
        return p_gem

    def E21_64_ModulusForPinionGearWithDistanceAToShaft(self):
        """ m__2_n = (2 * a * cos(beta)) / (1 + i - z_i)"""
        p_gem = (2 * self.a * sp.cos(self.beta)) / (1 + self.i - self.z_1)
        return p_gem

    def E21_65A_ModulusForCylindricalGearsWithHardenedFlanks(self):
        """ m__2_n = 1.85 * (((T_1eq * cos**2(beta)) / (z_1**2 * psi_d * sigma_(F)Vlim1))**(1/3) """
        p_gem = 1.85 * ((self.T_1eq * sp.cos(self.beta)**2) / (self.z_1 ** 2 * self.psi_d * self.sigma_FVLim1)) ** (1/3)
        return p_gem

    def E21_65B_ModulusForCylindricalGearsUnhardenedAndUntampered(self):
        """ m__3_n = (95 * cos(beta)) / z_1 *((T_1 * (u + 1)) / (psi_d * siga_HLim**2 * u))**(1/3)"""
        p_gem = (95 * sp.cos(self.beta)) / self.z_1 * ((self.T_1 * (self.u +1)) / (self.psi_d * self.sigma_HLim**2 * self.u))**(1/3)
        return p_gem

    def E21_66_NominalTorque(self):
        """ T_1Nom = F_t1 * d_r1 / 2"""
        p_gem = self.F_t1 * self.d_r1 / 2
        return p_gem

    def E21_67A_NominalPeripheralForce(self):
        """ F_t1,2 = F_bn1,2 * cos(alpha_r)"""
        p_gem = self.F_bn12 * sp.cos(self.alpha_r)
        return p_gem

    def E21_67B_NominalPeripheralForce(self): # in equation in book they use the notation of T_1,2 to refer to different possible circles from my understanding. I used the variable refreing to a specific circle in the calculations.
        """ F_t1,2 = 2 * T_1,2 / d_r12"""
        p_gem = 2 * self.T_1 / self.d_r1
        return p_gem

    def E21_68_ToothNormalForce(self):
        """ F_bn1,2 = F_t1,2 / cos(alpha_r)"""
        p_gem = self.F_t1 / sp.cos(self.alpha_r)
        return p_gem

    def E21_69_RadialForces(self):
        """ F_r1,2 = F_t1,2 * tan(alpha_r)"""
        p_gem = self.F_t1 * sp.tan(self.alpha_r)
        return p_gem

    def E21_70_NominalTransversalForces(self):
        """ F_t1,2 = 2 * T_nom1,2 / d_r1,2"""
        p_gem = 2 * self.T_1Nom / self.d_r1
        return p_gem

    def E21_71_RadialForce(self):
        """ F_r1,2 = F_t1,2 * tan(alpha_n) / cos(beta)"""
        p_gem = self.F_t1 * sp.tan(self.alpha_n) / sp.cos(self.beta)
        return p_gem

    def E21_72_AxialForces(self):
        """ F_a1,2 = F_t1,2 * tan(beta)"""
        p_gem = self.F_t1 * sp.tan(self.beta)
        return p_gem

    def E21_73_DynamicFactor(self):
        """ K_v = 1 + ( K_1 / (K_A * (F_t/b)) + K_2) * K_3"""
        p_gem = 1 + self.K_3 * (self.K_1 / (self.K_A * (self.F_t / self.b)) + self.K_2)
        return p_gem

    def E21_74A_LoadConcentrationFactors(self):
        """K_HBeta = 1 + 10 * F_betaY / (F_m / b)"""
        p_gem = 1 + 10 * self.F_betaY / (self.F_m / self.b)
        return p_gem

    def E21_74B_LoadConcentrationFactors(self):
        """ K_HBeta = 2 * sqrt(10 * F_betaY / (F_m / b)"""
        p_gem = 2 * sp.sqrt(10 * self.F_betaY / (self.F_m / self.b))
        return p_gem

    def E21_74C_LoadConcentrationFactors(self):
        """ K_VBeta = K_HBeta ** N_F"""
        p_gem = self.K_HBeta ** self.N_f
        return p_gem

    def E21_75_FlankLineDeviationDueToDeformation(self):
        """ 0.023 * (F_m / b) * [|0.7 + k__1 * (l *s / d_1 ** 2) * ( d_1 / d_as)**4| + 0.3] (b / d_1)**2"""
        p_gem = 0.023 * (self.F_m / self.b) * (abs(0.7 + self.K__1 * (self.l * self.s / self.d_1**2) * (self.d_1 / self.d_as)**4) + 0.3) * (self.b / self.d_1)**2
        return p_gem

    def E21_76A_ManufacturingDependantFlankDeviation(self):
        """ f_ma = c * f_HBeta"""
        p_gem = self.c * self._f_HBeta
        return p_gem

    def E21_76B_ManufacturingDependantFlankDeviation(self):
        """ f_ma = c * 4.16 * b**0.14 * q_H"""
        p_gem = self.c * 4.16 * self.b**0.14 * self.q_H
        return p_gem

    def E21_77_EffectiveDeviationFromFlank(self):
        """ F_betaX = f_ma + 1.33 * f_as"""
        p_gem = self.f_ma + 1.33 * self.f_as
        return p_gem

    def E21_78_EffectiveFlankLineDeviation(self):
        """ F_betaY = F_betaX - y_beta"""
        p_gem = self.F_betaX - self.y_beta
        return p_gem

    def E21_79_ExponentForFactorOfToothBase(self):
        """ N_v = (b/h)**2 / [1 + b/h + (b/h)**2]"""
        p_gem = (self.b / self.h)**2 / (1 + self.b / self.h + (self.b / self.h)**2)
        return p_gem

    def E21_80A_PerimeterFactor(self):
        """ K_HAlpha = epsilon_gama/2 * (0.9+(0.4*c_gama*(f_pe-y_alpha)/(F_tH/b))"""
        p_gem = self.epsilon_gama / 2 * (0.9 + (0.4 * self.c_gama * (self.f_pe - self.y_alpha)/(self.F_tH / self.b)))
        return p_gem

    def E21_80B_PerimeterFactor(self):
        """ K_HAlpha = 0.9 + 0.4 * sqrt(2*(epsilon_gama-1)/epsilon_gama) * (c_gama*(f_pe - y_alpha)/(F_tH / b))"""
        p_gem = 0.9 + 0.4 * sp.sqrt(2*(self.epsilon_gama - 1)/self.epsilon_gama) * (self.c_gama * (self.f_pe - self.y_alpha) / (self.F_tH / self.b))
        return p_gem

    def E21_81A_LoadInfluenceForToothBaseStrength(self):
        """ K_VTot = K_A * K_V * K_Vx * K_Vbeta"""
        p_gem = self.K_A * self.K_V * self.K_VX * self.K_VBeta
        return p_gem

    def E21_81B_LoadInfluenceForContactStrength(self):
        """ K_HTot = K_A * K_V K_HAlpha * K_HBeta"""
        p_gem = self.K_A * self.K_V * self.K_HAlpha * self.K_HBeta
        return  p_gem

    def E21_82_ApproxElocalToothStrength(self):
        """ sigma_V0 = F_t / (b * m_n) * Y_Va * Y_Ka * Y_epsilon * Y_beta"""
        p_gem = self.F_t / (self.b * self.m_n) * self.Y_Va * self.Y_Ka * self.Y_epsilon * self.Y_beta
        return p_gem

    def E21_83_ToothBaseStress(self):
        """ sigma_V1,2 = sigma_V01,2 * K_VTot"""
        p_gem = self.sigma_V012 * self.K_VTot
        return p_gem

    def E21_84A_PermissibleToothBaseStress(self):
        """ sigma_VG = sigma_VLim * Y_ST * Y_NT * Y_delatrelT * Y_RRelT * Y_X"""
        p_gem = self.sigma_VLim * self.Y_ST * self.Y_NT * self.Y_deltaRelT * self.Y_RRelT * self.Y_X
        return p_gem

    def E21_84B_PermissibleToothBaseStress(self):
        """ sigma_VG = 2 * sigma_VLim * Y_NT * Y_X"""
        p_gem = 2 * self.sigma_VLim * self.Y_NT * self.Y_X
        return p_gem

    def E21_85_SafetyOfToothBaseStrength(self):
        """ S_V1,2 = sigma_VG1,2 / sigma_V1,2"""
        p_gem = self.sigma_VG12 / self.sigma_V12
        return p_gem

    def E21_86A_MaximumStress(self):
        """ sigma_HMax = sqrt(F * E / (2 * pi * (1- nu)**2 * b * rho)"""
        p_gem = sp.sqrt(self.F * self.E / (2 * sp.pi * (1 - self.nu**2) * self.b * self.rho))
        return p_gem

    def E21_86B_MaximumStress(self):
        """ sigma_HMax  sqrt(0.175 * F * E / (b * rho))"""
        p_gem = sp.sqrt(0.175 * self.F * self.E / (self.b * self.rho))
        return  p_gem

    def E21_87A_StressInPoleC(self):
        """ sigma_HC = sqrt(0.175 * F_t * E * (u+1) * 2 / (b * d_1 * u * cos**2(alpha) * tan(alpha_r)))"""
        p_gem = sp.sqrt(0.175 * self.F_t * self.E * (self.u + 1) * 2 / (self.b * self.d_1 * self.u * sp.cos(self.alpha)**2 * sp.tan(self.alpha_r)))
        return p_gem

    def E21_87B_StressInPoleC(self):
        """ sigma_HC = sqrt(F_t * (u + 1) / ( b * d_1 * u)) * Z_H * Z_E"""
        p_gem = sp.sqrt(self.F_t * (self.u + 1) / (self.b * self.d_1 * self.u)) * self.Z_E * self.Z_H
        return p_gem

    def E21_88A_FlankStressInPoleC(self):
        """ sigma_H0 = sigma_HC * Z_epsilon * Z_beta """
        p_gem = self.sigma_HC * self.Z_epsilon * self.Z_beta
        return p_gem

    def E21_88B_FlankStressInPoleC(self):
        """ sigma_H0 = Z_H * Z_E * Z_epsilon * Z_beta * sqrt(F_t * (u+1) / (b * d_1 * u)"""
        p_gem = self.Z_H * self.Z_E * self.Z_epsilon * self.Z_beta * sp.sqrt(self.F_t * (self.u + 1) / (self.b * self.d_1 * self.u))
        return p_gem

    def E21_89_ContactStressOnRollCircle(self):
        """ Sigma_H = sigma_H0 * K_HTot"""
        p_gem = self.sigma_H0 * self.K_HTot
        return p_gem

    def E21_90A_AllowableStressFlank(self):
        """ sigma_HG = sigma_HLim * Z_L * Z_vir * Z_R * Z_M * Z_X"""
        p_gem = self.sigma_HLim * self.Z_L * self.Z_vir * self.Z_R * self.Z_M * self.Z_X
        return p_gem

    def E21_90B_SafetyFlankStrength(self):
        """ S_H1,2 = sigma_HG1,2 / sigma_H"""
        p_gem = self.sigma_HG12 / self.sigma_H
        return p_gem








