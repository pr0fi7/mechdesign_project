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


class Belt:
    
    
    MySymbolDict = {'F_t'       : '[N] the minimal frictional force needed',\
                    'F_N'       : '[N] normal force needed',\
                    'mu'        : '[] friction coeficient {table 16-1}',\
                    'muprime'   : '[] muprime corrected friction coefficient cfr equation 16.1',\
                    'alpha'     : '[] wedge angle of V belt, {table 16-13 & 16-14}',\
                    'F_1'       : '[N] the force in the pulling part of the belt',\
                    'F_2'       : '[N] the force in the "free" part of hte belt ',\
                    'm'         : '[] belt force ratio',\
                    'beta_1'    : '[] Circumferentialangle of the belt of the small disk',\
                    'kappa'     : '[] utilisation factor',\
                    'F_c'       : '[N] centrifugal force',\
                    'A_S'       : '[mm²] belt cross section ',\
                    'rho'       : '[kg/dm³] belt denistiy {tabel 16-1} ',\
                    'v'         : '[m/s] belt velosity',\
                    'F_a'       : '[N] radial load on the axle holding the pulley, not taking centrefugal force into account',\
                    'F_a0'      : '[N] total radial load on the axle holding the pulley',\
                    'k'         : '[] shaft load factor {tabel 16-5}',\
                    'psi'       : '[%] belt slip',\
                    'v_1'       : '[m/s] belt speed @pulley 1 (v_1 : d_1 * pi *n_1)',\
                    'v_2'       : '[m/s] belt speed @pulley 2 (v_2 : d_2 * pi *n_2)',\
                    'i'         : '[] transmission ratio (equivalent to gear ratio in case of gears)',\
                    'd_1'       : '[m] diameter pulling disk',\
                    'd_2'       : '[m] diameter pulled disk',\
                    't'         : '[mm] belt thickness ',\
                    'sigma_1'   : '[N/mm²] tension due to pulling force, side 1',\
                    'sigma_2'   : '[N/mm²] tension due to pulling force, side 2',\
                    'sigma_b'   : '[N/mm²] tension due to bending allong pulley',\
                    'sigma_c'   : '[N/mm²] tension due to centrifugal forces',\
                    'sigma_N'   : '[N/mm²] useable tension',\
                    'sigma_tot' : '[N/mm²] total tension',\
                    'sigma_max' : '[N/mm²] maximum permissible tension {tabel 16-1}',\
                    'sigma_S'   : '[N/mm²] additional tension due to crossing belts ',\
                    'E_b'       : '[N/mm²] ideal elasticity modulus {tabel 16-1}',\
                    'b'         : '[mm] belt width, table {tabel 16-9}',\
                    'e'         : '[mm] belt length',\
                    'P'         : '[W] Power to transmit (nominal value)',\
                    'Pprime'    : '[W] Power to transmit with integrated safety factor',\
                    'K_A'       : '[] correctionuse-factor {tabel 3-5}',\
                    'T_nom'     : '[Nm] nominal Torque ',\
                    'v_opt'     : '[m/s] optimal belt speed',\
                    'd_dg'      : '[mm] reference diameter of larger disk, note explenation equation 16-36',\
                    'd_dk'      : '[mm] reference diameter of smaller disk, note explenation equation 16-36',\
                    'z_g'       : '[] number of teeth on bigger disk, in case of toothed belt',\
                    'z_k'       : '[] number of teeht on smaller disk, in case of toothed belt',\
                    'beta_k'    : '[] Circumferentialangle of the belt of the small disk',\
                    "eprime"    : "[mm] e' :  the global shaft distance",\
                    "eprime_LL" : "[mm] lower limit on e' :  the global shaft distance",\
                    "eprime_UL" : "[mm] Upper Limit on e' :  the global shaft distance",\
                    'p'         : '[mm] pitch of the toothed belt',\
                    'L_d'       : '[mm] belt length, normed length',\
                    'L'         : '[mm] belt length, reference length',\
                    'Lprime'    : '[mm] belt length, theoretical length',\
                    'x'         : '[mm] tensioning length',\
                    'y'         : '[mm] mounting length',\
                    'bprime'    : '[mm] belt width, theoretical length, ',\
                    'Fprime_t'  : '[N/mm] specific load {tabel 16-8}',\
                    'P_N'       : '[W] nominal power of one belt or rib, {tabel 16-15} ',\
                    'U_z'       : '[W] ratio bonus when i>1, {tabel 16-16} ',\
                    'c_1'       : '[] angle factor when beta_1 < :  180°, {tabel 16-17}',\
                    'c_2'       : '[] length factor depending on L_d, {tabel 16-17}',\
                    'P_spec'    : '[W/mm] specific power toothed belt {tabel 16-20}',\
                    'T_max'     : '[Nm] maximum torque to be transmitted ',\
                    'T_spec'    : '[Nm/mm] specific transmissible torque {tabel 16-20}',\
                    'epsilon_1' : '[%] required stress in standstill {tabel 16-8}',\
                    'epsilon_2' : '[%] stress due to centrifugal force {tabel 16-10}',\
                    'k_1'       : '[] belt type factor {tabel 16-6}',\
                    'epsilon_tot' : '[%] total belt stress',\
                    'z'         : '[] number of belts',\
                    'z_i'       : '[] number of interacting teeth',\
                    'd_w'       : '[m] effective diameter',\
                    'd_d'       : '[m] reference diameter',\
                    'h_b'       : '[m] reference height ',\
                    'n'         : '[rpm] rotational speed',\
                    'zz'        : '[] number of disks the belt passes',\
                    'f_B'       : '[s-1] bending frequency maximum values{tabel 16-1;16-2;16-3}',\
                    'F_max'     : '[N] maximum permissible load toothed belt',\
                    'n_1'       : '[rpm] rotational speed pulling disk ',\
                    'n_2'       : '[rpm] rotational speed pulled disk',\
                    'T_1'       : '[Nm] torque @ pulling disk',\
                    'T_2'       : '[Nm] torque @ pulled disk',\
                    'eta'       : '[%] power efficiency, during caluculations 100% efficency is assumed',\
                    }

    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------


                               
    def E16_1_FrictionForce(self,MyHelp=False):
        """ F_t = mu * F_N  """
        F_t = self.mu * self.F_N
        return F_t    
        
    def E16_1B_Set2muprime(self):
        """ muprime = mu / sin(alpha/2 )
            remeber convertion from ° to rad
            be carefull when using this function
            you can easily create an unwanted loop e.g. 
            b.mu = b.E16_1B_SetMuprime  
            will create a loop because self.mu is inside of the equation itself
            solution is to add an explicit reïnitialisation to the symbol mu
            b.muprime = b.E16_1B_SetMuprime
            b.mu = b.muprime
        """
        muprime = self.mu / sp.sin( self.alpha/2)
        return muprime
        
    def E16_2_CircumferenceForce(self) :
        """ F_t = F_1 - F_2 """
        F_t = self.F_1 - self.F_2
        return F_t
    
    def E16_3_Eytelwein(self):
        """ m = exp (mu*beta1) """
        m = sp.exp(self.mu *  self.beta_1)
        return m
    
    def E16_4_CircumferenceForce(self):
        """ F_t_formula = F_1 * kappa """
        F_t = self.F_1 * self.kappa
        return F_t
        
    def E16_4B_kappa(self):
        """ kappa = (m -1) / m """
        kappa = (self.m - 1)/self.m
        return  kappa
    
    def E16_5_CentrifugalForce(self):
        """ F_c = A_S * 1E-6 * 1E3 * rho * v**2 """
        F_c = self.A_S * 1E-6 * 1E3 * self.rho * self.v**2
        return  F_c
    
    def E16_6_AxleLoad(self):  
        """ F_a = F_t *  sqrt( m**2 + 1 -2 * m * cos( beta_1 / 180 * pi ) ) / ( m - 1 ) """
        F_a = self.F_t *  sp.sqrt( self.m**2 + 1 -2 * self.m * sp.cos( self.beta_1) ) / ( self.m - 1 ) 
        return F_a
    
    def E16_6B_AxleLoad(self):
        """ F_a = k * F_t """
        F_a = self.k * self.F_t
        return F_a
    
    def E16_7_TotalAxleLoad(self):
        """ F_a0 = F_a + F_c """
        F_a0 = self.F_a + self.F_c
        return F_a0
        
    def E16_8_BeltSlip(self):
        """ psi = 100 * ( v_1 - v_2 ) / v_2 """
        psi = 100 * ( self.v_1 - self.v_2 ) / self.v_2
        return psi
        
    def E16_9_TransmissionRatio(self):
        """ i = ( d_2 + 1E-3 * t ) / (d_1 + 1E-3 * t ) * 100 / ( 100 - psi ) """
        i = ( self.d_2 + 1E-3 * self.t ) / (self.d_1 + 1E-3 * self.t ) * 100 / ( 100 - self.psi )
        return i
        
    def E16_10_TransmissionRatioApprox(self):
        """ i = d_2 / d_1 """
        i = self.d_2 / self.d_1
        return i
    
    def E16_11_TensionPulling(self):
        """ simga_1 = F_1 / A_S """
        sigma_1 = self.F_1 / self.A_S
        return sigma_1
        
    def E16_12_TensionBending(self):
        """ sigma_b = E_b * t / d_1 
            d is replaced with d_1, the smaller of the two disks
        """
        sigma_b = self.E_b * self.t / self.d_1        
        return sigma_b
        
    def E16_13A_TensionCentrifugal(self):
        """ sigma_c = F_c / A_s """
        sigma_c = self.F_c / self.A_S
        return sigma_c
        
    def E16_13B_TensionCentrifugal(self):
        """ sigma_c = rho v**2 """
        sigma_c = self.rho * self.v**2
        return sigma_c

    def E16_14_TensionTotal(self):
        """ sigma_tot = sigma_1 + sigma_b + sigma_c """
        sigma_tot = self.sigma_1 + self.sigma_b + self.sigma_c
        return sigma_tot
        
    def E16_14B1_TensionCrossingBelt(self):
        """ sigma_S = E_b * (b / e)** 2
            the value of E is set to E_b
        """
        sigma_S = self.E_b *  (self.b / self.e)**2 
        return sigma_S
    
    def E16_14B2_TensionHalfCrossingBelt(self):
        """ sigma_S = E_b * b * d_2 / (2 *e**2) 
            the value of E is set to E_b
        """
        sigma_S = self.E_b * self.b * self.d_2 / (2 * self.e**2)
        return sigma_S
        
    def E16_15A_TensionUseable(self):
        """ sigma_N = sigma_1 - sigma_2 """
        sigma_N = self.sigma_1 - self.sigma_2 
        return sigma_N
        
    def E16_15B_TensionUseable(self):
        """ simga_N = sigma_1 * kappa """
        sigma_N = self.sigma_1 * self.kappa
        return sigma_N
        
    def E16_17_TransmittablePower(self):
        """ P = (sigma_max - E_b * (t/d_1)-rho*(v**2)*E-3)*kappa*b*t*v 
            P will be calculated in W and not KW (as is inidicated in the table below the formula in the book)
        """
        P = (self.sigma_max - self.E_b * (self.t/self.d_1)-self.rho*(self.v**2)*1E-3)*self.kappa*self.b*self.t*self.v 
        return P
    
    def E16_18A_OptimalSpeed(self):
        """ v_opt = sqrt(1E3*(sigma_max-E_b*(t/d_1))/(3)*rho) """
        v_opt = sp.sqrt(1E3*(self.sigma_max-self.E_b*(self.t/self.d_1)) / (3*self.rho) )
        return v_opt
        
    def E16_18B_OptimalSpeed(self):
        """ v_opt = sqrt(1E3*(sigma_max-sigma_b)/(3*rho)) """
        v_opt = sp.sqrt(1E3*(self.sigma_max-self.sigma_b) / (3*self.rho)) 
        return v_opt
        
    def E16_19A_TransmissionRatio(self):
        """ i = d_dg / d_dk 
            applies to flat and V belts
        """
        i = self.d_dg / self.d_dk
        return i
        
    def E16_19B_TransmissionRatio(self):
        """ i = z_g/z_k 
            applies to toothed belts
        """
        i = self.z_g / self.z_k    
        return i
        
    def E16_20A_LargeDiskDiameter(self):
        """ d_dg = i * d_dk 
            applies to flat and V belts
        """
        d_dg = self.i * self.d_dk
        return  d_dg
        
    def E16_20B_LargeDiskDiameter(self):
        """ d_dg = i * p / pi * z_k 
            applies to toothed belt
        """
        d_dg = self.i * self.p / sp.pi * self.z_k
        return d_dg
        
    def E16_21A_ShaftDistanceLowerLimit(self):
        """ eprime = 0.7 * (d_dg + d_dk) 
            applies to flat and V belts
        """
        eprime_LL = 0.7 * (self.d_dg + self.d_dk)
        return eprime_LL
        
    def E16_21B_ShaftDistanceLowerLimit(self):
        """ eprime = 0.5 * (d_dg + d_dk) +15
            applies to flat and V belts
        """
        eprime_LL = 0.5 * (self.d_dg + self.d_dk) + 15
        return eprime_LL
        
    def E16_21_ShaftDistanceUpperLimit(self):
        """ eprime = 2 * (d_dg + d_dk) """
        eprime_UL = 2 * (self.d_dg + self.d_dk)
        return eprime_UL
        
    def E16_22_FinalShaftDistance(self):
        """ e = L_d/4 - pi/8 * (d_dg+d_dk) + sqrt( (L_d/4 - pi/8 * (d_dg+d_dk))**2 - ((d_dg-d_dk)**2)/8 ) """
        e = self.L_d/4 - sp.pi/8 * (self.d_dg+self.d_dk) + \
                 sp.sqrt( ((self.L_d/4 - sp.pi/8 * (self.d_dg+self.d_dk))**2)\
                          - ((self.d_dg-self.d_dk)**2)/8 )
        return e
             
    def E16_23_TheoreticalBeltLength(self):
        """ Lprime = 2 * eprime  + pi/2*(d_dg+d_dk) + ((d_dg-d_dl)**2)/(4*eprime) """
        Lprime = 2 * self.eprime\
                    + sp.pi/2*(self.d_dg+self.d_dk)\
                    + ((self.d_dg-self.d_dk)**2)/(4*self.eprime)
        return Lprime
                    
    def E16_24A_circumfranceAngle(self):
        """ beta_k = 2 * acos((d_dg - d_dk)/(2*e))
            applies to flat and V belts
        """
        beta_k = 2 * sp.acos((self.d_dg - self.d_dk)/(2*self.e))
        return  beta_k
        
    def E16_24B_circumfranceAngle(self):
        """ beta_k = 2 * acos(p/pi*(z_g - z_k)/(2*e))
            applies to toothed belt
        """
        beta_k = 2 * sp.acos(self.p/sp.pi*(self.z_g - self.z_k)/(2*self.e))
        return beta_k
    
    def E16_25A_TensioningLength(self):
        """ x = 0.03* L 
            applies to flat and V belts
        """
        x = 0.03* self.L 
        return x
        
    def E16_25B_TensioningLength(self):
        """ x = 0.005 * L 
            applies to toothed belt
        """
        x = 0.005 * self.L 
        return x
        
    def E16_26A_MountingLength(self):
        """ y = 0.015 * L
        applies to flat and V belts
        """
        y = 0.015 * self.L
        return y
        
    def E16_26B_MountingLength(self):
        """ y = 1 * p
        applies to toothed flat belt and V belts
        a minimum value of 1 is implemented, indictated to be between 1 and 2.5 times p
        """
        y = 1 * self.p
        return y
        
    def E16_27_FrictionForce(self):       
        """ F_t = Pprime / v """
        F_t = self.Pprime / self.v 
        return F_t
    
    def E16_28_BeltWidth(self):
        """ bprime = F_t /  Fprime_t 
            applies to flat belt
        """
        bprime = self.F_t / self.Fprime_t
        return bprime
        
    def E16_29_NumberOfBelts(self):
        """ z = Pprime / ( (P__N + U__z)*c_1*c_2)  
            applies to V belts an poly V belts     
        """
        z = self.Pprime / ( (self.P_N + self.U_z)*self.c_1*self.c_2) 
        return z
        
    def E16_30_NumberOfInteractingTeeth(self):
        """ z_i = z_k * beta_k / 360 """
        z_i = self.z_k * self.beta_k / (2*sp.pi)
        return z_i
        
    def E16_31_MinimumBeltWidth(self):
        """ b >= Pprime / (z_k*z_i*T_spec) 
            applies to Toothed Belt
        """
        b = self.Pprime / (self.z_k * self.z_i * self.T_spec)
        return b
        
    def E16_32_MinimumBeltWidth(self):
        """ b >= T_max / ( z_k * z_i * T_spec)
            applies to Toothed Belt
        """
        b = self.T_max / (self.z_k * self.z_i * self.T_spec)
        return b
        
    def E16_33_AxleLoad(self):
        """ F_a = sqrt(F_1**2 + F_2**2-2*F_1*F_2*cos(beta_k)) """
        F_a = sp.sqrt(self.F_1**2 + self.F_2**2 - 2*self.F_1*self.F_2*sp.cos(self.beta_k))
        return F_a
        
    def E16_33B_AxleLoad(self):
        """ F_a = k * F_t """
        F_a = self.k * self.F_t
        return F_a
        
    def E16_34_AxleLoad(self):
        """ F_aO =(epsilon_1+epsilon2)*k_1**bprime """
        F_a0 = (self.epsilon_1 + self.epsilon_2) * self.k_1 * self.bprime
        return F_a0
        
    def E16_35_AxleLoadValeus(self):
        """the limts on F_a0 can be found in expression 16_35"""
        return False
        
    def E16_36_BeltSpeed(self):
        """ v = d_w * pi * n """
        v = self.d_w * sp.pi * self.n
        return v
        
    def E16_36B_EffectiveDiameter(self):
        """ d_w = d + t 
            applies to Extremultus belt {tabel 16-6}
        """
        d_w = self.d_d * self.t 
        return d_w

    def E16_36C_EffectiveDiameter(self):
        """ d_w = d_d + 2*h_b 
            applies to poly V belt
        """
        d_w = self.d_d + 2* self.h_b
        return d_w
        
    def E16_37_BendingFrequency(self):
        """ f_B = v * zz / L_d """
        f_B = self.v * self.zz / self.L_d
        return f_B
        
    def E16_00A_ToothForce(self):
        """ T_max / (d_d / 2) 
            can be found after bending frequency (E16_37), in belt force section
            reference values in tabel {16-19c}
        """
        F_max = self.T_max / (self.d_d / 2)
        return F_max        


        

    
        
        
