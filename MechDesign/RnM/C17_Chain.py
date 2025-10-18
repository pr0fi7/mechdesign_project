
"""
This folder will contain - for student's reference - expressions from Roloff & Matek. 
These expressions are only for internal use and e.g. for studying/answering the exam. 
These expressions may never and under no exception be distributed or shared to anyone else 
and serve only for the purpose of students that follow the KULeuven courses on dimensioning of machine parts 
By using these files you commit to adhering to this condition ! 
"""

import sympy as sp
from . import MySymbol
from ..Units.Units import deg_, h_, m_, s_


class Chain:
    
    
    MySymbolDict = {'i'         : '[] gear ratio',\
                    'n_1'       : '[] speed of driving chain wheel',\
                    'n_2'       : '[] speed of driven chain wheel',\
                    'z_2'       : '[] nr of teeth of the driven chain wheel',\
                    'z_1'       : '[] nr of teeth of the driving chain wheel',\
                    'd_2'       : '[] pitch diameter of the driven chain wheel',\
                    'd_1'       : '[] pitch diameter of the driving chain wheel',\
                    'T_1'       : '[Nm] torque applied to the driving chain wheel',\
                    'tau'       : '[deg] pitch angle',\
                    'z'         : '[] nr of teeth',\
                    'd'         : '[mm] pitch diamter',\
                    'p'         : '[mm] chain pitch {table 17-1}',\
                    'd_v'       : '[mm] foot diameter',\
                    'dprime_1'  : '[mm] maximum rol diamter {table 17-1}',\
                    'd_a'       : '[mm] top circle diameter',\
                    'd_amax'    : '[mm] maximum top circle diameter',\
                    'd_amin'    : '[mm] minimum top circle diameter',\
                    'd_s'       : '[mm] freewheel diamter beneath the pitch circle',\
                    'g_1'       : '[mm] chain height, {table 17-1}',\
                    'r_4'       : '[mm] radius of the sprocket teeth {table 17-2}',\
                    'F'         : '[mm] free play dimension of the chain {table 17-2}',\
                    'P_D'       : '[W] diagram power',\
                    'P_1'       : '[W] input power',\
                    'P_2'       : '[W] output power',\
                    'eta'       : '[] chain drive efficiency ',\
                    'K_A'       : '[] working condtion factor {table 3-5}',\
                    'f_1'       : '[] power factor on number of teeth {table 17-5}',\
                    'f_2'       : '[] power factor on shaft to shaft distance {table 17-6}',\
                    'f_3'       : '[] power factor on shackle shape {exp 17.7}',\
                    'f_4'       : '[] power factor on number of chain wheels {exp 17.7}',\
                    'f_5'       : '[] power factor on lifetime expectancy ',\
                    'f_6'       : '[] power factor on working conditions {table 17-7}',\
                    'n'         : '[] number of chain wheels, the chain rolls on',\
                    'L_h'       : '[h] lifetime expectancy of the chain',\
                    'a'         : '[mm] shaft to shaft distance',\
                    'a_0'       : '[mm] desired shaft tot shaft distance',\
                    'X_0'       : '[] initial estimate of the number of shackels',\
                    'X'         : '[] number of shakels in the chain',\
                    'L'         : '[mm] total length of the chain',\
                    'f'         : '[mm] chain sagging distance',\
                    'f_rel'     : '[] relative chain sagging distance',\
                    'l_T'       : '[mm] length chain section under tension',\
                    'F_t'       : '[N] mathematical tangential chain force',\
                    'v'         : '[m/s] chain speed',\
                    'F_c'       : '[N] centrifugal force',\
                    'q'         : '[kg/m] chain weight / lenght',\
                    'F_s'       : '[N] supporting force horizontal chain',\
                    'F_G'       : '[N] gravity force on the chain',\
                    'g'         : '[m/s²] gravitationl constant, to be set by user',\
                    'F_sb'      : '[N] supporting force for top chain wheel',\
                    'F_so'      : '[N] supportinf force for bottom chain wheel',\
                    'Fprime_s'  : '[] relative load factor ',\
                    'psi'       : '[deg] angle between horizontal and chain contact points line',\
                    'epsilon_0' : '[deg] angle between chain contact points line and chain wheel center to center line',\
                    'delta'     : '[deg] angle between horizontal and chain wheel center to center line',\
                    'F_a'       : '[N] Support force for the shaft, standing still, horizontal chain',\
                    'F_Ab'      : '[N] Support force for the shaft, standing still, top chain wheel',\
                    'F_Ao'      : '[N] Support force for the shaft, standing still, bottom chain wheel',\
                    'F_res'     : '[N] Supoprt force for the shaft, while rotating, horizontal chain',\
                    
                    }

    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

    def E17_1A_GearRatioSpeed(self):
        """ i = n_1/n_2 """
        i = self.n_1 / self.n_2
        return i

    def E17_1B_GearRatioTeeth(self):
        """ i = z_2/z_1 """
        i = self.z_2 / self.z_1
        return i
    
    def E17_1C_GearRatioDiameter(self):
        """ i = d_2/d_1 """
        i = self.d_2 / self.d_1
        return i    
                               
    def E17_2_PitchAngle(self):
        """ tau = 360°/z """
        tau = 360 * deg_ / self.z 
        return tau
        
    def E17_3A_PitchDiamter(self):
        """ d = p / sin(tau/2) 
            with tau in [deg_]
        """
        d = self.p / sp.sin((self.tau/(180*deg_)*sp.pi) / 2)
        return d

    def E17_3B_PitchDiameter(self):
        """ d = p / sin (180°/z) """
        d = self.p / sp.sin(sp.pi/self.z)
        return d


    def E17_4_FoothDiameter(self):
        """ d_v = d - d'_1 """
        d_v = self.d - self.dprime_1
        return d_v
    
    def E17_5a_TopCircleDiameter(self):
        """ d_a  = d * cos (tau/2) + 0.8*d'_1 
            with tau in [deg_]
        """
        d_a = self.d * sp.cos(self.tau/(180*deg_)*sp.pi / 2) + 0.8 * self.dprime_1
        return d_a
    
    def E17_5b_MaxTopCircleDiameter(self):
        """ d_amax = d +1.25*p-d'_1 """
        d_amax = self.d+1.25*self.p - self.dprime_1
        return d_amax
    
    def E17_5c_MinTopCircleDiameter(self):
        """ d_amin = d + (1 -1.6/z)*p - d'_1 """
        d_amin = self.d + (1-1.6/self.z)*self.p -self.dprime_1
        return d_amin
    
    def E17_6a_FreeWheelDiameter(self):
        """ d_s = p/tan(tau/2) - 1.05*g_1 -2*r_4-1 """
        d_s = self.p / sp.tan(self.tau/(180*deg_)*sp.pi/2) - 1.05*self.g_1 -2*self.r_4 - 1 
        return d_s
    
    def E17_6b_FreeWheelDiameter(self):
        """ d_s = d - 2*F """
        d_s = self.d - 2*self.F 
        return d_s
    def E17_7_DiagramPower(self):
        """ P_D = (K_A*P_1*f_1)/f_2/f_3/f_4/f_5/f_6"""
        P_D = (self.K_A*self.P_1*self.f_1)/(self.f_2*self.f_3*self.f_4*self.f_5*self.f_6)
        return P_D
    
    def E17_7_hA_PowerEfficiency(self):
        """ P_1 = P_2/eta """
        P_1 = self.P_2 / self.eta
        return P_1
    
    def E17_7_hB_NrOfWheelsCorrection(self):
        """ f_4 = 0.9**(n-2) """
        f_4 = 0.9**(self.n-2) 
        return f_4
    
    def E17_7_hC_LifetimeExpectancyFactor(self):
        """ f_5 = (15000/L_h)**(1/3) """
        f_5 = (15000*h_/self.L_h)**(1/3)
        return f_5
    
    def E17_9_InitalNrOfShackels(self):
        """ X_0= 2*a_0/p + (z_1+z_2)/2 + ((z_2-Z1)/(2*pi))**2 * p / a_0 """
        X_0=2*self.a_0/self.p + (self.z_1+self.z_2)/2 + (self.p/self.a_0)*((self.z_2-self.z_1)/(2*sp.pi))**2
        return X_0
    
    def E17_10_ShaftDistance(self):
        """ a = (p/4) * ( (X-(z_1+z_2)/2) + sqrt( (X-(z_1+z_2)/2)**2 - 2*((z_2-z_1)/pi)**2   ) ) """
        a = (self.p/4) * ( (self.X-(self.z_1+self.z_2)/2) + sp.sqrt( (self.X-(self.z_1+self.z_2)/2)**2 - 2*((self.z_2-self.z_1)/sp.pi)**2   ) )
        return a
    
    def E17_12_NrOfShackels(self):
        """ X = L/p """
        X = self.L / self.p
        return X
    
    def E17_13_RelChainSagging(self):
        """ f_rel = f/l_T """
        f_rel = self.f / self.l_T
        return f_rel
    
    def E17_14A_TangentialForce(self):
        """ F_T = P_1 / v """
        F_T = self.P_1 / self.v 
        return F_T
    
    def E17_14B_TangentialForce(self):
        """ F_T = T_1 / (d_1/2) """
        F_T = self.T_1 / (self.d_1/2)
        return F_T
    
    def E17_14_h_Speed(self):
        """ v = d_1/2*n_1 
        in order to correctly calculate speed with units, no pi in the formula
        """
        v = self.d_1/2 * self.n_1
        return v
    
    def E17_15_CentrifugalForce(self):
        """ F_c = q* v**2 """
        F_c = self.q*self.v**2 
        return F_c
    
    def E17_16A_HorizontalSupportForce(self):
        """ F_s = F_G*l_T/(8*f) """
        F_s = self.F_G * self.l_T / (8*self.f )
        return F_s
        
    def E17_16B_HorizontalSupportForce(self):
        """ F_s = q*g*l_T/(8*f_rel)"""
        F_s = (self.q * self.l_T * self.g)/(8*self.f_rel)
        return F_s
    
    def E17_17_TopWheelSupportForce(self):
        """ F_sb = q*g*l_T*(Fprime_s+sin(psi)) """
        F_sb = self.q * self.g * self.l_T *(self.Fprime_s+ sp.sin(self.psi/(180*deg_)*sp.pi))
        return F_sb
    
    def E17_17_hA_psi(self):
        """ psi = delta - epsilon_0 """
        psi = self.delta - self.epsilon_0
        return psi
    
    def E17_17_hB_epsilon_0(self):
        """ epsilon_0 = asin((d_2-d_1)/(2*a)) """
        epsilon_0 = sp.asin((self.d_2-self.d_1)/(2*self.a ))
        return epsilon_0
    
    def E17_18_BottomWheelSupportForce(self):
        """ F_so = q*g*l_T*Fprime_s """
        F_so = self.q*self.g*self.l_T*self.Fprime_s
        return F_so
    
    def E17_19A_ShaftSupportForceHorizontal(self):
        """ F_a = F_t*K_A+2*F_s """
        F_a = self.F_t * self.K_A + 2*self.F_s 
        return F_a
    
    def E17_19B_ShaftSupportForceTopChainWheel(self):
        """ F_a = F_t*K_A+2*F_sb """
        F_a = self.F_t * self.K_A + 2*self.F_sb 
        return F_a
    
    def E17_19C_ShaftSupportForceBottomChainWheel(self):
        """ F_a = F_t*K_A+2*F_so """
        F_a = self.F_t * self.K_A + 2*self.F_so
        return F_a
        
    def E17_20_ResultingSupportForce(self):
        """ F_res = F_t * K_A +F_c +F_s """
        F_res = self.F_t * self.K_A +self.F_c+self.F_s 
        return F_res
    
    
    
        