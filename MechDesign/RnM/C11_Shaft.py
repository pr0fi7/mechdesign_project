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
from ..Units.Units import deg_


class Shaft:
    

    MySymbolDict = {'sigma_bToel' :'[N/mm²] permissible stress, equation {3.26}'  ,\
              'd' :'[mm] Minimum diameter',\
              'dprime' : '[mm] minimum diameter estimation according to figure {11-21}',\
              'T' :'[N*m] The greatest torque to be transmitted by the shaft',\
              'tau_max' :'[N/m**2] Permissible torsional stress according to equation {3.26}',\
              'd_u' :'[mm] Outside diameter',\
              'd_i' :'[mm] Inside diameter',\
              'k' :'[] Diameter ratio d_i/d/d_u, equation {11.3}',\
              'M_v' :'[N*m] Moment of comparison (imaginary bending moment)',\
              'M_b' :'[N*m] Bending moment for cross section in danger',\
              'sigma_max' :'[N/m**2] permissible stress for the load cases present, equation {11.1}',\
              'phi' :'[] Factor in the calculation of the load ratio; 𝜑 = 1,73',\
              'T_nom' :'[N*m] nominal torsion moment',\
              'P' :'[kW] Largest transferable nominal power',\
              'n' :'[min**-1] at the rated power 𝑃 corresponding (smallest) speed' ,\
              'T_eq' :'[N*m]equivalent torsion moment',\
              'K_A' :'[] Operating factor**2 to take into account the influence of dynamic processes according to table {3-5}',\
              'l' :'[mm] Shaft length',\
              'tau_t' :'[N/mm**2] Torsional stress',\
              'r' :'[mm] Radius of the shaft',\
              'G' :'[N/mm**2] Shear modulus',\
              'T_1' :'[N*mm] Torsion moment in Nmm eq{11.1}',\
              'I_p' :'[mm**4] Polar surface moment of inertia',\
              'f_res' :'[]Resulting deflection, table {11-6}',\
              'f_x' :'[] Deflection in x direction',\
              'f_y' :'[] Deflection in y direction',\
              'alpha_x' :'[°] Angular deflection in x direction',\
              'alpha_y' :'[°] Angluar deflection in y direction',\
              'alpha_res' :'[°] Resulting angular deflection {11-6}',\
              'f_A' :'[] Deflection f under the load F_A',\
              'F_A' :'[] Load F_A',\
              'f_B' :'[] Deflection f under the load F_B',\
              'F_B' :'[] Load F_A',\
              'E' :'[N/m**2] Youngs modulus',\
              'a' :'[]',\
              'a_1' :'[]',\
              'a_2' :'[]',\
              'a_3' :'[]',\
              'd_a1' :'[] Shaft diameter',\
              'd_a2' :'[] Shaft diameter',\
              'd_a3' :'[] Shaft diameter',\
              'd_b1' :'[] Shaft diameter',\
              'd_b2' :'[] Shaft diameter',\
              'd_b3' :'[] Shaft diameter',\
              'b_1' :'[]',\
              'b_2' :'[]',\
              'b_3' :'[]',\
              'f' :'[] Deflection',\
              'alpha' :'[]',\
              'alphaprime' :'[°]  The angle of inclination of the taps in the bearings',\
              'beta' :'[]',\
              'betaprime' :'[°]  The angle of inclination of the taps in the bearings',\
              'omega_k' :'[s**-1] Natural frequency (critical angular velocity)',\
              'c' :'[] Spring constant for elastic bending',\
              'm' :'[] Mass of the spinning disk',\
              'n_k' :'[min**-1] The bending-critical speed',\
              'n_kb' :'[min**-1] The bending-critical speed for bearing and drive shafts',\
              'c_t' :'[] Torsion spring constant from 𝑐𝑡 = 𝐼𝑝 ⋅ 𝐺 ∕l',\
              'J' :'[kg*m**2] Moment of inertia (mass moment 2nd degree) for solid cylinders (drive shafts, discs):𝐽 = (1 ∕ 8) ⋅ 𝑚 ⋅ (𝑑2)',\
              'n_kt' :'[] torsion critical rotational velocity',\
              'J_1' :'[] Moment of inertia (mass moment 2nd degree) as in eq {11.32}',\
              'J_2' :'[] Moment of inertia (mass moment 2nd degree) as in eq {11.32}',\
                  }
        
    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))
            
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------        
        
    def E11_1_MinDiamter(self):
        """ d > cbrt((32*M_b)/(pi*sigma_btoel))"""
        d = sp.cbrt((32*self.M_b)/(sp.pi*self.sigma_bToel))
        return d

    def E11_2_MinOuterDiameter(self):
        """ d_u > cbrt((32*M_b)/(pi*(1-k**4)*sigma_bToel))"""
        d_u = sp.cbrt((32*self.M_b)/(sp.pi*(1-self.k**4)*self.sigma_bToel))
        return d_u

    def E11_3_InnerDiameter(self):
        """ d_i < k*d_u"""  
        d_i = self.k*self.d_u
        return d_i      

    def E11_4_MinDiameter(self):
        """  d = cbrt((32*M_b)/(pi*sigma_bToel) 
        note: no specific variables for diameter at a given location, nor bending moment at given location are provided"""
        d = sp.cbrt((32*self.M_b)/(sp.pi*self.sigma_bToel))
        return d
                               
    def E11_5A_MinDiameter(self):
        """ d=cbrt((16T)/pi*tau_max) """
        d = sp.cbrt(16*self.T)/(sp.pi*self.tau_max)
        return d
    
    def E11_5B_MinDiameter(self):
        """ d= 1.72cbrt((T)/pi*tau_max * F_N) """
        d = 1.72*sp.cbrt(self.T/self.tau_max)
        return d
    
    def E11_6A_OutsideDiameter(self):
        """ d_u=cbrt(16T/(pi*(1-4k^4)*tautmax)) """
        d_u = sp.cbrt(16*self.T/(sp.pi*(1-(self.k)**4)*self.tau_max))
        return d_u
    
    def E11_6B_OutsideDiameter(self):
        """ d_u=1.72*cbrt(16T/((1-4k^4)*tautmax)) """
        d_u = 1.72*sp.cbrt(16*self.T/((1-(self.k)**4)*self.tau_max))
        return d_u
    
    def E11_7A_MomentOfComparison(self):
        """M_v=sqrt(M_b**2+.75(sigma_max/(phi*tau_max)*T)**2) """
        M_v=sp.sqrt(self.M_b**2+0.75*(self.sigma_max*self.T/(self.phi*self.tau_max))**2)
        return M_v
    
    def E11_7B_MomentOfComparison(self):
        """M_v=sqrt(M_b**2+(sigma_max) """
        M_v=sp.sqrt(self.M_b**2+(self.sigma_max*self.T/(2*self.tau_max))**2)
        return M_v
    
    def E11_8A_DiameterShaftCircularCross(self):
        """d=cbrt(32Mv/(pi*sigmabmax))"""
        d=sp.cbrt((32*self.M_v)/(sp.pi*self.sigma_max))
        return d
    
    def E11_8B_DiameterShaftCircularCross(self):
        """d=2.17*cbrt(Mv/sigmabmax)"""
        d=2.17*sp.cbrt(self.M_v/self.sigma_max)
        return d
    
    def E11_9A_OutsideDiameterHollowShaft(self):
        """d_u=cbrt(32Mv/(pi*(1-k**4)*sigmabmax))"""
        d_u=sp.cbrt(32*self.M_v/(sp.pi*(1-self.k**4)*self.sigma_max))
        return d_u
    
    def E11_9B_OutsideDiameterHollowShaft(self):
        """d_u=2.17*cbrt(Mv/((1-k**4)*sigmabmax))"""
        d_u=2.17*sp.cbrt(self.M_v/((1-self.k**4)*self.sigma_max))
        return d_u
    
    def E11_10A_TorsionalMoments(self):
        """T_nom=P/(2*pi*n)"""
        print('warning: this formale is an explicit change of units, not needed if you correctly implemented units')
        T_nom=self.P/(2*sp.pi*self.n)
        return T_nom
    
    def E11_10B_TorsionalMoments(self):
        """T_nom=9550*P/n"""
        print('warning: this formale is an explicit change of units, not needed if you correctly implemented units')
        T_nom=9550*self.P/self.n
        return T_nom
    
    def E11_11A_TorsionalMoments(self):
        """T_nom=K_A*T_nom"""
        T_eq=self.K_A*self.T_nom
        return T_eq
    
    def E11_11B_TorsionalMoments(self):
        """T_eq=9550*P/n"""
        T_eq=9550*self.K_A*self.P/self.n
        return T_eq
    
    def E11_12_DesingDiameter(self):
        """dprime_a = 2.7 * cbrt(T/((1-k**4)*tau_max)"""
        dprime = 2.7 * sp.cbrt(self.T/((1-self.k**4)*self.tau_max))
        return dprime
    
    def E11_13_DesignDiameter(self):
        """dprime = 2.7 cbrt(T/tau_t)"""
        dprime = 2.7*sp.cbrt(self.T/self.tau_t)
        return dprime
    
    def E11_14_DesignDiameter(self):
        """dprime = 3.5*cbrt(M_v/sigma_b) """
        dprime = 3.5*sp.cbrt(self.M_v / self.sigma_bToel)
        return dprime
    
    def E11_15_DesignDiameter(self):
        """dprime = 4.5*cbrt(M_v/sigma_b) """
        dprime = 4.5*sp.cbrt(self.M_v / self.sigma_bToel)
        return dprime
    
    def E11_16A_DesignDiameter(self):
        """dprime = 3.4*cbrt(M_b/sigma_b) 
        using M_b, not M_v"""
        dprime = 3.4*sp.cbrt(self.M_b / self.sigma_bToel)
        return dprime

    def E11_16B_DesignDiameter(self):
        """dprime = 3.4*cbrt(M_b/sigma_b) 
        using M_v, not M_b"""
        dprime = 3.4*sp.cbrt(self.M_v / self.sigma_bToel)
        return dprime
    
    def E11_17A_DesingDiameter(self):
        """dprime_a = 3.4 * cbrt(M_b/((1-k**4)*tau_max)
        using M_b, not M_v"""
        dprime = 2.7 * sp.cbrt(self.M_b/((1-self.k**4)*self.sigma_bToel))
        return dprime
    
    def E11_17B_DesingDiameter(self):
        """dprime_a = 3.4 * cbrt(M_b/((1-k**4)*tau_max)
        using M_v, not M_b"""
        dprime = 2.7 * sp.cbrt(self.M_v/((1-self.k**4)*self.sigma_bToel))
        return dprime
    
    def E11_18A_AngleOfTorsion(self):
        """ phi= 180°/pi*l*tao_t/(r*G)"""
        phi=180*deg_/sp.pi*self.l*self.tau_t/(self.r*self.G)
        return phi
    
    def E11_18B_AngleOfTorsion(self):
        """ phi= 180°/pi*T*l/(G*I_p)"""
        phi=180*deg_/sp.pi*self.T*self.l/(self.G*self.I_p)
        return phi
    
    def E11_19A_GlobalShaftDiameter(self):
        """ d=2.32*(T)**0.25"""
        d=2.32*self.T**0.25
        return d
    
    def E11_19B_GlobalShaftDiameter(self):
        """ d=129*(K_A*P/n)**0.25"""
        d=129*(self.K_A*self.P/self.n)**0.25
        return d
    
    def E11_20_TorsionAngle(self):
        """ phi=180/pi*32*T/(pi*G)*SUM(l/d**4)"""
        phi= 0 #180/sp.pi*32*self.T/(sp.pi*self.G)*self.SUM(self.l/self.d**4)
        print('!!!!!! SUM TO BE IMPLEMENTED by user')
        return phi
    
    def E11_21_ResultingDeflection(self):
        """ f_res=sqrt(fx**2+fy**2"""
        f_res= sp.sqrt(self.f_x**2+self.f_y**2)
        return f_res
    
    def E11_22_ResultingDeflectionAngle(self):
        """ a_res=sqrt(a_x**2+a_y**2"""
        alpha_res= sp.sqrt(self.alpha_x**2+self.alpha_y**2)
        return alpha_res
    
    def E11_23_DeflectionUnderF(self):
        """ f_A=6.79*F_A/E*SUM(a_1/d)"""
        (a,i)=sp.symbols('a,i')
        f_A= 0 #6.79*self.F_A/self.E * sp.summation((a(i-1)**3-a(i)**3)/a(i)**4 , (i,1,len(self.a)))  
        # !!!!!!!!!!!!!! "*SUM(self.a_1**3/self.d_a1**3+(self.a_2**3-self.a_1**3)/self.d_a2**3+...+(self.a_n**3-self.a_n-1**3)/self.d_an**3)"    #SUM TO BE IMPLEMENTED //a_n, d_an not defined
        print('!!!!!! SUM TO BE IMPLEMENTED by user')
        return f_A
    
    def E11_24_DeflectionUnderF(self):
        """ f_A=6.79*F_A/E*SUM(a_1/d)"""
        f_B= 0
        # !!!!!!!!!!!!!! "*SUM(self.b_1**3/self.d_b1**3+(self.b_2**3-self.b_1**3)/self.d_b2**3+...(self.b_n**3-self.b_n-1**3)/self.d_bn**3)"    #SUM TO BE IMPLEMENTED //b_n, d_bn not defined
        print('!!!!!! SUM TO BE IMPLEMENTED by user')
        return f_B
    
    def E11_25_DeflectionUnderF(self):
        """ f=f_A+a/l*(f_B-f_A)"""
        f=self.f_A+self.a/self.l*(self.f_B-self.f_A)
        return f
    
    def E11_26A_AngleOfInclination(self):
        """  """
        alphaprime=10.19*self.F_A/self.E, "*SUM(self.a_1**2/self.d_a1**4+(self.a_2**2-self.a_1**2)/self.d_a2**4+...+(self.a_n**2-self.a_n-1**2)/self.d_an**4)"  #SUM TO BE IMPLEMENTED //a_n, d_an not defined
        return alphaprime
    
    def E11_26B_AngleOfInclination(self):
        """  """
        betaprime=10.19*self.F_A/self.E, "*SUM(self.b_1**2/self.d_b1**4+(self.b_2**2-self.b_1**2)/self.d_b2**4+...+(self.b_n**2-self.b_n-1**2)/self.d_bn**4)"   #SUM TO BE IMPLEMENTED //b_n, d_bn not defined
        return betaprime
    
    def E11_27A_AngleOfInclination(self):
        """ alpha=alphaprime+(f_B-f_A)/l """
        alpha=self.alphaprime+(self.f_B-self.f_A)/self.l
        return alpha
    
    def E11_27B_AngleOfInclination(self):
        """  beta=betaprime-(f_B-f_A)/l """
        beta=self.betaprime-(self.f_B-self.f_A)/self.l
        return beta
    
    def E11_28_NaturalFrequency(self):
        """ omega_k= sp.sqrt(c/m) """
        omega_k= sp.sqrt(self.c/self.m)
        return omega_k
    
    def E11_29_NaturalFrequency(self):
        """ =946*sqrt(f^-1)"""
        n_k=946*sp.sqrt(1/self.f)
        return n_k
    
    def E11_30_BendingCritSpeed(self):
        """  """
        n_kb=self.k*946*sp.sqrt(1/self.f)
        return n_kb
    
    def E11_31_ShaftAngularSpeed(self):
        """  """
        omega_k="(sp.sqrt(SUM(1/self.omega_k0**2+1/self.omega_k1**2...1*omega_kn**2))**-1"   #SUM TO BE IMPLEMENTED Omega#'s are not defined (wk0,wk1 etc.)
        return omega_k
    
    def E11_32_EigenCircleFrequency(self):
        """ omega_k=sqrt(c_t/J) """
        omega_k=sp.sqrt(self.c_t/self.J)
        return omega_k
    
    def E11_33A_RotationalVelocity(self):
        """ n_kt=30/pi*sqrt(c_t/T) """
        n_kt=30/sp.pi*sp.sqrt(self.c_t/self.T)
        return n_kt
        
    def E11_33B_RotationalVelocity(self):
        """ n_kt=72.3*sp.sqrt(T/(phi*J)) """
        n_kt=72.3*sp.sqrt(self.T/(self.phi*self.J))
        return n_kt
    
    def E11_34_EigenCircleFrequency (self):
        """ omega_k=c_t*(1/J_1+1/J_2) """
        omega_k=self.c_t*(1/self.J_1+1/self.J_2)
        return omega_k
    
    def E11_35A_TorsionCritVelocity (self):
        """ n_kt=30/pi*sqrt(c_t*(1/J_1+1/J_2)) """
        n_kt=30/sp.pi*sp.sqrt(self.c_t*(1/self.J_1+1/self.J_2))
        return n_kt
    
    def E11_35B_TorsionCritVelocity (self):
        """ n_kt=72.3*sqrt(T/self.phi*(1/J_1+1/J_2)) """
        n_kt=72.3*sp.sqrt(self.T/self.phi*(1/self.J_1+1/self.J_2))
        return n_kt
    
    
    
    
    

    





    
    
        
        
