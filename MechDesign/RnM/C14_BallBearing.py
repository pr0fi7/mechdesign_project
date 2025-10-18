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


class BallBearing:
    
    
    MySymbolDict={  'C'     :'[kN]requireddynamicloadcapacity',\
                    'P'     :'[kN]dynamicload',\
                    'p'     :'[]lifeexpextancyexponent',\
                    'f_L'   :'[]standarddynamicloadvalue{figure14-35}',\
                    'f_n'   :'[]speedfactor{tabel14-4',\
                    'n'     :'[rpm]rotationalspeed',\
                    'L_10'  :'[1E6rotatons]nominallifeexpectency',\
                    'L_10h' :'[h]nominallifeexpectancy',\
                    'P_0'   :'[kN]equivalentstaticload',\
                    'C_0'   :'[kN]staticloadbearingcapacity',\
                    'S_0'   :'[]staticloadbearingsafetyfactor',\
                    'X_0'   :'[]staticradialloadfactor{tabel14-3b}',\
                    'F_r0'  :'[kN]staticradialload',\
                    'F_a0'  :'[kN]staticaxialload',\
                    'Y_0'   :'[]staticaxialloadfactor{tabel14-3b}',\
                    'X'     :'[]radialloadfactor{tabels14-3aand14-2}',\
                    'F_r'   :'[kN]radialload',\
                    'Y'     :'[]axialloadfactor{tabels14-3aand14-2}',\
                    'F_a'   :'[kN]axialload',\
                    'P_list':'[kN]thisisaLISTofnumericalvaluesofdynamicload',\
                    'n_list':'[rpm]thisisaLISTofnumericalvaluesofrotationalspeeds',\
                    'q_list':'[%]thisisaLISTofnumericalvaluesof%thespecificloadapplies',\
                    'n_m'   :'[rpm]thisistheaveragerpm',\
                    'P_min' :'[kN]minimumdynamicloadincaseoflinearloadchange',\
                    'P_max' :'[kN]maximumdynamicloadincaseoflinearloadchange',\
                    'L_na'  :'[1e6revolutions]modifiedlifeexpectancy',\
                    'a_1'   :'[]breakdownprobabilityfactor{figure13-19B}',\
                    'a_ISO' :'[]useconditionfactor{tabel14-12}',\
                    'L_nah' :'[h]modifiedlifeexpectancy',\
                    'L_nahList':'[h]thisisaLISTofnumericalvaluesofthel_nah',\
                    'C_100' :'[kN]dynamicloadfactorwithlifeexpectancyof100km',\
                    's_slag':'[m]strokelength',\
                    'n_slag':'[min-1]strokefrequency(dubbelstrokespermin)',\
                    'v_m'   :'[m/min]averagespeed',\
                    'v_list':'[m/min]thisisaLISTofnumericalvaluesofaveragespeeds',\
                    'F_x'   :'[kN]loadonthelinearguiderail,xdirection',\
                    'F_y'   :'[kN]loadonthelinearguiderail,xdirection',\
                    }


    # limit the possible attributes, avoiding typo's from students
    __slots__ = [str(k) for k in MySymbolDict.keys()]

    # initialise the symbols
    def __init__(self):
        for k in self.MySymbolDict.keys():
            setattr(self, str(k), MySymbol(str(k), self.MySymbolDict[k]))
            
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------


#
                               
    def E14_1A_DynamicLoad(self):
        """ C = P * F_l / f_n """
        C = self.P *self.f_L / self.f_n
        return C
    
    def E14_1B_DynamicLoad(self):
        """ C = P *(60*n*L_10h/1e6)**(1/p) """
        C = self.P * (60*self.n*self.L_10h/1e6)**(1/self.p)
        return C
    
    def E14_2_StaticLoad(self):
        """ C_0 = P_0 * S_0 """
        C_0 = self.P_0 * self.S_0 
        return C_0
    
    def E14_3_StaticSafetyFactor(self):
        """ S_0 = C_0 / P_0 """
        S_0 = self.C_0 / self.P_0 
        return S_0
    
    def E14_4_EquivalentStaticLoad(self):
        """ P_0 = X_0 * F_r0 + Y_0 * F_a0 """
        P_0 = self.X_0 * self.F_r0 + self.Y_0 * self.F_a0
        return P_0
    
    def E14_5A1_NominalLifeExpetancy(self):
        """ L_10 = (C / P)**p """
        L_10 = (self.C / self.P)**self.p
        return L_10
    def E14_5A2_NominalLifeExpectancy(self):
        """ L_10h = 1e6*L_10/(60*n) """
        L_10h = 1e6*self.L_10/(60*self.n)
        return L_10h
    
    def E14_5B_standardDynamicLoadValue(self):
        """ f_L = C/ P * f_n """
        f_L = self.C / self.P * self.f_n
        return f_L
    
    def E14_6_EquivalentDynamicLoad(self):
        """ P = X * F_r + Y * F_a """
        P = self.X * self.F_r +  self.Y * self.F_a 
        return P
    
    def E14_7_CyclicLoad(self):
        """ P = ( sum(Pn**p*n/nm*qn/100))**(1/p)
        this function will not work when using symbolic expressions,
        a for loop will loop trhough the lists
        """
        P = 0 
        if self.P_list is list and self.n_list is list and self.q_list is list:
            for P_i, n_i, q_i in zip(self.P_list, self.n_list, self.q_list):
                P += (P_i**self.p) * n_i / self.n_m * q_i / 100
            P = P**(1/self.p)
        return P
    
    def E14_8_AverageSpeed(self):
        """ n_= sum (ni * qi /100)
        this function will not work when using symbolic expressions, 
        a for loop will loop through the lists
        """
        n_m = 0
        if self.n_list is list and self.q_list is list:
            for n_i, q_i in zip(self.n_list, self.q_list):
                n_m += n_i * q_i / 100
        return n_m
    
    def E14_9_CyclicLoad(self):
        """ P = ( sum(Pn**p*qn/100))**(1/p)
        this function will not work when using symbolic expressions,
        a for loop will loop through the lists
        """
        P = 0 
        if self.P_list is list and self.q_list is list:
            for P_i, q_i in zip(self.P_list, self.q_list):
                P += (P_i**self.p ) * q_i / 100
            P = P**(1/self.p)
        return P
    
    def E14_11A_ModifiedLifeExpectancy(self):
        """ L_na = a_1 * a_ISO * L_10 """
        L_na = self.a_1 * self.a_ISO * self.L_10
        return L_na
    
    def E14_11B_ModifiedLifeExpectancy(self):
        """ L_nah = a_1 * a_ISO * L_10h """
        L_nah = self.a_1 * self.a_ISO * self.L_10h
        return L_nah
    
    def E14_12_CyclicLoadL_nah(self):
        """ L_nah = 100 / (sum (qi/Lnahi))
        this function will not work when using symbolic expressions,
        a for loop will loop through the lists
        """
        L_nah = 0
        if self.L_nahList is list:
            for q_i, L_nah_i in zip(self.q_list, self.L_nahList):
                L_nah += q_i/L_nah_i
            L_nah = 100 / L_nah
        return L_nah
    
    def E14_13A_LinearSlideNominalLifeExpectancy(self):
        """ L_10 = (C_100 / P)**p * 1E5 """
        L_10 = (self.C_100 / self.P )**self.P * 1e5
        return L_10
    
    def E14_13B1_LinearSlideNominaLifeExpectancy(self):
        """ L_10h = L_10 / ( 2 * s_slag * n_slag * 60) """
        L_10h = self.L_10 / (2 * self.s_slag * self.n_slag * 60)
        return L_10h
    
    def E14_13B2_LinearSlideNominaLifeExpectancy(self):
        """ L_10h = L_10 / ( 60 * v_m) """
        L_10h = self.L_10 / (60 * self.v_m)
        return L_10h
    
    def E14_14_LinearSlideEquivalentLoad(self):
        """ P = abs(F_x) + abs(F_y) """
        P = sp.Abs(self.F_x) + sp.Abs(self.F_y)
        return P
        
    def E14_15_LinearSlideDynamicLoad(self):
        """ P = (sum (abs(P)**1/p * vi/v_m *q_i / 100))
        this function will not work when using symbolic expressions,
        a for loop will loop through the lists
        """
        P = 0
        if self.v_list is list and self.q_list is list:
            for P_i, v_i, q_i in zip(self.P_list, self.v_list, self.q_list):
                P += ((sp.Abs(P_i))**self.p) * sp.Abs(v_i) / self.v_m * q_i / 100 
            P = P**(1/self.p)
        return P
    
    def E14_16_LinearSlideAverageSpeed(self):
        """ v_m  = sum (abs(v_i)*q_1/100)
        this function will not work when using symbolic expressions,
        a for loop will loop through the lists
        """        
        v_m = 0
        if self.v_list is list and self.q_list is list:
            for v_i, q_i in zip(self.v_list, self.q_list):
                v_m += sp.Abs(v_i) * q_i / 100
        return v_m
    
    def E14_17_LinearSlideSafetyFactor(self):
        """ S_0 = C_0 / P_0 """
        S_0 = self.C_0 / self.P_0 
        return S_0
        