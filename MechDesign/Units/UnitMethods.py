# -*- coding: utf-8 -*-
"""
this file contains several methods, helpfull in managing units
per definition units have trailing '_'
"""

import sympy as sp
from .Units import *



# -----------------------------------------------------------------------------

def PrintAvailableUnits():
    list2print = UnitStrList[1::]
    print(list2print)


def All_to_SI(MyExpression):
    """
    set: 
       mm,cm,dm,mu_m to m; 
       rpm to rad/s; 
       W to kg*m**2/s**3, 
       Pa to N/mm**2 
       h to s
    """
    All_to_SI = { mu_m_   :    1e-6*m_\
                 ,mm_     :    1e-3*m_\
                 ,cm_     :    1e-2*m_\
                 ,dm_     :    1e-1*m_\
                 ,rpm_    :    1/60*2*sp.pi/s_\
                 ,W_      :    kg_*m_**2/s_**3\
                 ,Pa_     :    N_/mm_**2\
                 ,h_      :    3600*s_\
                 ,deg_    :    sp.pi/180\
                     }
            
    MyExp = MyExpression.evalf(subs=All_to_SI)
    return MyExp


def RemoveUnits(MyExpression):
    """ replace all unit symbols with 1 (units = symbols ending in '_' eg N_) """
    listOfUnitsUsed = {}
    for tt in MyExpression.free_symbols:
        if tt.name[-1] == "_":
            listOfUnitsUsed[tt] = 1
    MyExp = MyExpression.evalf(subs=listOfUnitsUsed)
    return MyExp

# def SetPower_to_int(MyExpression):
#     """ when a floating a power of '1.0' is found it is convert to int """
#     #f = MyExpression.xreplace({n:int(n) for n in MyExpression.atoms(sp.Number)})
#     display(MyExpression)
#     e = MyExpression.replace(
#         lambda x: x.is_Float and x==int(x),
#         lambda x: reps.setdfault(x,Dummy()))

#     return e


def N_mm2_to_Pa(MyExpression):
    """ Convert N_/mm_**2 to Pa """
    Npmm2_to_Pa = {N_/mm_**2 : Pa_}
    MyExp = MyExpression.evalf(subs=Npmm2_to_Pa)
    return MyExp


def kgm_s2_to_N(MyExpression):
    """ convert kg*m/s**2 to N """
    kgmps2_to_N = {kg_*m_/s_**2 : N_}
    MyExp = MyExpression.evalf(subs=kgmps2_to_N)
    return MyExp


def m_to_mm(MyExpression):
    """ convert m to mm """
    m_to_mm = {m_ : 1e3*mm_}
    MyExp = MyExpression.evalf(subs=m_to_mm)
    return MyExp

def m_to_mu_m(MyExpression):
    """ convert m to mu_mm """
    m_to_mu_m = {m_ : 1e6*mu_m_}
    MyExp = MyExpression.evalf(subs=m_to_mu_m)
    return MyExp

def mm_to_mu_m(MyExpression):
    """ convert mm to mu_mm """
    mm_to_mu_m = {mm_ : 1e3*mu_m_}
    MyExp = MyExpression.evalf(subs=mm_to_mu_m )
    return MyExp


def rpm_to_rad_s(MyExpression):   
    """ convert rpm to rad/s """
    rpm_to_radps = {rpm_ : 1/60*2*sp.pi/s_}
    MyExp = MyExpression.evalf(subs=rpm_to_radps)
    return MyExp


def N_to_kgm_s2(MyExpression):   
    """ convert N to kg*m/s**2 """
    N_to_kgmps2 = {N_ : kg_*m_/s_**2}
    MyExp = MyExpression.evalf(subs=N_to_kgmps2)
    return MyExp 

