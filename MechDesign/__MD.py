# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 06:19:23 2022
@author: u0032080
TODO 
    create License file
    remove units from files, include units for specific expressions with constants
    

This folder will contain - for student's reference - expressions from Roloff & Matek. 
These expressions are only for internal use and e.g. for studying/answering the exam. 
These expressions may never and under no exception be distributed or shared to anyone else 
and serve only for the purpose of students that follow the KULeuven courses on dimensioning of machine parts 
By using these files you commit to adhering to this condition ! 
    
"""


import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display

from .RnM.C2_Tolerance import Tolerance as Tolerance
from .RnM.C3_StrengthAndStress import StrengthAndStress as StrengthAndStress  # still incomplete
from .RnM.C11_Shaft import Shaft as Shaft
from .RnM.C12_ShaftConnection import ShaftConnection as ShaftConnection
from .RnM.C14_BallBearing import BallBearing as BallBearing
from .RnM.C16_Belt import Belt as Belt
from .RnM.C20_GearDesign import GearDesign as GearDesign
from .RnM.C21_InvoluteGear import InvoluteGear as InvoluteGear


# creating units for symbols
( kg_, m_, dm_, cm_, mm_, mu_m_, N_, rpm_ , s_, W_, Pa_)= sp.symbols(\
 'kg_, m_, dm_, cm_, mm_, mu_m_, N_ ,rpm_ , s_ ,W_ ,Pa_',positive=True)



# -----------------------------------------------------------------------------
"convert units"


def All_to_SI(MyExpression):
    """
    set: 
       mm to m; 
       rpm to rad/s; 
       W to kg*m**2/s**3, 
       Pa to N/mm**2 
    """
    All_to_SI = { mu_m_   :    1e-6*m_\
                 ,mm_     :    1e-3*m_\
                 ,cm_     :    1e-2*m_\
                 ,dm_     :    1e-1*m_\
                 ,rpm_    :    1/60*2*sp.pi/s_\
                 ,W_      :    kg_*m_**2/s_**3\
                 ,Pa_    :     N_/mm_**2}
            
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


def Npmm2_to_Pa(MyExpression):
    """ Convert N_/mm_**2 to Pa """
    Npmm2_to_Pa = {N_/mm_**2 : Pa_}
    MyExp = MyExpression.evalf(subs=Npmm2_to_Pa)
    return MyExp


def kgmps2_to_N(MyExpression):
    """ convert kg*m/s**2 to N """
    kgmps2_to_N = {kg_*m_/s_**2 : N_}
    MyExp = MyExpression.evalf(subs=kgmps2_to_N)
    return MyExp


def m_to_mm(MyExpression):
    """ convert m to mm """
    m_to_mm = {m_ : 1e3*mm_}
    MyExp = MyExpression.evalf(subs=m_to_mm)
    return MyExp


def rpm_to_radps(MyExpression):   
    """ convert rpm to rad/s """
    rpm_to_radps = {rpm_ : 1/60*2*sp.pi/s_}
    MyExp = MyExpression.evalf(subs=rpm_to_radps)
    return MyExp


def N_to_kgmps2(MyExpression):   
    """ convert N to kg*m/s**2 """
    N_to_kgmps2 = {N_ : kg_*m_/s_**2}
    MyExp = MyExpression.evalf(subs=N_to_kgmps2)
    return MyExp


def EqPrint(variable, MyExpression):
    """ 
    display an expression as an equation
        variable:  text string of the left hand side
        MyExpression: symbolic expression
    """
    Equation = sp.Eq(sp.Symbol(variable), MyExpression)
    display(Equation)
    return Equation



#-----------------------------------------------------------------------------
def evaluate(MyExpression, inputs):
    ''' 
    evaluate a symbolic expression at set points, using lambda functions
    MyExpression:   symbolic expression
                    can still hold unit symbols, 
                    ALL UNIT SYMBOLS WILL BE SET TO '1'
          inputs:   dictionary with key, value pairs, max 3 pairs
             key:   holds the symbol of to be replaced     
                    alternatively, it can hold the corresponding string to be converted to symbol
          values:   single value, list, or array
                    CANNOT HOLD UNIT SYMBOLS, 
                    double check that units are corresponding to the units of MyExpression
    '''
    # step 1 get list of args
    listOfParameters = []
    listOfValues = []
    for key, value in inputs.items():
        key = sp.Symbol(key) if type(key) == str else key
        listOfParameters.append(key)
        listOfValues.append(value)    
    # remove units from expression
    MyExpression = RemoveUnits(MyExpression)
    # lambdify the expression
    MyLamExp = sp.lambdify(listOfParameters, MyExpression, "numpy")
    # evaluate the expression, 
    if len(listOfValues) == 1:
        A = np.array(listOfValues).reshape(-1)
        values = MyLamExp(A)
        return values
    elif len(listOfValues) == 2:
        A, B = np.meshgrid(np.array(listOfValues[0]), np.array(listOfValues[1]), indexing="ij")
        values = MyLamExp(A, B)
        return values, A, B
    elif len(listOfValues) == 3:
        A, B, C = np.meshgrid(np.array(listOfValues[0]).reshape(-1), np.array(listOfValues[1]).reshape(-1), np.array(listOfValues[2]).reshape(-1), indexing="ij")
        values = MyLamExp(A, B, C)
        return values, A, B, C


def substitute(MyExp,ElInstance,**kwargs):
    ''' 
    susbstitutes symbols with the values from a specific instance of a design component 
    unknown symbols will be ignored
    MyExp:         Symbolic expression
    ElInstance:    Instance of a mechanical design, holding specifc values/symbols/expressions
    pairs of symbol and values can also be given, they take priority over the ElInstance values
    '''
    # get a list of symbols used in the expression
    subst = {} # empty dictionary that will hold the 'symbol:value' pairs.
    for t in kwargs:
        print('kwarg '+t)
        print(kwargs[t])
    for t in MyExp.free_symbols:
        print(t.name)
        if t.name in kwargs:  # if a symbol's name can be found in the **kwargs, prioritise this over the ElInstance value
            print('in kwargs')    
            subst[t] = kwargs[t.name]
        elif hasattr(ElInstance, t.name): # else if the symbol can be found in the ElInstance 
            value = getattr(ElInstance, t.name)
            subst[t] = value
    print(subst)
    MyNewExp = MyExp.evalf(subs = subst)
    return MyNewExp


def reorderEquation(Lhand, ExpIn, Var2Ext):
    """
    Lhand:      a text with the sympbol on the left hand side, written as a text e.g. 'F_t'
    ExpIn:      the symbolic expression for the value on the left had side
    Var2Ext:    a ext with the name of the symbol to be extracted
    """
    # from the expression, create an equation
    Eq2Use = sp.Eq(sp.symbols(Lhand),ExpIn)
    # solve the equation for a given symbol
    Exp_result = sp.solve(Eq2Use,sp.symbols(Var2Ext))
    if len(Exp_result) != 1:
        raise Warning("!!! multiple symbolic solutions found !!! " + str(len(Exp_result)) + " answers found")
    else: 
        return Exp_result[0] 
  
    
def myHelp(Expression):
    """
    Print out the Expression and the definition of all it's variables
    Expression: symbolic expression
    """    
    # print the expression
    display(Expression)     
    for t in Expression.free_symbols: print('    '+t.name+' : '+t.description)    


def set_all_to_positive(Expression):
    """
    all symbols in the expression will be set to positive, 
    this can help when simplification of an expressions fails
    use with caution

    """
    temp, temp2 = sp.posify(Expression)
    temp = temp.subs(temp2)
    return temp
