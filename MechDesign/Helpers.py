# -*- coding: utf-8 -*-
"""
this file holds helper classes and function, not intended to be viewed by end-users
"""

import sympy as sp
import numpy as np
from IPython.display import display, Math


#from .Units.Units import *
from  .Units import UnitMethods as UM

    
# -----------------------------------------------------------------------------
class MySymbol(sp.Dummy):
    """
    subclass of standard sympy symbols, but adding description
    the symbols will all start with a leading '_', 
        using 'print' will show this
        using 'display' in a notebook will hide the leading '_'
    
    
    example
    =======
    # make new objects with description
    x = MySymbol('x')
    x.description = 'Distance (m)'
    t = MySymbol('t', 'Time (s)')
    print( x.description, t.description)
    
    # test
    expr = (x*t + 2*t)/t
    print (simplify(expr))  
    
    """
    # limit the attributes of this type of symbols, in order to avoid typo's
    __slots__=['description','def_value','comment','range']
    def __new__(cls, name, description=''):
        obj = super().__new__(cls,name,canonical=False,real=True)
        return obj
    def __init__(self,name,description=''):
        # self.symbol = sp.symbol( name)
        self.description = description
        self.comment = None
        self.def_value= None
        self.range = None



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
        key = MySymbol(key) if type(key) == str else key
        listOfParameters.append(key.name)
        listOfValues.append(value)    
    # remove units from expression
    MyExpression = UM.RemoveUnits(MyExpression)
    # create a list of effectively used symbols in MyExpresion, a lis fixes the order
    listEffectiveSymbols = list(MyExpression.free_symbols) #can be different order!!
    listNamesEffectiveSymbols = [s.name for s in listEffectiveSymbols]
    listEffectiveSymbols = [listEffectiveSymbols[listNamesEffectiveSymbols.index(n)] for n in listOfParameters]
    
    # Create lambda expression
    MyLamExp = sp.lambdify(listEffectiveSymbols,MyExpression,"numpy")
    

    # evaluate the expression, 
    if len(listOfValues) == 1:
        A = np.array(listOfValues).reshape(-1)
        values = MyLamExp(A)
        return values
    elif len(listOfValues) == 2:
        # A, B = np.meshgrid(np.array(listOfValues[0]), np.array(listOfValues[1]), indexing="ij")
        A, B = np.meshgrid(np.array(inputs[listEffectiveSymbols[0].name]).reshape(-1),\
                           np.array(inputs[listEffectiveSymbols[1].name]).reshape(-1),\
                               indexing="ij")
        values = MyLamExp(A, B)
        return values, A, B
    elif len(listOfValues) == 3:
        #A, B, C = np.meshgrid(np.array(listOfValues[0]).reshape(-1), np.array(listOfValues[1]).reshape(-1), np.array(listOfValues[2]).reshape(-1), indexing="ij")
        A, B, C = np.meshgrid(np.array(inputs[listEffectiveSymbols[0].name]).reshape(-1),\
                              np.array(inputs[listEffectiveSymbols[1].name]).reshape(-1),\
                              np.array(inputs[listEffectiveSymbols[2].name]).reshape(-1),\
                               indexing="ij")
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
    # for t in kwargs:
    #     print('kwarg '+t)
    #     print(kwargs[t])
    for t in MyExp.free_symbols:
        # print(t.name)
        if t.name in kwargs:  # if a symbol's name can be found in the **kwargs, prioritise this over the ElInstance value
            # print('in kwargs')    
            subst[t] = kwargs[t.name]
        elif hasattr(ElInstance, t.name): # else if the symbol can be found in the ElInstance 
            value = getattr(ElInstance, t.name)
            subst[t] = value
    # print(subst)
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
    

def set_all_to_positive(Expression):
    """
    all symbols in the expression will be set to positive, 
    this can help when simplification of an expressions fails
    use with caution

    """
    temp, temp2 = sp.posify(Expression)
    temp = temp.subs(temp2)
    return temp


    

# -----------------------------------------------------------------------------
def EqPrint(variable, MyExpression):
    """ 
    display an expression as an equation
        variable:  text string of the left hand side
        MyExpression: symbolic expression
    """
    if isinstance(MyExpression,float):
        f = float(f'{float(MyExpression):.4g}')
    elif isinstance(MyExpression,int):
        f = MyExpression
    else:
        f = MyExpression.xreplace({n :float(f'{float(n):.4g}') for n in MyExpression.atoms(sp.Number) if int(float(n))!=n })
    
    Equation = sp.Eq(sp.Symbol(variable), f)
    display(Equation)
    return Equation
    
    
def MyHelp(Expression):
    """
    Print out the Expression and the definition of all it's variables
    Expression: symbolic expression
    """    
    # print the expression
    f = Expression.xreplace({n:float(f'{float(n):.4g}') for n in Expression.atoms(sp.Number) if int(float(n))!=n})
    display(f)
    # l = sp.latex(Expression)
    # display(Math(l))
    for t in Expression.free_symbols: 
        if not t.name[-1]=='_':  # if it is not a symbol
            try:
                print(('    '+t.name+' : '+t.description))
                if not t.comment==None: print('    '+' '*len(t.name)+'     '+t.comment)
                if not t.def_value==None: print('    '+' '*len(t.name)+'     '+'default value: '+str(t.def_value))
                if not t.range==None: print('    '+' '*len(t.name)+'     '+'range: '+';'.join(str(s) for s in t.range))
            except: print('this is not a correctly created "MySymbol" symbol')





