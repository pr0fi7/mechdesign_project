# -*- coding: utf-8 -*-
"""
the objective of this file is to manage the default units provided

do not add any other methods, variables, ...
use:
    
from .Units import *


"""
import sympy as sp

UnitStrList = ['UnitStrList', 'kg_', 'm_', 'dm_', 'cm_', 'mm_', 'mu_m_', 'N_' ,'rpm_' , 's_' ,'W_' ,'Pa_', 'deg_', 'h_']

for u in UnitStrList:
    if u !='UnitStrList':
        locals()[u]=sp.symbols(u,positive=True,real=True) 
  
    
__all__ = UnitStrList



