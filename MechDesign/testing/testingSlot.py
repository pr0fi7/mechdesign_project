# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 12:08:34 2023

@author: u0032080
"""

class testingSlot():
    __slots__ = ['a','b']
    
    

rr = testingSlot()
rr.a = 5
print(rr.a)
rr.c = 8