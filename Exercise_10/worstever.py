# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from math import *
from scipy import *
import numpy

x=numpy.zeros(5)
y=numpy.zeros(5)
x[0]=0
x[1]=1.8
x[2]=2.35
x[3]=3.0
x[4]=5.0
y[0]=0.01
y[1]=0.15
y[2]=2.5
y[3]=0.1
y[4]=0.002

x=0.01
z0=((y1-y0)/(x1-x0))*(x-x0)+y0
print z0
p=0.078/((x-2)**4+sin(x-3)**8)
print p

x=1.7
z0=((y1-y0)/(x1-x0))*(x-x0)+y0
print z0
p=0.078/((x-2)**4+sin(x-3)**8)
print p

x=1.81
z1=((y2-y1)/(x2-x1))*(x-x1)+y1
print z1
p=0.078/((x-2)**4+sin(x-3)**8)
print p

x=2.34
z1=((y2-y1)/(x2-x1))*(x-x1)+y1
print z1
p=0.078/((x-2)**4+sin(x-3)**8)
print p

x=2.36
z2=((y3-y2)/(x3-x2))*(x-x2)+y2
print z2
p=0.078/((x-2)**4+sin(x-3)**8)
print p

x=2.99
z2=((y3-y2)/(x3-x2))*(x-x2)+y2
print z2
p=0.078/((x-2)**4+sin(x-3)**8)
print p

x=3.01
z3=((y4-y3)/(x4-x3))*(x-x3)+y3
print z3
p=0.078/((x-2)**4+sin(x-3)**8)
print p

x=4.99
z3=((y4-y3)/(x4-x3))*(x-x3)+y3
print z3
p=0.078/((x-2)**4+sin(x-3)**8)
print p