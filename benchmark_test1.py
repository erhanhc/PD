# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 11:14:35 2021

@author: Aaron
"""

from pd_funcs import preprocess, iterate

import pandas, numpy
columns = ['coordx','coordy','volume','dispx','dispy','SED','Dilatation','D1','D2','S1','S2','pforcex','pforcey','bforcex','bforcey','massvecx','massvecy','pforceoldx','pforceoldy','velhalfoldx','velhalfoldy','velx','vely','velhalfx','velhalfy']
ndivx=100
ndivy=50
totalpd = ndivx*ndivy
#Collocation points with all attributes equal to 0. with datatype float32
df = pandas.DataFrame(numpy.zeros((totalpd,len(columns)),dtype='float32'),columns=columns)
#Geometry
dt = 1.
L=1.
W=0.5
delta=L/ndivx
horizon = 3.015 * delta
thickness=delta
current_id=0
for i in range(0,ndivx):
    for j in range(0,ndivy):
        df.loc[current_id].coordx=-0.5*L+(delta/2)+(i)*delta
        df.loc[current_id].coordy=-0.5*W+(delta/2)+(j)*delta
        df.loc[current_id].volume = delta**2 * thickness
        current_id=current_id+1
#Material Properties
E = 200
nu = 1/3
kappa = E/3/(1-2*nu)
mu = E/2/(1+nu)
a = 0.5 * (kappa - 2*mu)
b = 6*mu/numpy.pi/thickness/horizon**4
d = 2 / numpy.pi/thickness/horizon**3
#Preprocess for SCF's 
condition_list = ['uniaxial stretch x','uniaxial stretch y','simple shear in x-y']
disp_grad=0.001
preprocess(df,condition_list,disp_grad,horizon,delta,a,b,d,mu)
#Iteration with ADR
iterate(df,'uniaxial tensile loading',1,500,E*1e-3,horizon,delta,a,b,d)
