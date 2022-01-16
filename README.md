# PD Modelling of Composite Laminates for Damage Analysis
**Currently implementing the theory given in Peridynamic Theory and Its Application by Madenci and Oterkus** 

***Whats been done so far?***
1. Ordinary state based PD formulation is studied following the given instructions/codes given in the refence book. 
2. ADR process is applied but not verified for quasi-static solutions. 
3. Time integeration is added. Under verification. 
4. Pandas-HDF5 object and data formats are implemented using Ordinary state based formulation and isotropic material model. 

***TO-DO List***
1. Reformulation of the theory given in the reference book and its references. 
2. Verification and clear formulation of ADR process. Problems still exist in certain terms. ***-> This requires computation of stiffness terms Kij using small displacement assumption. Volume corr. factor for the current material point is taken as =1.***
3. Time integration schemes to be implemented in the code. ***-> Forward and Backward differencing scheme for velocity and displacement fields.***
4. Composite lamina and lamination theory to be studied, formulated and implemented in the code.
5. Damage prediction chapter to be studied, formulated and implemented in the code. 

***Files in the Repository***
1. pd_funcs.py: Consists of predefined functions that are used for numerical operations that are necessary for preprocess, solution and outputs. 
2. benchmark_test1.py. Tension test of a isotropic plate using ADR. (Not verified yet)
3. Docs: Reformulations that are given in the reference book and schemas for the generated codes. 
3. benchmark_test2.py: Tension test of a isotropic plate using time integration.
EHC.