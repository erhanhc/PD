## Example Script for ISOTROPIC PLATE UNIAXIAL TENSION CASE
from numpy import array as arr
from numpy import where
import matplotlib.pyplot as plt
from numpy import arange
from numpy import zeros,eye
from numpy import pi
from numpy import exp
from numpy.linalg import inv
from numpy import matmul
from numpy.linalg import norm
from pandas import DataFrame
from time import perf_counter
## PLOTTING FLAGS
plot_body = False

# Geometry
##--INPUTS--
L=100.#mm
unit_length = 1.#mm
ndivx=L/unit_length
thickness = unit_length#mm

#--------
W=L/2
ndivy=ndivx*(W/L)
delta = L/ndivx*3.
volume = L*W*thickness/ndivx/ndivy

xy=arr([[i,j] for i in arange(0.,L+L/ndivx,L/ndivx) for j in arange(0,W+W/ndivy,W/ndivy)],float)
if plot_body:
    fig,ax = plt.subplots()
    ax.scatter(xy[:,0],xy[:,1],color='blue')
    fig
# Neighbors 

NeighborList = [[] for i in range(len(xy))]
i=0

NumberOfBonds=0
for current in xy:
    distances = norm(current-xy,axis=1)
    neighbor_index = where((distances<=delta) * (distances!=0.0))[0]
    NumberOfBonds+=len(neighbor_index)
    NeighborList[i]=neighbor_index.tolist()
    i+=1
## Linear Peridynamic Solid Initialization
WeightedVolume = zeros(shape=len(xy),dtype=float)
for i in range(len(xy)):
    current = xy[i]
    others = xy[NeighborList[i]]
    xis = others - current
    #Gaussian = e**(-|xi|**2/d**2)
    InfluenceWeights = exp( -norm(xis,axis=1)**2/delta**2)
    #mi = Sum ( w |xi|**2 * volume)
    WeightedVolume[i] = (InfluenceWeights*volume*norm(xis,axis=1)**2).sum()

##Material Properties 
Kappa = 1.6 * 10**11 #Pascal = N/m2
Kappa = Kappa *(1/1000)**2#N/mm2
Mu = 79.3 * 10**11 #Pascal = N/m2
Mu = Mu *(1/1000)**2#N/mm2
rho = 8050#kg/m3
rho=rho*(1/1000)**3#N/mm2

##Initialize
Accel = zeros(shape=xy.shape,dtype=float)
Velocity = zeros(shape=xy.shape,dtype=float)
Deformation = zeros(shape=xy.shape,dtype=float)
ForceInternal = zeros(shape=xy.shape,dtype=float)
ForceExternal = zeros(shape=xy.shape,dtype=float)

Left_Application = where((xy[:,0]>=L-delta))[0]
Right_Application = where((xy[:,0]<=delta))[0]
Total_Applied_Force = 1.#N
ForceExternal[Left_Application,0] = Total_Applied_Force/(len(Left_Application)+1)/volume
ForceExternal[Right_Application,0] = -Total_Applied_Force/(len(Right_Application)+1)/volume


MassMatrix= eye(N=len(xy)) * volume * rho#kg
InvMassMatrix = inv(MassMatrix)#diag(1/(volume*rho))

Dilatation = zeros(shape=xy.shape[0],dtype=float)
StrainEnergyDensity = zeros(shape=xy.shape[0],dtype=float)

#Time integration

#Critical Time estimation Safety Factor 
SafetyFactor = 0.7

##CriticalTimeStep Calculation
CriticalTimeSteps = zeros(shape=xy.shape[0],dtype=float)
for i in range(len(xy)):
    current=xy[i]
    others = xy[NeighborList[i]]
    xis = others-current
    CpEffs = 18 * Kappa / norm(xis,axis=1) / pi / delta**4 #N/mm2/mm/mm4
    CriticalTimeSteps[i]=(2*rho / (CpEffs*volume).sum()/1000)**0.5
DeltaTime=CriticalTimeSteps.min()*SafetyFactor#seconds

#-----
Time=0.0
n=0
while Time<=10.:
    #Updating Times
    Time1Step = Time+DeltaTime
    TimeHalfStep = 0.5*(Time1Step+Time)
    #First Partial Velocity Update
    VelocityHalfStep = Velocity + Accel*(TimeHalfStep - Time)
    #Updating Nodal Disp
    Deformation1Step = Deformation + VelocityHalfStep * DeltaTime
    
    ## Linear Peridynamic Solid Internal Force Calculation 

    ForceInternal[:,:]=0.
    Dilatation[:] = 0.
    StrainEnergyDensity[:]=0.

    #Dilatation Calculation
    for i in range(len(xy)):
        current = xy[i]
        others = xy[NeighborList[i]]
        xis = others - current
        InfluenceWeights = exp(-norm(xis,axis=1)**2/delta**2)
        etas = Deformation1Step[NeighborList[i],:] - Deformation1Step[i,:]
        ys = xis+etas
        BondExtensionScalars = norm(ys,axis=1)- norm(xis,axis=1)
        Dilatation[i] = 3 * (InfluenceWeights * norm(xis,axis=1) * BondExtensionScalars).sum() * volume / WeightedVolume[i]
    StrainEnergyDensity[i] += Kappa /2 * Dilatation[i]**2
    for i in range(len(xy)):
        current = xy[i]
        others = xy[NeighborList[i]]
        xis = others - current
        InfluenceWeights = exp(-norm(xis,axis=1)**2/delta**2)
        etas = Deformation1Step[NeighborList[i],:] - Deformation1Step[i,:]
        ys = xis+etas
       
        BondExtensionScalars = norm(ys,axis=1) - norm(xis,axis=1)
        BondExtensionScalarsIsotropic = Dilatation[i] * norm(xis,axis=1) / 3
        BondExtensionScalarsDeviatoric = BondExtensionScalars - BondExtensionScalarsIsotropic        
        
        InternalForceScalarsIsotropic = 3 * Kappa * Dilatation[i] * InfluenceWeights * norm(xis,axis=1) / WeightedVolume[i]
        InternalForceScalarsDeviatoric = (15 * Mu / WeightedVolume[i]) * InfluenceWeights * BondExtensionScalarsDeviatoric
        
        InternalForceScalars = InternalForceScalarsIsotropic + InternalForceScalarsDeviatoric
        InternalForceScalars = InternalForceScalars.reshape((InternalForceScalars.shape[0],1))
        
        DeformationDirectionVectors = ys/norm(ys,axis=1).reshape((norm(ys,axis=1).shape[0],1))
        
        ForceInternal[i,:]+= ( InternalForceScalars * DeformationDirectionVectors * volume).sum(axis=0)
        ForceInternal[NeighborList[i],:]-= InternalForceScalars * DeformationDirectionVectors * volume            
        
        StrainEnergyDensity[i] += (15 * Mu / WeightedVolume[i] ) / 2 * (InfluenceWeights * BondExtensionScalarsDeviatoric**2 * volume).sum()
    #----------
    #Computing Acceleration
    Accel = matmul(InvMassMatrix,ForceInternal*volume+ForceExternal*volume)
    Accel = Accel/1000
    #Second Partial Velocity Update
    Velocity = VelocityHalfStep + (Time1Step-TimeHalfStep) * Accel
    Deformation = Deformation1Step
    
    WorkExternals = (ForceExternal*Deformation).sum(axis=1)
    StrainEnergys = StrainEnergyDensity*volume
    PotentialEnergys = (StrainEnergys - WorkExternals) 
    KineticEnergys = matmul(MassMatrix,Velocity**2).sum(axis=1)/2
    KineticEnergys = KineticEnergys/1000
    Lagrangians = KineticEnergys-PotentialEnergys
    Time = Time1Step
    y = xy+Deformation
    print(y[where((xy[:,0]==50.)*(xy[:,1]==25.)),0])
    print(Lagrangians.sum())
    print(Time)
    # if not bool(n%10000):
    #     df = DataFrame()
    #     df[['X','Y']] = xy
    #     df[['DispX','DispY']] = Deformation
    #     df[['VelX','VelY']] = Velocity
    #     df[['AccelX','AccelY']] = Accel
    #     df[['InternalForceX','InternalForceY']] = ForceInternal
    #     df[['ExternalForceX','ExternalForceY']] = ForceExternal
    #     df['StrainEnergy'] = StrainEnergys
    #     df['WorkPotential'] = WorkExternals
    #     df['PotentialEnergy'] = PotentialEnergys
    #     df['KineticEnergy'] = KineticEnergys
    #     df['Lagrangian'] = Lagrangians
    #     df.to_csv('out%d.csv'%n,index_label='ID',float_format='%1.5E')
    # with open('pd_out.txt','a') as file:
    #     file.writelines('Time:%1.7f, Iter:%d, TotalKineticEnergy:%1.5E, TotalPotentialEnergy:%1.5E, Lagrangian:%1.5E,\n' %(Time,n,KineticEnergys.sum(),PotentialEnergys.sum(),Lagrangians.sum()))
    n+=1