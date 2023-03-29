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
## PLOTTING FLAGS
plot_body = False

# Geometry
##--INPUTS--
L=100.
ndivx=100./1.
thickness = 1.#mm

#--------
W=L/2
ndivy=ndivx/2
delta = L/ndivx*3.
volume = L/ndivx*W/ndivy*thickness

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
    distances = (current - xy)**2
    distances = distances[:,0] + distances[:,1]
    distances = distances**0.5
    neighbor_index = where((distances<=delta) * (distances!=0.0))[0]
    NumberOfBonds+=len(neighbor_index)
    NeighborList[i]=neighbor_index.tolist()
    i+=1

## Linear Peridynamic Solid Initialization
WeightedVolume = zeros(shape=len(xy),dtype=float)
for i in range(len(xy)):
    current = xy[i]
    for j in NeighborList[i]:
        other = xy[j]
        xi = other - current           
        InfluenceWeight= exp(-norm(xi)**2/delta**2)
        WeightedVolume[i]+=InfluenceWeight*norm(xi)**2*volume

##Material Properties 
Kappa = 1.6 * 10**11 #Pascal = N/m2
Kappa = Kappa *(1/1000)**2
Mu = 79.3 * 10**9 #Pascal = N/m2
Mu = Mu *(1/1000)**2
rho = 8050#kg/m3
rho=rho*(1/1000)**3

##Initialize
Accel = zeros(shape=xy.shape,dtype=float)
Velocity = zeros(shape=xy.shape,dtype=float)
Deformation = zeros(shape=xy.shape,dtype=float)
ForceInternal = zeros(shape=xy.shape,dtype=float)
ForceExternal = zeros(shape=xy.shape,dtype=float)

Left_Application = where((xy[:,0]>=100.-delta))[0]
Right_Application = where((xy[:,0]<=0.+delta))[0]
Total_Applied_Force = 1.
ForceExternal[Left_Application,0] = Total_Applied_Force/(len(Left_Application)+1)
ForceExternal[Right_Application,0] = -Total_Applied_Force/(len(Right_Application)+1)


MassMatrix= eye(N=len(xy)) * volume * rho
InvMassMatrix = inv(MassMatrix)

Dilatation = zeros(shape=xy.shape[0],dtype=float)
StrainEnergyDensity = zeros(shape=xy.shape[0],dtype=float)

#Time integration

#Critical Time estimation Safety Factor 
SafetyFactor = 0.7

##CriticalTimeStep Calculation
CriticalTimeSteps = zeros(shape=xy.shape[0],dtype=float)
for i in range(len(xy)):
    current=xy[i]
    Summation=0.
    for j in NeighborList[i]:
        other = xy[j]
        xi=other-current
        CpEff = 18 * Kappa /norm(xi) / pi / delta**4 #N/mm2/mm/mm4
        Summation+=volume*CpEff#N/mm2/mm/mm
    CriticalTimeSteps[i]=(2*rho / Summation)**0.5
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
        for j in NeighborList[i]:
            other = xy[j]
            xi = other-current
            InfluenceWeight = exp(-norm(xi)**2/delta**2)
            eta = Deformation1Step[j,:] - Deformation1Step[i,:]
            y = xi+eta
            BondExtensionScalar = norm(y)- norm(xi)
            Dilatation[i]+= 3 / WeightedVolume[i] * InfluenceWeight * norm(xi) * BondExtensionScalar * volume
    StrainEnergyDensity += Kappa /2 * Dilatation**2
    
    for i in range(len(xy)):
        current = xy[i]
        for j in NeighborList[i]:
            other = xy[j]
            xi = other-current
            InfluenceWeight = exp(-norm(xi)**2/delta**2)
            eta = Deformation1Step[j,:] - Deformation1Step[i,:]
            y = xi+eta
            
            BondExtensionScalar = norm(y) - norm(xi)
            
            BondExtensionScalarIsotropic = Dilatation[i] * norm(xi) / 3
            BondExtensionScalarDeviatoric = BondExtensionScalar - BondExtensionScalarIsotropic
            
            InternalForceScalarIsotropic = 3 * Kappa * Dilatation[i] * InfluenceWeight * norm(xi) / WeightedVolume[i]
            InternalForceScalarDeviatoric = (15 * Mu / WeightedVolume[i]) * InfluenceWeight * BondExtensionScalarDeviatoric
            InternalForceScalar = InternalForceScalarIsotropic + InternalForceScalarDeviatoric
            
            DeformationDirectionVector = y/norm(y)
            
            ForceInternal[i,:]+= InternalForceScalar*DeformationDirectionVector*volume
            ForceInternal[j,:]-= InternalForceScalar*DeformationDirectionVector*volume
            StrainEnergyDensity[i] += (15 * Mu / WeightedVolume[i] ) / 2 * InfluenceWeight * BondExtensionScalarDeviatoric**2 * volume
    
    #----------
    #Computing Acceleration
    Accel = matmul(InvMassMatrix,ForceInternal*volume+ForceExternal*volume)
    #Second Partial Velocity Update
    Velocity = VelocityHalfStep + (Time1Step-TimeHalfStep) * Accel
    Deformation = Deformation1Step
    
    WorkExternals = (ForceExternal*Deformation).sum(axis=1)
    StrainEnergys = StrainEnergyDensity*volume
    PotentialEnergys = (StrainEnergys - WorkExternals) 
    KineticEnergys = matmul(MassMatrix,Velocity**2).sum(axis=1)*rho*volume
    Lagrangians = KineticEnergys-PotentialEnergys
    Time = Time1Step
    y = xy+Deformation
    # print(y[where((xy[:,0]==50.)*(xy[:,1]==25.)),0])
    # print(Lagrangians.sum())
    # print(Time)
    if not bool(n%10000):
        df = DataFrame()
        df[['X','Y']] = xy
        df[['DispX','DispY']] = Deformation
        df[['VelX','VelY']] = Velocity
        df[['AccelX','AccelY']] = Accel
        df[['InternalForceX','InternalForceY']] = ForceInternal
        df[['ExternalForceX','ExternalForceY']] = ForceExternal
        df['StrainEnergy'] = StrainEnergys
        df['WorkPotential'] = WorkExternals
        df['PotentialEnergy'] = PotentialEnergys
        df['KineticEnergy'] = KineticEnergys
        df['Lagrangian'] = Lagrangians
        df.to_csv('out%d.csv'%n,index_label='ID',float_format='%1.5E')
    with open('pd_out.txt','a') as file:
        file.writelines('Time:%1.7f, Iter:%d, TotalKineticEnergy:%1.5E, TotalPotentialEnergy:%1.5E, Lagrangian:%1.5E,\n' %(Time,n,KineticEnergys.sum(),PotentialEnergys.sum(),Lagrangians.sum()))
    n+=1
