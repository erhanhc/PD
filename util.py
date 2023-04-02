## PD Classes to handle 
from numpy import where,exp,pi,zeros,eye,matmul
from numpy.linalg import norm,inv

class MaterialPoint:
    def __init__(self,ReferenceCoordinates):
        self.ReferenceCoordinates = ReferenceCoordinates

class MaterialPointInABody:

    def __init__(self,label,Body):
        self.Body = Body
        self.label = label
        self.ReferenceCoordinates = lambda: self.Body.ReferenceCoordinates[self.label]
        self.volume = lambda: self.Body.Volumes[self.label]
        self.GetNeighbors()
        self.Initialize()
        self.WeightedVolume = (self.InfluenceWeights()*norm(self.ReferenceRelativePositionVectors(),axis=1)**2*self.volume()).sum()

    def GetNeighbors(self):
        self.NeighborList = where(
            (norm(self.ReferenceCoordinates()-self.Body.ReferenceCoordinates,axis=1)<=self.Body.delta) *(norm(self.ReferenceCoordinates()-self.Body.ReferenceCoordinates,axis=1)!=0.) 
            )[0]

    def Initialize(self):
        for name in self.Body.config1:
            setattr(
                self,
                name,
                lambda : getattr(self.Body, name)[self.label])

    def ReferenceRelativePositionVectors(self):
        return self.Body.ReferenceCoordinates[self.NeighborList]-self.ReferenceCoordinates()
    def DeformedRelativePositionVectors(self):
        return self.ReferenceRelativePositionVectors() + self.Body.Deformation[self.NeighborList]-self.Deformation()
    def BondExtensionScalars(self):
        return norm(self.DeformedRelativePositionVectors() - self.ReferenceRelativePositionVectors(),axis=1)
    def BondExtensionScalarsVolumetric(self):
        return self.DilatationScalar()*norm(self.ReferenceRelativePositionVectors(),axis=1)/3
    def BondExtensionScalarsDeviatoric(self):
        return self.BondExtensionScalars()-self.BondExtensionScalarsVolumetric()
    def InfluenceWeights(self):
        return exp(-1*(norm(self.ReferenceRelativePositionVectors(),axis=1))**2/self.Body.delta**2)
    def DilatationScalar(self):
        return 3 * ( self.InfluenceWeights() * norm(self.ReferenceRelativePositionVectors(),axis=1) * self.BondExtensionScalars()).sum() * self.volume() / self.WeightedVolume
    def StrainEnergyDensityScalar(self):
        return self.Body.ConstitutiveModel.StrainEnergyDensityScalar(self)
    def DeformedRelativeDeformationDirectionVectors(self):
        return self.DeformedRelativePositionVectors()/norm(self.DeformedRelativePositionVectors(),axis=1).reshape((self.DeformedRelativePositionVectors().shape[0],1))
    def InternalForceScalars(self):
        return self.Body.ConstitutiveModel.InternalForceScalars(self)
    def InternalForceScalarsIsotropic(self):
        return self.Body.ConstitutiveModel.InternalForceScalarsIsotropic(self)
    def InternalForceScalarsDeviatoric(self):
        return self.Body.ConstitutiveModel.InternalForceScalarsDeviatoric(self)
    def InternalForceVectorsHost(self):
        return (self.InternalForceScalars().reshape((self.InternalForceScalars().shape[0],1))*self.DeformedRelativeDeformationDirectionVectors()).sum(axis=0)*self.volume()
    def InternalForceVectorsNeighbors(self):
        return (self.InternalForceScalars().reshape((self.InternalForceScalars().shape[0],1))*self.DeformedRelativeDeformationDirectionVectors())*self.Body.Volumes[self.NeighborList].reshape((self.Body.Volumes[self.NeighborList].shape[0],1))
    def UpdateInternalForceVectors(self):
        self.Body.InternalForce[self.label,:] += self.InternalForceVectorsHost()
        self.Body.InternalForce[self.NeighborList,:] -= self.InternalForceVectorsNeighbors()

class Body:
    config1=[
        'Deformation',
        'Velocity',
        'Acceleration',
        'InternalForce',
        'ExternalForce'
        ]

    def __init__(self,ConstitutiveModel,ReferenceCoordinates,Volumes,delta):

        self.delta=delta
        self.ReferenceCoordinates=ReferenceCoordinates
        self.Volumes = Volumes
        self.ConstitutiveModel=ConstitutiveModel
        [self.__setattr__(name, zeros(shape=self.ReferenceCoordinates.shape))
            for name in self.config1            
            ]
        self.MaterialPointSetup()
        self.MassMatrix = eye(N=len(self.ReferenceCoordinates)) * self.Volumes * self.ConstitutiveModel.Rho#kg
        self.InverseMassMatrix = inv(self.MassMatrix)
        
    def MaterialPointSetup(self):
        self.MaterialPointList = [
            MaterialPointInABody(label, self)
            for label in range(len(self.ReferenceCoordinates))
            ]
        
class LinearPeridynamicElastic:
    def __init__(self,Kappa,Mu,Rho):
        self.Kappa = Kappa
        self.Mu = Mu
        self.Rho = Rho
    
    def StrainEnergyDensityScalar(self,MaterialPointInABody):
        return self.Kappa * MaterialPointInABody.DilatationScalar()**2 / 2 + (15 * self.Mu / MaterialPointInABody.WeightedVolume ) / 2 * (MaterialPointInABody.InfluenceWeights() * MaterialPointInABody.BondExtensionScalarsDeviatoric()**2 * MaterialPointInABody.volume()).sum()
    def InternalForceScalarsIsotropic(self,MaterialPointInABody):
        return 3 * self.Kappa * MaterialPointInABody.DilatationScalar() * MaterialPointInABody.InfluenceWeights() * norm(MaterialPointInABody.ReferenceRelativePositionVectors(),axis=1) / MaterialPointInABody.WeightedVolume
    def InternalForceScalarsDeviatoric(self,MaterialPointInABody):
        return (15 * self.Mu / MaterialPointInABody.WeightedVolume) * MaterialPointInABody.InfluenceWeights() * MaterialPointInABody.BondExtensionScalarsDeviatoric()
    def InternalForceScalars(self,MaterialPointInABody):
        return self.InternalForceScalarsIsotropic(MaterialPointInABody) + self.InternalForceScalarsDeviatoric(MaterialPointInABody)

class ExplicitIntegrator:
    def __init__(self,Body,max_steps=10000):
        self.Body=Body
        self.Time = 0.0
        self.step=0
        self.max_steps=max_steps
        self.SetDeltaTime()
    def SetDeltaTime(self):
        CTS=1.
        for MP in self.Body.MaterialPointList:
            CpEffs = 18 * self.Body.ConstitutiveModel.Kappa / norm(MP.ReferenceRelativePositionVectors(),axis=1) / pi / self.Body.delta**4
            CTSi = (2*self.Body.ConstitutiveModel.Rho / (CpEffs*self.Body.Volumes[MP.NeighborList]).sum()/1000)**0.5
            if CTSi<CTS:
                CTS=CTSi
        self.DeltaTime = CTS*0.7
    
    def integrate(self):
        Time1Step = self.Time+self.DeltaTime
        TimeHalfStep = 0.5*(Time1Step+self.Time)
        VelocityHalfStep = self.Body.Velocity + self.Body.Acceleration*(TimeHalfStep - self.Time)
        self.Body.Deformation += VelocityHalfStep * self.DeltaTime
        self.Body.InternalForce[:,:]=0.
        [MP.UpdateInternalForceVectors() for MP in self.Body.MaterialPointList]
        self.Body.Acceleration = matmul(self.Body.InverseMassMatrix,(self.Body.InternalForce+self.Body.ExternalForce)*self.Body.Volumes)
        self.Body.Velocity = VelocityHalfStep + (Time1Step-TimeHalfStep)*self.Body.Acceleration
        self.Time+=self.DeltaTime
        self.step+=1
if __name__=='__main__':
    from numpy import arange,array
    from time import perf_counter
    Kappa = 1.6 * 10**11 #Pascal = N/m2
    Kappa = Kappa *(1/1000)**2#N/mm2
    Mu = 79.3 * 10**11 #Pascal = N/m2
    Mu = Mu *(1/1000)**2#N/mm2
    Rho = 8050#kg/m3
    Rho=Rho*(1/1000)**3#N/mm2
    metal = LinearPeridynamicElastic(Kappa,Mu,Rho)
    Plate=array(object=[[x,y,z] for x in arange(10) for y in arange(10) for z in arange(10)])
    Volumes = array(object = [[1.] for i in range(len(Plate))])
    example=Body(metal, Plate, Volumes,3.)
    Left_Application = where((Plate[:,0]>=10.-2.))[0]
    Right_Application = where((Plate[:,0]<=2.))[0]
    Total_Applied_Force = 1.#N
    example.ExternalForce[Left_Application,0] = Total_Applied_Force/(len(Left_Application)+1)
    example.ExternalForce[Right_Application,0] = -Total_Applied_Force/(len(Right_Application)+1)
    integrator = ExplicitIntegrator(example)
    while integrator.step<integrator.max_steps:
        integrator.integrate()
        print(integrator.step)