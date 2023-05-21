## PD Classes to handle
with_cupy=False
if with_cupy: 
    from cupy import where,exp,pi,zeros,eye,matmul
    from cupy.linalg import norm,inv
else:
    from numpy import where,exp,pi,zeros,eye,matmul,dot
    from numpy.linalg import norm,inv
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

class MaterialPoint:
    def __init__(self,ReferenceCoordinates):
        self.ReferenceCoordinates = ReferenceCoordinates

class MaterialPointInABody:

    def __init__(self,label,Body):
        self.Body = Body
        self.label = label
        self.ReferenceCoordinates = self.Body.ReferenceCoordinates[self.label]
        self.volume = self.Body.Volumes[self.label]
        self.GetNeighbors()
        
        
    
    def ReferenceRelativePositionVectors(self):
        return self.Body.ReferenceCoordinates[self.NeighborList]-self.ReferenceCoordinates
    def DeformedRelativePositionVectors(self):
        return self.ReferenceRelativePositionVectors() + self.Body.Deformation[self.NeighborList]-self.Deformation()
    def StretchScalars(self):
        return ( norm (  self.DeformedRelativePositionVectors(), axis=1 ) - norm( self.ReferenceRelativePositionVectors(),axis=1))/norm(self.ReferenceRelativePositionVectors(),axis=1)
    def LambdaScalars(self):
        return (self.ReferenceRelativePositionVectors() * self.DeformedRelativePositionVectors()).sum(axis=1) / norm(self.ReferenceRelativePositionVectors(),axis=1) / norm(self.DeformedRelativePositionVectors(),axis=1)
    def DilatationScalar(self):
        return self.Body.ConstitutiveModel.d * self.Body.delta * (self.StretchScalars() * self.LambdaScalars() * self.Body.Volumes[self.NeighborList].flatten()).sum()
    def StrainEnergyDensityScalar(self):
        return self.Body.ConstitutiveModel.a * self.DilatationScalar()**2 + self.Body.ConstitutiveModel.b * (self.Body.delta/norm(self.ReferenceRelativePositionVectors(),axis=1) * (norm(self.DeformedRelativePositionVectors(),axis=1)-norm(self.ReferenceRelativePositionVectors(),axis=1))**2*self.Body.Volumes[self.NeighborList]).sum()
    def DeformedRelativeDeformationDirectionVectors(self):
        return self.DeformedRelativePositionVectors()/norm(self.DeformedRelativePositionVectors(),axis=1).reshape((self.DeformedRelativePositionVectors().shape[0],1))
    def InternalForceScalars(self):
        return ( 2 * self.Body.ConstitutiveModel.a * self.Body.ConstitutiveModel.d * self.Body.delta * self.LambdaScalars() * self.DilatationScalar() / norm(self.ReferenceRelativePositionVectors(),axis=1) + 2 * self.Body.ConstitutiveModel.b * self.Body.delta * self.StretchScalars())
    def InternalForceVectorsHost(self):
        return (self.DeformedRelativeDeformationDirectionVectors() * self.InternalForceScalars().reshape((self.InternalForceScalars().shape[0],1)) * self.Body.Volumes[self.NeighborList].reshape((self.Body.Volumes[self.NeighborList].shape[0],1))).sum(axis=0)
    def InternalForceVectorsNeighbors(self):
        return self.DeformedRelativeDeformationDirectionVectors() * self.InternalForceScalars().reshape((self.InternalForceScalars().shape[0],1)) * self.volume
    def UpdateInternalForceVectors(self):
        self.Body.InternalForce[self.label,:] += self.InternalForceVectorsHost()
        self.Body.InternalForce[self.NeighborList,:] -= self.InternalForceVectorsNeighbors()
    
    def GetNeighbors(self):
        self.NeighborList = where(
            (norm(self.ReferenceCoordinates-self.Body.ReferenceCoordinates,axis=1)<=self.Body.delta) *(norm(self.ReferenceCoordinates-self.Body.ReferenceCoordinates,axis=1)!=0.) 
            )[0]

    def Deformation(self):
        return self.Body.Deformation[self.label]
            
class IsotropicMaterial3D:
    def __init__(self,Kappa,Mu,Rho,delta):
        self.Kappa = Kappa
        self.Mu = Mu
        self.Rho = Rho
        self.a = 1 / 2 * ( Kappa - 5 * Mu / 3)
        self.b = 15 * Mu / 2 / pi / delta**5
        self.d = 9 / 4 / pi / delta**4

class IsotropicMaterialPlaneStress:
    def __init__(self,Kappa,Mu,Rho,delta,thickness):
        self.Kappa = Kappa
        self.Mu = Mu
        self.Rho = Rho
        self.a = 1 / 2 * (Kappa-2*Mu)
        self.b = 6 * Mu / pi / thickness / delta**4
        self.d = 2 / pi/ thickness / delta**3

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
        
def main_func():
    if with_cupy:
        from cupy import arange,array
    else:
        from numpy import arange,array
    from time import perf_counter
    from matplotlib import pyplot as plt
    Kappa = 140 * 10**9 #Pascal = N/m2
    Kappa = Kappa *(1/1000)**2#N/mm2
    Mu = 80 * 10**9 #Pascal = N/m2
    Mu = Mu *(1/1000)**2#N/mm2
    Rho = 8050#kg/m3
    Rho=Rho*(1/1000)**3#kg/mm3
    L=100.
    W=L/2.
    t = 3.
    dx = 0.5*4
    delta = dx*3.1
    metal = IsotropicMaterialPlaneStress(Kappa,Mu,Rho,delta,t)
    Plate=array([[x,y] for x in arange(dx/2,L+dx/2,dx) for y in arange(dx/2,W+dx/2,dx)])
    Volumes = array([[L*W*t/len(Plate)] for i in range(len(Plate))])
    example=Body(metal, Plate, Volumes,delta)
    Left_Application = where((Plate[:,0]>=L-3*dx))[0]
    Right_Application = where((Plate[:,0]<=3*dx))[0]
    Total_Applied_Force = 1*W*t
    example.ExternalForce[Left_Application,0] = Total_Applied_Force/(len(Left_Application))
    example.ExternalForce[Right_Application,0] = -Total_Applied_Force/(len(Right_Application))
    integrator = ExplicitIntegrator(example)
    return example,integrator,Left_Application,Right_Application,Total_Applied_Force
if __name__=='__main__':
    example,integrator,Left_Application,Right_Application,Total_Applied_Force=main_func()
    from matplotlib import pyplot as plt
    integrator.DeltaTime=1e-5
    while integrator.step<10000:
        integrator.integrate()
        if not integrator.step%250:
            print(integrator.step)
            fig,ax = plt.subplots(2,2)
            fig.set_size_inches(18,9)
            ax1 = ax[0][0]
            im1=ax1.scatter(example.ReferenceCoordinates[:,0] + example.Deformation[:,0],example.ReferenceCoordinates[:,1] + example.Deformation[:,1],c=norm(example.Deformation,axis=1))
            fig.colorbar(im1,ax=ax1)
            ax1.set_title('Deformed Configuration')
            ax1.set_xlabel('Deformed Coordinates - x')
            ax1.set_ylabel('Deformed Coordinates - y')
            ax2=ax[0][1]
            im2=ax2.scatter(example.ReferenceCoordinates[:,0],example.ReferenceCoordinates[:,1],c=example.ExternalForce[:,0])
            fig.colorbar(im2,ax=ax2)
            ax2.set_title('Applied External Loads x Magnitude')
            ax3=ax[1][0]
            im3=ax3.scatter(example.ReferenceCoordinates[:,0],example.ReferenceCoordinates[:,1],c=example.InternalForce[:,0])
            fig.colorbar(im3,ax=ax3)
            ax3.set_title('Internal Force x Magnitude')
            ax4=ax[1][1]
            im4=ax4.scatter(example.ReferenceCoordinates[:,0],example.ReferenceCoordinates[:,1],c=norm(example.Velocity,axis=1))
            fig.colorbar(im4,ax=ax4)
            ax4.set_title('Speed')
            fig
            fig.savefig(f'Madenci_Oterkus_Plate_Example_v2_{integrator.step}')
            example.ExternalForce[Left_Application,0] += Total_Applied_Force/(len(Left_Application))
            example.ExternalForce[Right_Application,0] -= Total_Applied_Force/(len(Right_Application))
