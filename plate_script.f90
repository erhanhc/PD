!Initialize all arrays
!Initialize PD and Geometric Properties
!SCF's
!   Set Coniditions
!   Compute S1, S2, D1, D2
!Apply Initial Conditions
!Apply Boundary Conditions
!Time integration tt=1
!   Preprocess
!   Compute PD Forces
!   Integrate over time 


!TIME INTEGRATION
!Preprocess loop
!for each material point
!   for each neighbor
!       compute nlength
!       compute strecth
!       apply volume correction
!       compute Lambda
!       compute pforcex pforcey
!       compute massvecx massvecy
include 'pd_routines.f90'
program main
    implicit none
    integer ndivx, ndivy, total_point,max_member,i,max_iter,tt
    real *8 L, W, delta, horizon, thickness, volume
    parameter(ndivx = 100)
    parameter(ndivy = 50)
    parameter(total_point = ndivx * ndivy)
    parameter(L = 0.25d0)
    parameter(W = 0.125d0)
    parameter(delta = L/ndivx)
    parameter(horizon = 3.015d0*delta)
    parameter(max_member = total_point*50)
    parameter(max_iter = 1000)
    real *8 coord(total_point,2), disp(total_point,2), bforce(total_point,2),pforce(total_point,2), pforceold(total_point,2), vel(total_point,2), velhalf(total_point,2), velhalfold(total_point,2), accel(total_point,2)
    real *8 DSCF(total_point,2), SSCF(total_point,2), Theta(total_point,1)
    integer numfam(total_point,1), pointfam(total_point,1), nodefam(max_member,1)
    character(len=60) condition
    character(len=60),dimension(5,1) :: condition_list
    real *8 applied
    real *8 E, nu, kappa, mu , a, b, d, pi
    real *8 idist, nlength, stretch, dens
    applied = 1.5d5
    thickness = 0.005d0
    volume = delta*delta*thickness
    pi = dacos(-1.0d0)
    dens = 7850.0d1
    ! --MATERIAL PROPERTIES--
    E = 200.0d9
    nu = 1.0d0/3.0d0
    kappa = E/2.0d0/(1.0d0-nu)
    mu = E/2.0d0/(1.0d0+nu)
    a = 0.5d0 * (kappa - 2.0d0 * mu)
    b = 6.0d0 * mu / pi / thickness / horizon**4
    d = 2.0d0 / pi / thickness / horizon**3 
    print*, 'Mass Density ',dens, '[kg/m3]'
    print*,'Youngs Modulus',E,' [Pa]'
    print*,'Poisson Ratio ',nu
    print*,'Bulk Modulus ',kappa, '[Pa]'
    print*,'Shear Modulus ',mu, '[Pa]'
    print*, 'Parameter a ', a, '[Pa]'
    print*,'Parameter b ',b, '[N/m7]'
    print*,'Parameter d ',d, '[1/m4]'
    do i = 1, total_point
        pointfam(i,1)=0
        numfam(i,1)=0
        coord(i,1)=0.0
        coord(i,2)=0.0
    enddo

    condition_list(1,1) = 'uniaxial stretch x'
    condition_list(2,1) = 'uniaxial stretch y'
    condition_list(3,1) = 'simple shear in x-y'
    condition_list(4,1) = "uniaxial tensile loading"
    condition_list(5,1) = 'displacement uniaxial tensile'

    do i = 1, max_member
        nodefam(i,1)=0
    enddo
    call set_coord(coord,L,W,ndivx,ndivy,delta,total_point)
    call set_neigbors(coord,numfam,pointfam,nodefam,horizon,total_point,max_member)
    call preprocess(condition_list, horizon, delta, volume, d, b, a, coord, disp, vel, pforce, bforce, numfam, pointfam,nodefam,total_point, max_member,DSCF,SSCF,mu)
    
    do i = 1, total_point
        disp(i,1) = 0.0
        disp(i,2) = 0.0
        bforce(i,1) = 0.0
        bforce(i,2) = 0.0
        pforce(i,1) = 0.0
        pforce(i,2) = 0.0
        vel(i,1) = 0.0
        vel(i,2) = 0.0
        velhalf(i,1) = 0.0
        velhalf(i,2) = 0.0
        velhalfold(i,1) = 0.0
        velhalfold(i,2) = 0.0
        pforceold(i,1) = 0.0
        pforceold(i,2) = 0.0
        accel(i,1) = 0.0
        accel(i,1) = 0.0
    enddo
    condition = 'displacement uniaxial tensile'
    applied = 1.5d5
    call set_conditions(condition, coord, disp, bforce, applied, delta, total_point,vel,pforce, horizon)
    ! call preprocess_with_SCF(horizon, delta, volume, d, b, a, coord, disp, numfam, pointfam,nodefam,total_point, max_member,DSCF,SSCF,Theta)
    ! call iterate(max_iter,horizon, delta, volume, d, b, a, coord, disp, numfam, pointfam,nodefam,total_point, max_member, DSCF, SSCF, Theta, pforce, bforce, pforceold, vel, velhalf, velhalfold)
    call time_integration(horizon, delta, volume, d, b, a, dens, coord, disp, vel, accel, pforce, bforce, numfam, pointfam,nodefam,total_point, max_member,DSCF,SSCF,Theta,100000000,1.0d-6)
end program main

