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
    integer ndivx, ndivy, total_point,max_member,i
    real *8 L, W, delta, horizon, thickness, volume
    parameter(ndivx = 100)
    parameter(ndivy = 50)
    parameter(total_point = ndivx * ndivy)
    parameter(L = 1.0)
    parameter(W = 0.5)
    parameter(delta = L/ndivx)
    parameter(horizon = 3*delta)
    parameter(max_member = total_point*36)
    real *8 coord(total_point,2), disp(total_point,2), bforce(total_point,2)
    integer numfam(total_point,1), pointfam(total_point,1), nodefam(max_member,1)
    character(len=60) condition
    character(len=60),dimension(4,1) :: condition_list
    real *8 applied
    real *8 E, nu, kappa, mu , a, b, d, pi
    real *8 idist, nlength, stretch
    thickness = delta
    volume = delta*delta*thickness
    pi = dacos(-1.0d0)
    ! --MATERIAL PROPERTIES--
    E = 200.0
    nu = 1/3
    kappa = E/3/(1-2*nu)
    mu = E/2/(1+nu)
    a = 0.5 * (kappa - 2 * mu)
    b = 6 * mu / pi / thickness / horizon**4
    d = 2 / pi / thickness / horizon**3 

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

    do i = 1, max_member
        nodefam(i,1)=0
    enddo
    call set_coord(coord,L,W,ndivx,ndivy,delta,total_point)
    call set_neigbors(coord,numfam,pointfam,nodefam,horizon,total_point,max_member)
    call preprocess(condition_list, horizon, delta, volume, d, b, idist, nlength, stretch, coord, disp, numfam, pointfam,nodefam,total_point, max_member)
    
end program main

