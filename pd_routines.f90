subroutine set_coord(coord,L,W,ndivx,ndivy,delta,total_point)
    implicit none
    !This subroutine gets the plate size and discritizes based on ndivx and ndivy and 
    !sets collocation point coord's and their global material point id's. 
    !current_point  :   placeholder for id setting for the current x,y location
    !i,j            :   division increments in the do loops for ndivx, ndivy
    !coord          :   array that relates the coordinate locations to material point id's
    !L              :   Length of the plate
    !W              :   Width of the plate
    !ndivx, ndivy   :   Divisions in x and y
    !delta          :   Spacing delta = L/ndivx
    !total_point    :   Total number of points = ndivx * ndivy
    logical echo
    parameter(echo = .FALSE.)
    integer current_point, i, j, ndivx, ndivy, total_point
    real *8 coord(total_point, 2), L, W, delta
    current_point = 1
    do i = 1, ndivx
        do j = 1, ndivy
            coord(current_point,1) = -0.5 * L + ( delta / 2.0) + (i-1) * delta
            coord(current_point,2) = -0.5 * W + ( delta / 2.0) + (i-j) * delta
            if (echo) then
                print*, 'ID ',current_point,' x= ',coord(current_point,1), ' y= ',coord(current_point,2)
            endif
            current_point = current_point + 1
        enddo
    enddo
    return
end


subroutine set_neigbors(coord,numfam,pointfam,nodefam,horizon,total_point,max_member)
    implicit none
    !This subroutine gets all material point coordinate locations and looks for proximity and
    !sets neighborhood relations
    !current_point  :   host point id
    !other_point    :   neighbor candidate
    !nodefam        :   neighborhood relation (bond) list where neighbor id's are stored
    !pointfam       :   array that relates host material point id to starting index in nodefam (neighborhood list)
    !numfam         :   array that relates host material point id to its neighbor count
    logical echo
    parameter(echo = .FALSE.)
    integer current_point, other_point 
    integer, intent(in) :: max_member, total_point
    integer,intent(inout) :: numfam(total_point,1), pointfam(total_point,1), nodefam(max_member,1)
    real *8 , intent(in) ::  horizon
    real *8 distance
    real*8, dimension(total_point,2), intent(in) :: coord
    do current_point = 1, total_point
        if (current_point.eq.1) then
            pointfam(current_point,1) = 1
        else
            pointfam(current_point,1) = pointfam(current_point-1,1) + numfam(current_point-1,1)
        endif
        do other_point = 1, total_point
            if (current_point /= other_point) then
                distance = dsqrt((coord(current_point,1)-coord(other_point,1))**2+(coord(current_point,2)-coord(other_point,2))**2)
                if (distance.le.horizon) then
                    numfam(current_point,1) = numfam(current_point,1) + 1
                    nodefam(pointfam(current_point,1)+numfam(current_point,1)-1,1) = other_point
                    if (echo) then
                        print*, current_point, other_point, distance, numfam(current_point,1)
                    endif
                endif
            endif
        enddo
    enddo
    return
end

subroutine set_conditions(condition, coord, disp, bforce, applied, delta, total_point)
    implicit none
    ! This subroutine sets specific boundary conditions for surface correction factor calculations.
    ! Uniaxial tensile loading is added for the benchmarking.
    logical echo
    parameter(echo = .FALSE.)
    integer, intent(in) :: total_point
    integer current_point
    character(len=60), intent(in) :: condition
    real *8 ,intent(in) :: coord(total_point,2), applied, delta
    real *8 , intent(inout) :: bforce(total_point,2), disp(total_point,2)
    real *8 left_bound, right_bound, disp_grad
    do current_point=1, total_point

        ! -----------UNIAXIAL STRETCH IN X DIRECTION-----------
        if (condition == 'uniaxial stretch x') then
            disp_grad = applied
            disp(current_point,1) = disp_grad * coord(current_point,1)
            disp(current_point,2) = 0.0d0 
            if (echo) then
                print*, condition, current_point, 'at', coord(current_point,1), 'disp is', disp(current_point,1)
            endif
        !-------------------------------------------------------

        ! -----------UNIAXIAL STRETCH IN Y DIRECTION-----------
        else if (condition =='uniaxial stretch y') then
            disp_grad = applied
            disp(current_point,1) = 0.0
            disp(current_point,2) = disp_grad * coord(current_point,2)
            if (echo) then
                print*, current_point, 'at y=', coord(current_point,2), 'dispy is', disp(current_point,2)
            endif
        !-------------------------------------------------------

        ! -----------SIMPLE SHEAR STRECTH CASE-----------
        else if (condition == 'simple shear in x-y') then
            disp_grad = applied
            disp(current_point,1) = 0.5 * disp_grad * coord(current_point,1)
            disp(current_point,2) = -0.5 * disp_grad * coord(current_point,2)
            if (echo) then
                print*, current_point, 'at x ', coord(current_point,1), 'disp x is', disp(current_point,1)
                print*, current_point, 'at y ', coord(current_point,2), 'disp y is', disp(current_point,2)
            endif
        !-------------------------------------------------------

        ! -----------UNIAXIAL TENSILE LOADING IN X DIRECTION-----------
        else if (condition == 'uniaxial tensile loading') then
            disp(current_point,1) = 0.0d0
            disp(current_point,2) = 0.0d0
            left_bound = minval(coord(:,1))
            right_bound = maxval(coord(:,1))
            if (coord(current_point,1) == right_bound) then
                bforce(current_point,1) = -1.0d0 * applied / delta
                bforce(current_point,2) = 0.0d0
            else if (coord(current_point,1) == left_bound) then
                bforce(current_point,1) = 1.0d0 * applied / delta
                bforce(current_point,2) = 0.0d0
            else 
                bforce(current_point,1) = 0.0d0
                bforce(current_point,2) = 0.0d0
            endif
            if (echo) then
                print*, current_point, 'applied x force', bforce(current_point,1)
                print*, current_point, 'applied y force', bforce(current_point,2)
            endif
        endif
        !-------------------------------------------------------
            
    enddo
    return
end

subroutine compute_fac_volume_corr(fac_vol, idist, horizon, delta)
    implicit none
    logical echo
    parameter(echo = .FALSE.)
    real *8 ,intent(in) :: idist, horizon, delta
    real *8 fac_vol
    if (idist < horizon - delta / 2) then
        fac_vol = 1.
    else if (idist <= horizon) then
        fac_vol = (horizon + delta / 2 - idist) / delta
    else
        fac_vol = 0.
    endif
    if (echo) then
            print*, 'volume correction fac:', fac_vol,'idist:',idist, 'range:', horizon - delta / 2
    endif
    return
end subroutine compute_fac_volume_corr

subroutine dot_product(vec_1,vec_2, product)
    real*8, intent(in) :: vec_1(1,2), vec_2(1,2)
    real*8 product 
    product = 0.
    do i = 1,2
        product = product + vec_1(1,i)*vec_2(1,i)
    enddo
return
end subroutine dot_product

subroutine compute_kinematics(idist, nlength,stretch,Lambda, x_k, x_j, u_k, u_j) 
    implicit none
    logical echo
    parameter(echo = .FALSE.)
    real*8, intent(in) ::  x_k(1,2), x_j(1,2), u_k(1,2), u_j(1,2)
    real*8  nlength, idist, stretch, Lambda
    real*8 y_k(1,2), y_j(1,2)
    y_k = x_k + u_k
    y_j = x_j + u_j
    call dot_product(x_j-x_k,x_j-x_k,idist)
    idist = dsqrt(idist)
    call dot_product(y_j-y_k,y_j-y_k,nlength)
    nlength = dsqrt(nlength)
    call dot_product((y_j-y_k)/nlength,(x_j-x_k)/idist,Lambda)
    stretch = (nlength - idist) / idist
    if (echo) then
        print*, 'idist:',idist, 'nlength:',nlength, 'stretch:', stretch, 'Lambda:', Lambda
    endif
    return
end subroutine compute_kinematics

subroutine compute_Theta_k_j(Theta_k_j, stretch_k_j, Lambda_k_j, volume_j)
    implicit none
    real*8, intent(in) ::  stretch_k_j, Lambda_k_j, volume_j
    real*8 Theta_k_j
    Theta_k_j = stretch_k_j * Lambda_k_j * volume_j
    return
end subroutine compute_Theta_k_j

subroutine compute_Distort_k_j(Distort_k_j,idist,nlength,volume_j)
    implicit none
    real*8, intent(in) :: idist,nlength,volume_j
    real*8 Distort_k_j
    Distort_k_j = (nlength-idist)**2 * volume_j /idist
    return
end subroutine compute_Distort_k_j

subroutine set_surface_correction_factors(current_point,condition,disp_grad,SED,Theta_k,mu,total_point,DSCF,SSCF)
    implicit none
    logical echo
    parameter(echo = .FALSE.)
    real*8, intent(in) :: disp_grad , SED , Theta_k , mu
    real*8 DSCF(total_point,2), SSCF(total_point,2)
    integer, intent(in) :: current_point,total_point
    character(len=60), intent(in) :: condition
    real*8 d_1 , d_2 , s_1 , s_2
    
    ! -----------UNIAXIAL STRETCH IN X DIRECTION-----------
    if (condition == 'uniaxial stretch x') then
        d_1 = disp_grad / Theta_k
        if (echo) then
            print*, d_1, 'Theta Correction D1 - xdir'
        endif
        DSCF(current_point,1)=d_1
    !-------------------------------------------------------
    ! -----------UNIAXIAL STRETCH IN Y DIRECTION-----------
    else if (condition =='uniaxial stretch y') then
        d_2 = disp_grad / Theta_k
        if (echo) then
            print*, d_2, 'Theta Correction D2 - ydir'
        endif
        DSCF(current_point,2)=d_2
    !-------------------------------------------------------

    ! -----------SIMPLE SHEAR STRECTH CASE-----------
    else if (condition == 'simple shear in x-y') then
        s_1 = 0.5*disp_grad**2*mu / SED
        s_2 = 0.5*disp_grad**2*mu / SED
        if (echo) then
            print*, s_1, 'Strain Energy Correction  S1 - xdir'
            print*, s_2, 'Strain Energy Correction  S2 - ydir'
        endif
        SSCF(current_point,1) = s_1
        SSCF(current_point,2) = s_2
    endif
    !-------------------------------------------------------
    return
end subroutine set_surface_correction_factors


subroutine preprocess(condition_list, horizon, delta, volume, d, b, a, coord, disp, numfam, pointfam,nodefam,total_point, max_member,DSCF,SSCF,mu)
    implicit none
    logical echo
    parameter(echo = .FALSE.)
    real*8, intent(in) :: horizon, delta, volume, d, b, a,mu
    integer i
    real*8 idist, nlength, stretch, fac_vol, applied
    real*8 coord(total_point,2), disp(total_point,2), bforce(total_point,2), DSCF(total_point,2), SSCF(total_point,2)
    integer numfam(total_point,1), pointfam(total_point,1), nodefam(max_member,1)
    integer current_point, other_point, total_point, max_member
    character(len=60) condition
    character(len=60),dimension(4,1) , intent(in) :: condition_list

    real*8 Lambda, Theta_k, Theta_k_j, Distort_k_j, Distort_k, SED_k
    do i = 1, 3
        condition = condition_list(i,1)
        applied = 0.001
        call set_conditions(condition, coord, disp, bforce, applied, delta, total_point)
        do current_point = 1, total_point
            Distort_k = 0.
            Theta_k = 0.
            SED_k = 0.
            do other_point = 1, numfam(current_point,1)
                idist = 0.
                nlength = 0.
                stretch = 0.
                Lambda = 0.
                Theta_k_j = 0.
                Distort_k_j = 0.
                call compute_kinematics(idist, nlength, stretch, Lambda, coord(current_point,:), coord(nodefam(pointfam(current_point,1)+other_point-1,1),:) , disp(current_point,:), disp(nodefam(pointfam(current_point,1)+other_point-1,1),:))
                call compute_fac_volume_corr(fac_vol, idist, horizon, delta)
                call compute_Theta_k_j(Theta_k_j, stretch, Lambda, volume*fac_vol)
                Theta_k = Theta_k + d * horizon * Theta_k_j
                call compute_Distort_k_j(Distort_k_j,idist,nlength,volume*fac_vol)
                Distort_k = Distort_k + b * horizon * Distort_k_j
                if (echo) then
                    print*, 'x_k:', coord(current_point,:), 'y_k:', coord(current_point,:)+disp(current_point,:)
                    print*, 'x_j:', coord(nodefam(pointfam(current_point,1)+other_point-1,1),:), 'y_j:', coord(nodefam(pointfam(current_point,1)+other_point-1,1),:)+disp(nodefam(pointfam(current_point,1)+other_point-1,1),:)                
                    print*, 'Theta_k_j:', Theta_k_j, 'Theta_k:', Theta_k
                    print*, 'Distort_k_j:', Distort_k_j, 'Distort_k:',Distort_k
                endif
            enddo
            SED_k = a * Theta_k **2 + Distort_k
            call set_surface_correction_factors(current_point,condition,applied,SED_k,Theta_k,mu,total_point,DSCF,SSCF)
            if (echo) then
                print*, 'SED:', SED_k
            endif
        enddo
    enddo
    return
end subroutine preprocess

subroutine surface_correction_vector(current_DSCF,current_SSCF,current_coord,other_DSCF,other_SSCF,other_coord,idist,Gb,Gd)
    implicit none
    !from page. 74
    logical echo
    parameter(echo = .FALSE.)
    real*8, intent(in) :: idist
    real*8, intent(in) :: current_DSCF(1,2),current_SSCF(1,2),current_coord(1,2),other_DSCF(1,2),other_SSCF(1,2),other_coord(1,2)
    real*8 gdx, gdy, gbx, gby, nx, ny
    real*8 Gd, Gb
    gdx = (current_DSCF(1,1) + other_DSCF(1,1)) * 0.5
    gdy = (current_DSCF(1,2) + other_DSCF(1,2)) * 0.5
    gbx = (current_SSCF(1,1) + other_SSCF(1,1)) * 0.5
    gby = (current_SSCF(1,2) + other_SSCF(1,2)) * 0.5
    nx = (other_coord(1,1) - current_coord(1,1)) / idist
    ny = (other_coord(1,2) - current_coord(1,2)) / idist
    Gd = ( ( nx / gdx )**2 + ( ny / gdy )**2 )** (-0.5)
    Gb = ( ( nx / gbx )**2 + ( ny / gby )**2 )** (-0.5)
    if (echo) then
        print*, 'Current Node:',current_coord
        print*, 'Neighbor Node:',other_coord
        print*, 'Initial Distance:',idist
        print*, 'Current Node DSCF:',current_DSCF
        print*, 'Current Node SSCF:',current_SSCF
        print*, 'Neighbor Node DSCF:',other_DSCF
        print*, 'Neighbor Node SSCF:',other_SSCF
        print*, 'Normal Vector:',nx,ny
        print*, 'Gd:', Gd
        print*, 'Gb:', Gb
    endif
end subroutine surface_correction_vector

subroutine compute_PD_forces(horizon, volume, d, b, a, fac_vol, idist, nlength, stretch, Gd, Gb, Lambda, Theta_current, Theta_other, coord_current, coord_other, disp_current, disp_other, pforce)
    implicit none
    logical echo
    parameter(echo = .FALSE.)
    real*8, intent(in) :: horizon, volume, d, b, a, fac_vol, idist, nlength, Gd, Gb, Lambda, stretch
    real*8, intent(in) :: Theta_current(1,1), Theta_other(1,1), coord_current(1,2), coord_other(1,2), disp_current(1,2), disp_other(1,2)
    real*8 pforce(1,2)
    real*8 AA, BB, delyx, delyy, tkjx, tkjy, tjkx, tjky
    AA = 2 * horizon * (d * Gd * Lambda / idist * (a * Theta_current(1,1)) + b * Gb * stretch )
    BB = 2 * horizon * (d * Gd * Lambda / idist * (a * Theta_other(1,1)) + b * Gb * stretch )
    delyx = ( (disp_other(1,1)+coord_other(1,1)) - ( disp_current(1,1)+coord_current(1,1) ) ) 
    delyy = ( (disp_other(1,2)+coord_other(1,2)) - ( disp_current(1,2)+coord_current(1,2) ) ) 
    tkjx = delyx * 0.5 * A  / nlength
    tkjy = delyy * 0.5 * A  / nlength
    tjkx = -delyx * 0.5 * B  / nlength
    tjky = -delyy * 0.5 * B  / nlength
    pforce(1,1) = pforce(1,1) + (tkjx - tjkx) * volume * fac_vol
    pforce(1,2) = pforce(1,2) + (tkjy - tjky) * volume * fac_vol
end subroutine compute_PD_forces

subroutine iterate(max_iter,horizon, delta, volume, d, b, a, coord, disp, numfam, pointfam,nodefam,total_point, max_member, DSCF, SSCF, Theta, pforce, bforce, pforceold, vel, velhalf, velhalfold)
    implicit none
    logical echo
    integer, intent(in) :: max_iter,total_point, max_member
    parameter(echo = .FALSE.)
    real*8, intent(in) :: horizon, delta, volume, d, b, a
    real*8 idist, nlength, stretch, fac_vol
    real*8 coord(total_point,2), disp(total_point,2), bforce(total_point,2), DSCF(total_point,2), SSCF(total_point,2), Theta(total_point,1), pforce(total_point,2)
    integer ,intent(in) :: numfam(total_point,1), pointfam(total_point,1), nodefam(max_member,1)
    integer current_point, other_point, tt
    real*8 Lambda
    real*8 Gd, Gb, Kijx, Kijy
    real*8 massvec(total_point,2)
    real*8 vel(total_point,2), velhalf(total_point,2), velhalfold(total_point,2), pforceold(total_point,2)
    !page 142
    do tt=1,max_iter
        call preprocess_with_SCF(horizon, delta, volume, d, b, a, coord, disp, numfam, pointfam,nodefam,total_point, max_member,DSCF,SSCF,Theta)
        do current_point = 1, total_point
            massvec(current_point,1) = 0.0
            massvec(current_point,2) = 0.0
            pforce(current_point,1) = 0.0
            pforce(current_point,2) = 0.0
            do other_point = 1, numfam(current_point,1)
                idist = 0.
                nlength = 0.
                stretch = 0.
                Lambda = 0.
                Gd = 0.
                Gb = 0.
                Kijx = 0.
                Kijy = 0.
                call compute_kinematics(idist, nlength, stretch, Lambda, coord(current_point,:), coord(nodefam(pointfam(current_point,1)+other_point-1,1),:) , disp(current_point,:), disp(nodefam(pointfam(current_point,1)+other_point-1,1),:))
                call compute_fac_volume_corr(fac_vol, idist, horizon, delta)
                call surface_correction_vector(DSCF(current_point,:),SSCF(current_point,:),coord(current_point,:),DSCF(nodefam(pointfam(current_point,1)+other_point-1,1),:),SSCF(nodefam(pointfam(current_point,1)+other_point-1,1),:),coord(nodefam(pointfam(current_point,1)+other_point-1,1),:),idist,Gb,Gd)
                call compute_PD_forces(horizon, volume, d, b, a, fac_vol, idist, nlength, stretch, Gd, Gb, Lambda, Theta(current_point,:), Theta(nodefam(pointfam(current_point,1)+other_point-1,1),:), coord(current_point,:), coord(nodefam(pointfam(current_point,1)+other_point-1,1),:), disp(current_point,:), disp(nodefam(pointfam(current_point,1)+other_point-1,1),:), pforce(current_point,:))
                call compute_ADR_Kij(coord(current_point,:),coord(nodefam(pointfam(current_point,1)+other_point-1,1),:), horizon, idist, a, b, d, Gb, Gd, fac_vol * volume, fac_vol * volume, Kijx, Kijy)
                massvec(current_point,1) = massvec(current_point,1) +  0.25 * 5 * Kijx
                massvec(current_point,2) = massvec(current_point,2) + 0.25 * 5 * Kijy
                if (echo) then
                    print*, 'Node: ',current_point,' massvecx = ',massvec(current_point,1), ' massvecy = ',massvec(current_point,2)
                endif 
            enddo
        if (echo) then
            print*, 'PD Force on ',current_point,' in x-dir:',pforce(current_point,1)
            print*, 'PD Force on ',current_point,' in y-dir:',pforce(current_point,2)
        endif
        enddo
        call ADR(tt, total_point, disp, vel, velhalf, velhalfold, pforce, bforce, pforceold, massvec)
        print*, 'ux: ',disp(2500,1), ' uy: ', disp(2500,2)
    enddo
end subroutine iterate

subroutine ADR(tt, total_point, disp, vel, velhalf, velhalfold, pforce, bforce, pforceold, massvec )
    implicit none
    logical echo
    parameter(echo = .FALSE.)
    integer, intent(in) :: tt, total_point
    integer current_point
    real*8 disp(total_point,2), vel(total_point,2), velhalf(total_point,2), velhalfold(total_point,2)
    real*8 pforce(total_point,2), bforce(total_point,2), pforceold(total_point,2), massvec(total_point,2)
    real*8 cn_num, cn_denom, cn
    if (tt.eq.1) then 
        do current_point = 1, total_point
            !starting ADR process with tt=1, page 137
            vel(current_point,1) = 0.
            vel(current_point,2) = 0.
            velhalf(current_point,1) = ( pforce(current_point,1) + bforce(current_point,1) ) / massvec(current_point,1)
            velhalf(current_point,2) = ( pforce(current_point,2) + bforce(current_point,2) ) / massvec(current_point,2)
            velhalfold(current_point,1) = velhalf(current_point,1)
            velhalfold(current_point,2) = velhalf(current_point,2)
            pforceold(current_point,1) = pforce(current_point,1)
            pforceold(current_point,2) = pforce(current_point,2)
            disp(current_point,1) =  disp(current_point,1) + velhalf(current_point,1) 
            disp(current_point,2) =  disp(current_point,2) + velhalf(current_point,2) 
            if (echo) then
                print*, total_point
                print*, 'displacement for ',current_point,' :', disp(current_point,:)
                print*, 'vel_n+1/2 for ',current_point,' :', velhalf(current_point,:)
                print*, 'pforce for ',current_point,' :', pforce(current_point,:)
                print*, 'bforce for ',current_point,' :', bforce(current_point,:)
            endif
        enddo
    elseif (tt.gt.1) then
        cn_num = 0.0
        cn_denom = 0.0
        cn = 0.0
        do current_point = 1,total_point
            call compute_ADR_cn(velhalfold(current_point,:), disp(current_point,:), pforce(current_point,:), pforceold(current_point,:), massvec(current_point,:), cn_num, cn_denom )
        enddo
        if (cn_denom.ne.0.0) then
            if ((cn_num/cn_denom).gt.0.0) then
                cn = 2.0 * sqrt(cn_num / cn_denom)
            else
                print*, 'Damping Coefficient, cn < 0. at iteration: ',tt,' ',cn_num,' / ',cn_denom
                cn = 0.0
            endif
        else
            print*, 'Damping Coefficient Denominator, cn_denom = 0. at iteration: ',tt
            cn = 0.0
        endif
        if (cn.gt.2.0) then
            print*, 'Damping Coefficient, cn_denom > 2. at iteration: ',tt
            cn = 1.9
        endif
        if (echo) then
            print*,'Damping Coefficient, cn = ',cn
        endif
        do current_point = 1,total_point
            !page 137
            velhalf(current_point,1) = ( ( 2 - cn ) * velhalfold(current_point,1) + 2 * ( pforce(current_point,1) + bforce(current_point,1) / massvec(current_point,1) )) / ( 2 + cn)
            velhalf(current_point,1) = ( ( 2 - cn ) * velhalfold(current_point,2) + 2 * ( pforce(current_point,2) + bforce(current_point,2) / massvec(current_point,2) )) / ( 2 + cn)
            vel(current_point,1) = 0.5 * ( velhalf(current_point,1) + velhalf(current_point,1) )
            vel(current_point,2) = 0.5 * ( velhalf(current_point,2) + velhalf(current_point,2) )
            disp(current_point,1) = disp(current_point,1) + velhalf(current_point,1)
            disp(current_point,2) = disp(current_point,2) + velhalf(current_point,2)
            velhalfold(current_point,1) = velhalf(current_point,1)
            velhalfold(current_point,2) = velhalf(current_point,2)
            pforceold(current_point,1) = pforce(current_point,1)
            pforceold(current_point,2) = pforce(current_point,2)
        enddo
    endif
end subroutine ADR

subroutine compute_ADR_cn(velhalfold, disp, pforce, pforceold, massvec, cn_num, cn_denom )
    implicit none
    real*8, intent(in) :: velhalfold(1,2), disp(1,2), pforce(1,2), pforceold(1,2), massvec(1,2)
    real*8 cn_num, cn_denom
    logical echo
    parameter(echo = .FALSE.)
    if (echo) then
        print*, '--ADR Damping Coeff--'
        print*, 'Displacement Vector: ',disp(1,:)
        print*, 'PD Bond Force Vector: ',pforce(1,:)
        print*, 'Stable Mass Vector: ',massvec(1,:)
        print*, 'PD Bond Force Vector n-1: ',pforceold(1,:)
        print*, 'Velocity Vector n-1/2 : ',velhalfold(1,:)
        print*, '---------------------'
    endif
    if (velhalfold(1,1).ne.0.0) then
        cn_num = cn_num + ( disp(1,1) * ( -1.0 * ( pforce(1,1) / massvec(1,1) - pforceold(1,1) / massvec(1,1) ) / velhalfold(1,1) ) * disp(1,1) )
    else
        print*, 'velhalfold(1,1) = 0'
    endif
    if (velhalfold(1,2).ne.0.0) then
        cn_num = cn_num + ( disp(1,2) * ( -1.0 * ( pforce(1,2) / massvec(1,2) - pforceold(1,2) / massvec(1,2) ) / velhalfold(1,2) ) * disp(1,2) )
    else 
        print*, 'velhalfold(1,2) = 0'
    endif
    cn_denom = cn_denom + disp(1,1) * disp(1,1) + disp(1,2) * disp(1,2)
end subroutine compute_ADR_cn

subroutine compute_ADR_Kij(coord_current,coord_other, horizon, idist, a, b, d, Gb, Gd, volume_current, volume_other, Kijx, Kijy)
    ! page 138
    implicit none
    logical echo
    parameter(echo=.FALSE.)
    real*8, intent(in) :: horizon, idist, a, b, d, volume_current, volume_other, Gb, Gd
    real*8 Kijx, Kijy
    real *8 coord_current(1,2),coord_other(1,2) 
    real*8 x_proj, y_proj
    real*8 ex(1,2), ey(1,2)
    ex(1,1) = 1.0
    ex(1,2) = 0.0
    ey(1,1) = 0.0
    ey(1,2) = 1.0
    call dot_product((coord_other(1,:) - coord_current(1,:)),ex,x_proj)
    call dot_product((coord_other(1,:) - coord_current(1,:)),ey,y_proj)
    Kijx = abs(x_proj) * 4 * horizon / idist / idist * ( 0.5 * a * (Gd * d)**2 * horizon / abs(idist) * ( volume_current + volume_other) + Gb * b )
    Kijy = abs(y_proj) * 4 * horizon / idist / idist * ( 0.5 * a * (Gd * d)**2 * horizon / abs(idist) * ( volume_current + volume_other) + Gb * b )
    if (echo) then
        print*, 'Kijx = ',Kijx
        print*, 'Kijy = ',Kijy
    endif
end subroutine compute_ADR_Kij

subroutine preprocess_with_SCF(horizon, delta, volume, d, b, a, coord, disp, numfam, pointfam,nodefam,total_point, max_member,DSCF,SSCF,Theta)
    implicit none
    logical echo
    parameter(echo = .FALSE.)
    real*8, intent(in) :: horizon, delta, volume, d, b, a
    integer i
    real*8 idist, nlength, stretch, fac_vol, applied
    real*8 coord(total_point,2), disp(total_point,2), DSCF(total_point,2), SSCF(total_point,2), Theta(total_point,1)
    integer numfam(total_point,1), pointfam(total_point,1), nodefam(max_member,1)
    integer current_point, other_point, total_point, max_member
    real*8 Lambda, Theta_k, Theta_k_j, Distort_k_j, Distort_k, SED_k
    real*8 Gd, Gb
    do current_point = 1, total_point
        Distort_k = 0.
        Theta_k = 0.
        SED_k = 0.
        do other_point = 1, numfam(current_point,1)
            idist = 0.
            nlength = 0.
            stretch = 0.
            Lambda = 0.
            Theta_k_j = 0.
            Distort_k_j = 0.
            Gd = 0.
            Gb = 0.
            call compute_kinematics(idist, nlength, stretch, Lambda, coord(current_point,:), coord(nodefam(pointfam(current_point,1)+other_point-1,1),:) , disp(current_point,:), disp(nodefam(pointfam(current_point,1)+other_point-1,1),:))
            call compute_fac_volume_corr(fac_vol, idist, horizon, delta)
            call surface_correction_vector(DSCF(current_point,:),SSCF(current_point,:),coord(current_point,:),DSCF(nodefam(pointfam(current_point,1)+other_point-1,1),:),SSCF(nodefam(pointfam(current_point,1)+other_point-1,1),:),coord(nodefam(pointfam(current_point,1)+other_point-1,1),:),idist,Gb,Gd)
            call compute_Theta_k_j(Theta_k_j, stretch, Lambda, volume*fac_vol)
            Theta_k = Theta_k + d * horizon * Gd * Theta_k_j
            call compute_Distort_k_j(Distort_k_j,idist,nlength,volume*fac_vol)
            Distort_k = Distort_k + b * horizon * Gb * Distort_k_j
        enddo 
        SED_k = a * Theta_k **2 + Distort_k
        Theta(current_point,1) = Theta_k
    enddo
end subroutine preprocess_with_SCF
