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
            print*, left_bound, right_bound
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
    parameter(echo = .TRUE.)
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
    stretch = nlength / idist
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
    real-8 Distort_k_j
    
    return
end subroutine compute_Distort_k_j

subroutine preprocess(condition_list, horizon, delta, volume, d, b, idist, nlength, stretch, coord, disp, numfam, pointfam,nodefam,total_point, max_member)
    implicit none
    logical echo
    parameter(echo = .TRUE.)
    real*8, intent(in) :: horizon, delta, volume, d, b
    integer i
    real*8 idist, nlength, stretch, fac_vol, applied
    real*8 coord(total_point,2), disp(total_point,2), bforce(total_point,1)
    integer numfam(total_point,1), pointfam(total_point,1), nodefam(max_member,1)
    integer current_point, other_point,  total_point, max_member
    character(len=60) condition
    character(len=60),dimension(4,1) , intent(in) :: condition_list
    real*8 Lambda, Theta_k, Theta_k_j, Distort_k_j, Distort_k
    do i = 1, 3
        condition = condition_list(i,1)
        applied = 0.001
        call set_conditions(condition, coord, disp, bforce, applied, delta, total_point)
        do current_point = 1, total_point
            do other_point = 1, numfam(current_point,1)
                idist = 0.
                nlength = 0.
                stretch = 0.
                Lambda = 0.
                Theta_k_j = 0.
                Theta_k = 0.
                Distort_k_j = 0.
                Distort_k = 0.
                call compute_kinematics(idist, nlength, stretch, Lambda, coord(current_point,:), coord(nodefam(pointfam(current_point,1)+other_point-1,1),:) , disp(current_point,2), disp(nodefam(pointfam(current_point,1)+other_point-1,1),2))
                call compute_fac_volume_corr(fac_vol, idist, horizon, delta)
                call compute_Theta_k_j(Theta_k_j, stretch, Lambda, volume*fac_vol)
                Theta_k = Theta_k + d * horizon * Theta_k_j
                call compute_Distort_k_j(Distort_k_j,idist,nlength,volume*fac_vol)
                Distort_k = Distort_k + b * horizon * Distort_k_j
                if (echo) then
                    print*, 'x_k:', coord(current_point,:), 'y_k:', coord(current_point,:)+disp(current_point,:)
                    print*, 'x_j:', coord(nodefam(pointfam(current_point,1)+other_point-1,1),:), 'y_j:', coord(nodefam(pointfam(current_point,1)+other_point-1,1),:)+disp(nodefam(pointfam(current_point,1)+other_point-1,1),:)                
                    print*, 'Theta_k_j:', Theta_k_j, 'Theta_k:', Theta_k
                endif
            
            enddo
        enddo
    enddo
    return
end subroutine preprocess
