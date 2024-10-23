subroutine normalizedSpatialHet(SH, array)
    implicit none
    integer :: y1, y2, x1, x2, y_end, x_end
    integer :: count
    real(kind=8), intent(out) :: SH ! note that kind=8 for double precision is specific to gfortran!!
    real(kind=8), dimension(:, :), intent(in) :: array
    real(kind=8), dimension(:, :), allocatable :: hstack_array
    real(kind=8), dimension(:, :), allocatable :: v_hstack_array
    real(kind=8) :: array_mean, subarray_mean
    real(kind=8) :: G
    integer, dimension(2) :: dim
    integer :: n, m
    G = 0

    dim = shape(array)
    n = dim(1)
    m = dim(2)

    allocate(hstack_array(n, 2*m))
    allocate(v_hstack_array(2*n, 2*m))
    
    array_mean = sum(array)/max(1, size(array))
    ! reshape does horizontal stacking, so vertical stack requires transpose, horiz stack, transpose
    hstack_array = reshape([array, array], [n, 2*m])
    v_hstack_array = transpose(reshape([transpose(hstack_array), transpose(hstack_array)], [2*n, 2*m]))

    do y1 = 1, n
        do y2 = 1, n
            do x1 = 1, m
                do x2 = 1, m

                    y_end = y1 + y2-1
                    x_end = x1 + x2-1

                    subarray_mean = 0.0
                    count = (y_end - y1 + 1) * (x_end - x1 + 1)
                    subarray_mean = sum(v_hstack_array(y1:y_end, x1:x_end)) / count
                    G = G + abs(subarray_mean - array_mean)

                end do
            end do
        end do
    end do

    SH = G*2/(array_mean*n*m*(3*n*m-n-m-1))
    deallocate(hstack_array)
    deallocate(v_hstack_array)
end subroutine normalizedSpatialHet


! subroutine to calculate spatial heterogeneity using monte carlo sampling
subroutine monteCarloSpatialHet(SH, array, n_permutes)
    implicit none
    integer :: y1, y2, x1, x2
    integer :: i, j, k, count
    real(kind=8), intent(out) :: SH ! note that kind=8 for double precision is specific to gfortran!!
    integer, intent(in) :: n_permutes
    real(kind=8), dimension(:, :), intent(in) :: array
    real(kind=8), dimension(:, :), allocatable :: hstack_array
    real(kind=8), dimension(:, :), allocatable :: v_hstack_array
    integer, dimension(:, :), allocatable :: permutation_array
    real(kind=8) :: array_mean, subarray_mean
    real(kind=8) :: G
    integer, dimension(2) :: dim
    integer :: d1, d2
    integer :: total_permutes
    real(kind=8) :: scaling_factor

    G = 0
    dim = shape(array)
    d1 = dim(1)
    d2 = dim(2)
    total_permutes = d1*d1*d2*d2

    if (n_permutes > total_permutes) then
        print *, "Number of permutations is greater than total possible permutations. Exiting."
        return
    end if

    allocate(hstack_array(d1, 2*d2))
    allocate(v_hstack_array(2*d1, 2*d2))
    allocate(permutation_array(4, n_permutes))
    
    array_mean = sum(array)/max(1, size(array))
    ! reshape does horizontal stacking, so vertical stack requires transpose, horiz stack, transpose
    hstack_array = reshape([array, array], [d1, 2*d2])
    v_hstack_array = transpose(reshape([transpose(hstack_array), transpose(hstack_array)], [2*d1, 2*d2]))

    do i = 1, n_permutes
        call rectangleBounds(x1, x2, 1, d1)
        call rectangleBounds(y1, y2, 1, d2)

        subarray_mean = 0.0
        count = 0
        do k = y1, y2
            do j = x1, x2
                subarray_mean = subarray_mean + v_hstack_array(j, k)
                count = count + 1
            end do
        end do
        subarray_mean = subarray_mean / count

        G = G + abs(subarray_mean - array_mean)
        
    end do

    SH = G*2/(array_mean*d1*d2*(3*d1*d2-d1-d2-1))
    scaling_factor = real(total_permutes, kind=8)/real(n_permutes, kind=8)
    SH = SH*scaling_factor

    deallocate(hstack_array)
    deallocate(v_hstack_array)

end subroutine monteCarloSpatialHet

subroutine rectangleBounds(l_bound, r_bound, min, max)
    ! note that the right bound for normalized spatial heterogeneity is
    ! independent of the left bound (i.e., x1 or y1)
    implicit none 
    integer :: l_bound, r_bound, min, max
    
    call randint(l_bound, min, max)
    call randint(r_bound, min, max)

    if (r_bound < l_bound) then 
        r_bound = r_bound + max
    end if

end subroutine rectangleBounds

subroutine randint(j, min, max)
    implicit none
    real :: r
    integer, intent(out) :: j 
    integer, intent(in) :: min, max
    call random_number(r) ! float value between 0 and 1
    !j = min + FLOOR((max-min)*r) ! scale value by min max range and round to int

    j = min + FLOOR((1+max-min)*r) ! scale value by min max range and round to int
end subroutine randint
    