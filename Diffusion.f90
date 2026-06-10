program cfd
    implicit none
    
    integer :: i,j,m,n,iter
    real (kind=8):: r_x,r_y,L_x,L_y,error,e,a_E,a_W,a_N,a_S,a_P,gamma,rho
    real (kind=8), dimension(42,42) :: T,T_new
    real (kind=8), dimension(40) :: dx
    real (kind=8), dimension(40) :: dy
    real (kind=8), dimension(41) :: delta_x
    real (kind=8), dimension(41) :: delta_y

    
    L_x=10.0     ! Length of domain in x direction
    L_y=10.0     ! Length of domain in y direction

    m=40      ! m = no of divisions in x
    n=40      ! n = no of divisions in y
    iter=0
    
    gamma=10.0

    error=1e-3
    e=0.0

        
    print*, "-------------------------------------------------"
 
    print*," ___________                                    | "
    print*," \\       //    C omputational                  | "    
    print*,"  \\     //     F low                           | "
    print*,"   \\   //      T urbulence                     | "
    print*,"    \\ //       C ombustion                     | "
    print*,"     \*/        L ab                            | "
    print*,"     ***                                        | "
    print*,"    *****                                       | "
    print*,"     ***                                        | "
    print*,"      *                                         | " 

    print*, "-------------------------------------------------"


    print*, "Enter expansion ratio in x direction:"
    read*, r_x     ! r_x = stretching factor or expansion ration in x direction


    print*,"Enter expansion ratio in y direction:"
    read*, r_y     ! r_y = stretching factor or expansion ration in y direction

    !"---------------------------------------------------"

    ! Calculating  dx and dy 

    if (r_x==1 .or. r_y==1) then

        do i=1,m
            dx(i)=L_x/m
        end  do
    end if

    
    if(r_y==1) then

        do i=1,n
            dy(i)=L_y/n
        end do
    end if

    if(r_x>1 .or. r_x<1) then
        
        dx(1)=L_x*((r_x-1)/(r_x**m-1))

        do i=1,m 
            dx(i)=dx(1)*(r_x**(i-1))
        end do
    end if

    if (r_y>1 .or. r_y<1) then

        dy(1)=L_y*((r_y-1)/(r_y**n-1))
        
        do i=1,n
            dy(i)=dy(1)*(r_y**(i-1))
        end do

    end if

    !"---------------------------------------------------"

    ! calculating delta_x and delta_y

    delta_x(1)=dx(1)/2
    delta_x(m+1)=dx(m)/2

    do i=2,m
        delta_x(i)=(dx(i)+dx(i-1))/2
    end do


    delta_y(1)=dy(1)/2
    delta_y(n+1)=dy(n)/2


    !"---------------------------------------------------"

    ! initialize thr T with zeros

    do i=1,m+2
        do j=1,n+2
            T(i,j)=0
        end do
    end do

    ! "---------------------------------------------------"

    ! Boundary conditions 

    do i=1,m+2
        T(1,i)=10
        T(n+2,i)=30
    end do

    do i=1,n+2
        T(i,1)=20
        T(m+2,i)=40
    end do


    ! "---------------------------------------------------"

    do i=2,n
        delta_y(i)=(dy(i)+dy(i-1))/2
    end do

    do while (.True.) 
        e=0.0
        
        do i=2,m+1
            do j=2,n+1

                a_E=(gamma*dy(j))/delta_x(i)
                a_W=(gamma*dy(j))/delta_x(i-1)
                a_N=(gamma*dx(i))/delta_y(j)
                a_S=(gamma*dx(i))/delta_y(j-1)
                a_P=(a_E+a_W+a_N+a_S)

                T_new(i,j)=(a_E*T(i,j+1)+a_W*T(i,j-1)+a_N*T(i+1,j)+a_S*T(i-1,j))/a_P

                e=e+abs(a_E*T(i,j+1)+a_W*T(i,j-1)+a_N*T(i+1,j)+a_S*T(i-1,j)-a_P*T(i,j))

            end do
        end do

        do i=1,m+2
            T_new(1,i)=10
            T_new(n+2,i)=30
        end do

        do i=1,n+2
            T_new(i,1)=20
            T_new(m+2,i)=40
        end do

        T=T_new

        print*,e

        iter=iter+1

        if (e<error) then 
            exit
        end if

    end do

    print*,"No of iteration required for convergence:",iter
    ! ---------------------------------------------------

end program cfd

