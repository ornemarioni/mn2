MODULE energia
    IMPLICIT NONE
    CONTAINS
    
    
!Calcula la energia potencial de cada una de las particulas

    SUBROUTINE epot(eps,m,x,y,z,n,ep)
    
    USE OMP_LIB
    
    INTEGER, INTENT (IN):: n
    REAL, INTENT (IN)   :: x(n), y(n), z(n), m(n), eps
    REAL, INTENT (OUT)  :: ep(n)
!f2py INTENT (HIDE)     :: n
!f2py INTENT (IN)       :: x, y, z , m, eps
!f2py INTENT (OUT)      :: ep
    REAL, PARAMETER     :: G = 4.299e-6
    REAL                :: aux, dist, dx, dy, dz
    INTEGER             :: i, j      
    
    
    DO i = 1, n  
       ep(i) = 0.
    END DO
    
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (n,x,y,z,eps,m,ep) &
!$OMP PRIVATE(i,j,dx,dy,dz,dist,aux)
!$OMP DO SCHEDULE(DYNAMIC)
    
    DO i = 1, n
        DO j = 1, n
            dx = x(j)-x(i)
            dy = y(j)-y(i)
            dz = z(j)-z(i)
            dist = sqrt(dx**2 + dy**2 + dz**2 + eps**2)
            
            IF (i /= j) THEN     
                   aux = G*m(i)*m(j) / dist
                   ep(i) = aux + ep(i) 
            END IF
        END DO
    END DO
    
!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE
END MODULE