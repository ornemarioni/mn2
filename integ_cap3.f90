MODULE integradores
    IMPLICIT NONE
    CONTAINS

!------------------------------------------------------------------------
!aceleracion
    SUBROUTINE aceleracion(eps,m,x,y,z,n,ax,ay,az)
    
    USE OMP_LIB
    
    INTEGER, INTENT (IN):: n
    REAL, INTENT (IN)   :: x(n), y(n), z(n), m(n), eps
    REAL, INTENT (OUT)  :: ax(n), ay(n), az(n)
!f2py INTENT (HIDE)     :: n
!f2py INTENT (IN)       :: x, y, z , m, eps
!f2py INTENT (OUT)      :: ax, ay, az
    REAL, PARAMETER     :: G = 6.674e-11
    REAL                :: acx, acy, acz, dx, dy, dz, dist
    INTEGER             :: i, j      
   

    DO i = 1, n  
        ax(i) = 0.
        ay(i) = 0.
        az(i) = 0.
    END DO
      
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (n,x,y,z,eps,m,ax,ay,az) &
!$OMP PRIVATE(i,j,dx,dy,dz,dist,acx,acy,acz)
!$OMP DO SCHEDULE(DYNAMIC)

    DO i = 1, n
    
        DO j = 1, n
            dx = x(j)-x(i)
            dy = y(j)-y(i)
            dz = z(j)-z(i)

            dist = sqrt(dx**2 + dy**2 + dz**2 + eps**2)

            IF (i /= j) THEN     
                acx = G * m(j) * dx / dist**3
                ax(i) = acx + ax(i)

                acy = G * m(j) * dy / dist**3
                ay(i) = acy + ay(i)

                acz = G * m(j) * dz / dist**3
                az(i) = acz + az(i)

            END IF
        END DO
    END DO

!$OMP END DO
!$OMP END PARALLEL

    END SUBROUTINE

!----------------------------------------------------------------------

!Integrador tipo EULER
    SUBROUTINE euler(eps,dt,nit,m,x,y,z,vx,vy,vz,ax,ay,az,n)
    !euler(n,h,nit,m,eps,x,y,z,vx,vy,vz,ax,ay,az)
    
    USE OMP_LIB
      
    INTEGER              :: i,j
    INTEGER, INTENT (IN) :: n, nit
    REAL, INTENT (IN)    :: dt, eps, m(n)
    REAL, INTENT (INOUT) :: x(n), y(n), z(n)
    REAL, INTENT (INOUT) :: vx(n), vy(n), vz(n)
    REAL, INTENT (INOUT) :: ax(n), ay(n), az(n)

!f2py INTENT (HIDE)  :: n
!f2py INTENT(IN,OUT) :: x(n), y(n), z(n)
!f2py INTENT(IN,OUT) :: vx(n), vy(n), vz(n)
!f2py INTENT(IN,OUT) :: ax(n), ay(n), az(n)

   
    DO j = 1,nit

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (eps,dt,nit,m,x,y,z,vx,vy,vz,ax,ay,az,n,j) &
!$OMP PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)

        DO i = 1,n
            vx(i) = vx(i) + ax(i)*dt*3.15e7
            vy(i) = vy(i) + ay(i)*dt*3.15e7
            vz(i) = vz(i) + az(i)*dt*3.15e7
            
            x(i) = x(i) + vx(i)*dt*3.15e7
            y(i) = y(i) + vy(i)*dt*3.15e7
            z(i) = z(i) + vz(i)*dt*3.15e7
        END DO
        
!$OMP END DO
!$OMP END PARALLEL 

        CALL aceleracion(eps,m,x,y,z,n,ax,ay,az)
      
    END DO

    END SUBROUTINE

END MODULE