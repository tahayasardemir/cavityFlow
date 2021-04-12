	  program CavityFlow
c..Taha Yaşar Demir / 1881978
c..CE-580 HomeWork #7
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx) 
	  common/tho/ a(mx),b(mx),c(mx),d(mx)
	  real tolerance
	  integer iter
	  open(22,file="erroradi.dat")
	  call grid
	  call init
	  Ers = 1.
	  Erv = 1.
	  iter= 0.
	  tolerance = 1e-7
	  do while(Erv.gt.tolerance)
c	  do m=1,25000
	  	iter = iter + 1
	  	call boundary
	  	call evalcoef
	  	call ADI
	  	call psor
	  	call velocity
	  	call error(iter)
c	  	print*, iter,Erv,Ers
	  enddo
	  call output
	  call dragcal
	  close(22)
	  close(33)
	  stop
	  end

c-----------------------------------------------------------------------
	  subroutine grid
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

	  rL   = 0.01 !m
	  N    = 101
	  dydx = rL/(N-1) 

	  do j=1,N 
	  	x(1,j) = 0.
	  	do i=2,N-1
	  		x(i,j) = x(i-1,j) + dydx
	  	enddo
	  	x(N,j) = rL
	  enddo


	  do i=1,N
	 	y(i,1) = 0.
	 	do j=2,N-1
	 		y(i,j) = y(i,j-1) + dydx
	 	enddo
	 	y(i,N) = rL
	  enddo


	  return
	  end
c-----------------------------------------------------------------------
	  subroutine init
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

	  time = 0.
	  u_p  = 0.01 ! 0.01-0.02-0.05-0.1-1
	  visc = 1E-6
	  dt   = 0.1*(dydx/u_p)
	  call boundary
	  do i=2,N-1
	  	do j=2,N-1
	  		u(i,j)   = 0.
	  		v(i,j)   = 0.
	  		psi(i,j) = 0.
	  		zeta(i,j)= (v(i,j)-v(i-1,j))/dydx - (u(i,j)-u(i-1,j))/dydx
	  	enddo
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine boundary
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

 	  do i=2,N-1
 	  	psi(i,1) = 0.
 	  	psi(i,N) = 0.
 	  	u(i,1)   = 0.
 	  	u(i,N)   = u_p
 	  	v(i,1)   = 0.
 	  	v(i,N)   = 0.
 	  	zeta(i,1)= 2*(psi(i,1)-psi(i,2))/(dydx**2)
 	  	zeta(i,N)= (2*(psi(i,N)-psi(i,N-1))/(dydx**2)) - 2*u_p/dydx
 	  enddo

 	  do j=1,N
 	  	psi(1,j) = 0.
 	  	psi(N,j) = 0.
 	  	u(1,j)   = 0.
 	  	u(N,j)   = 0.
 	  	v(1,j)   = 0.
 	  	v(N,j)   = 0.
 	  	zeta(1,j)= 2*(psi(1,j)-psi(2,j))/(dydx**2)
 	  	zeta(N,j)= 2*(psi(N,j)-psi(N-1,j))/(dydx**2)
 	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine evalcoef
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

	  do i=1,N
	  	do j=1,N
	  		if (u(i,j).gt.0.) then
	  			epsx(i,j) = 1.
	  		else
	  			epsx(i,j) =-1.
	  		endif
	  		if (v(i,j).gt.0.) then
	  			epsy(i,j) = 1.
	  		else
	  			epsy(i,j) =-1.
	  		endif
	  		crx(i,j) = u(i,j)*dt/dydx
	  		cry(i,j) = v(i,j)*dt/dydx
	  		dx       = visc*dt/(dydx**2)
	  		dy       = visc*dt/(dydx**2)
	  	enddo
	  enddo
	  return
	  end
c-----------------------------------------------------------------------
	  subroutine ADI
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

c...x_sweep 
	  do j=2,N-1
	  	call x_sweep(j)
	  enddo

c...y_sweep
	  do i=2,N-1
	  	call y_sweep(i)
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine x_sweep(j)
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)


	  do i=2,N-1
	  	if(i.eq.2) then
	  		a(i) =  0.
	  		b(i) =  1 + dx + 0.5*epsx(i,j)*crx(i,j)
	  		c(i) = -0.5*(dx - 0.5*(1-epsx(i,j))*crx(i+1,j))
	  		d(i) =  0.5*(dx + 0.5*(1+epsx(i,j))*crx(i-1,j))*zeta(i-1,j)
     &            + 0.5*(dy - 0.5*(1-epsy(i,j))*cry(i,j+1))*zeta(i,j+1)
     &            + (1-dy-0.5*epsy(i,j)*cry(i,j))*zeta(i,j)
     &            + 0.5*(dy + 0.5*(1+epsy(i,j))*cry(i,j-1))*zeta(i,j-1) 
	  	elseif(i.eq.N-1) then 
	  		a(i) = -0.5*(dx + 0.5*(1+epsx(i,j))*crx(i-1,j))
	  		b(i) =  1 + dx + 0.5*epsx(i,j)*crx(i,j)
	  		c(i) =  0.
	  		d(i) =  0.5*(dx - 0.5*(1-epsx(i,j))*crx(i+1,j))*zeta(i+1,j)
     &            + 0.5*(dy - 0.5*(1-epsy(i,j))*cry(i,j+1))*zeta(i,j+1)
     &            + (1-dy-0.5*epsy(i,j)*cry(i,j))*zeta(i,j)
     &            + 0.5*(dy + 0.5*(1+epsy(i,j))*cry(i,j-1))*zeta(i,j-1) 
	  	else
	  		a(i) = -0.5*(dx + 0.5*(1+epsx(i,j))*crx(i-1,j))
			b(i) =  1 + dx + 0.5*epsx(i,j)*crx(i,j)
	  		c(i) = -0.5*(dx - 0.5*(1-epsx(i,j))*crx(i+1,j))
			d(i) =  0.5*(dy - 0.5*(1-epsy(i,j))*cry(i,j+1))*zeta(i,j+1)
     &            + (1-dy-0.5*epsy(i,j)*cry(i,j))*zeta(i,j)
     &            + 0.5*(dy + 0.5*(1+epsy(i,j))*cry(i,j-1))*zeta(i,j-1)

	  	endif
	  enddo

	  call THOMAS(2,N-1,a,b,c,d)

	  do i=2,N-1 ! extract the solution from thomas algorithm
	  	zeta(i,j) = d(i)
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine y_sweep(i)
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

	  do j=2,N-1
	  	if(j.eq.2) then
	  		a(j) =  0.
	  		b(j) =  1 + dy + 0.5*(epsy(i,j))*cry(i,j)
	  		c(j) = -0.5*(dy - 0.5*(1-epsy(i,j))*cry(i,j+1))
	  		d(j) =  0.5*(dy + 0.5*(1+epsy(i,j))*cry(i,j-1))*zeta(i,j-1)
     &             +0.5*(dx - 0.5*(1-epsx(i,j))*crx(i+1,j))*zeta(i+1,j)
     &             + (1-dx-0.5*epsx(i,j)*crx(i,j))*zeta(i,j)
     &             +0.5*(dx + 0.5*(1+epsx(i,j))*crx(i-1,j))*zeta(i-1,j)
	  	elseif(j.eq.N-1) then
	  		a(j) = -0.5*(dy + 0.5*(1+epsy(i,j))*cry(i,j-1))
	  		b(j) =  1 + dy + 0.5*(epsy(i,j))*cry(i,j)
	  		c(j) =  0.
	  		d(j) =  0.5*(dy - 0.5*(1-epsy(i,j))*cry(i,j+1))*zeta(i,j+1)
     &             +0.5*(dx - 0.5*(1-epsx(i,j))*crx(i+1,j))*zeta(i+1,j)
     &             + (1-dx-0.5*epsx(i,j)*crx(i,j))*zeta(i,j)
     &             +0.5*(dx + 0.5*(1+epsx(i,j))*crx(i-1,j))*zeta(i-1,j)
	  	else
	  		a(j) = -0.5*(dy + 0.5*(1+epsy(i,j))*cry(i,j-1))
	  		b(j) =  1 + dy + 0.5*(epsy(i,j))*cry(i,j)
	  		c(j) = -0.5*(dy - 0.5*(1-epsy(i,j))*cry(i,j+1))
	  		d(j) =  0.5*(dx - 0.5*(1-epsx(i,j))*crx(i+1,j))*zeta(i+1,j)
     &             + (1-dx-0.5*epsx(i,j)*crx(i,j))*zeta(i,j)
     &             +0.5*(dx + 0.5*(1+epsx(i,j))*crx(i-1,j))*zeta(i-1,j)
	  	endif
	  enddo

	  call THOMAS(2,N-1,a,b,c,d)

	  do j=2,N-1
	  	zeta(i,j) = d(j)
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine psor 
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)
	  real omega,sum
	  open(33,file="errorpsor.dat")
	  omega = 1.8 ! over-relaxation parameter
	  sum   = 0.
	  do k=1,20
	  do j=2,N-1
	  	do i=2,N-1
	  		R(i,j) = 0.25*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)
     &                     - 4*psi(i,j) + (dydx**2)*zeta(i,j))
	  		psi(i,j) = psi(i,j) + omega*R(i,j)
	  		sum = sum + abs(R(i,j))
	  	enddo
	  enddo
	  sum = sum/((N-2)**2)
	  write(33,*) k,sum
	  enddo


	  return
	  end
c-----------------------------------------------------------------------
	  subroutine velocity
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

	  do i=2,N-1
	  	do j=2,N-1
	  		u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dydx)
	  		v(i,j) =-(psi(i+1,j) - psi(i-1,j))/(2*dydx)
	  	enddo
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine error(iteration)
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)
	  real adier , psorer, zetaold(mx,mx), psiold(mx,mx), sum
	  integer iteration

	  adier = 0.
	  psorer= 0.
	  sum   = 0.
	  do i=1,N
	  	sum = sum + abs(zeta(i,N))
	  enddo
	  sum = sum/N
c.. L2 normalization is used for vorticity and stream function values
	  if(iteration.eq.1) then ! store the previous zeta and psi values
	  	do j=1,N
	  		do i=1,N
	  			zetaold(i,j) = zeta(i,j)
	  			psiold(i,j)  = psi(i,j)
	  		enddo
	  	enddo
	  	print*, 
	  else
	  	do j=1,N
	  		do i=1,N
	  		adier=adier+abs((zeta(i,j)-zetaold(i,j))/sum)
	  		zetaold(i,j) = zeta(i,j) ! update the old values for next iteration
	  		psorer       = psorer+abs(psi(i,j)-psiold(i,j))
	  		psiold(i,j)  = psi(i,j)
	  		enddo
	  	enddo

	  	Erv = adier/(N**2) ! Vorticity transport equation error
	  	Ers = psorer/(N**2)! Stream function solution errror

	  	write(22,*) iteration,Erv,Ers

	  endif

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine output
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx)
	  common/tho/ a(mx),b(mx),c(mx),d(mx)

	  open(11,file='var.tec',form='formatted')
	  write(11,*) ' variables="x","y","zeta","psi","u","v" '
	  write(11,*) ' zone i=',N, 'j=',N
	  do j=1,N
	  	do i=1,N
	  		write(11,'(8E12.4)') x(i,j),y(i,j),zeta(i,j),psi(i,j),
     +                          u(i,j),v(i,j)
	  	enddo
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine dragcal
	  parameter (mx=101)
	  common/grd/ x(mx,mx),y(mx,mx),dydx,N,rL,epsx(mx,mx),epsy(mx,mx)
	  common/flw/ u(mx,mx),v(mx,mx),psi(mx,mx),zeta(mx,mx)
	  common/par/ u_p, visc, dt, time, crx(mx,mx), cry(mx,mx),dx,dy
	  common/err/ Erv, Ers, R(mx,mx) 
	  common/tho/ a(mx),b(mx),c(mx),d(mx)
	  real mu,rho,shear,drag

	  rho   = 1000.
	  mu    = visc * rho
	  drag  = 0.
	  do i=2,N
	  	shear = mu*(3*u(i,N)-4*u(i,N-1)+u(i,N-2))/(2*dydx)
	 	drag  = drag + shear*dydx
	  enddo
	  print*, "2-D drag force = ", drag
	  return
	  end
c-----------------------------------------------------------------------
      subroutine THOMAS(il,iu,aa,bb,cc,ff)
c............................................................
c  Solution of a tridiagonal system of n equations of the form
c  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = R(i)  for i=il,iu
c  the solution X(i) is stored in F(i)
c  A(il-1) and C(iu+1) are not used
c  A,Bb,C,R are arrays to bbe provided bby the user
c............................................................
      parameter (mx=101)
      dimension aa(mx),bb(mx),cc(mx),ff(mx),tmp(mx)

      tmp(il)=cc(il)/bb(il)
      ff(il)=ff(il)/bb(il)
      ilp1 = il+1
      do i=ilp1,iu
         z=1./(bb(i)-aa(i)*tmp(i-1))
         tmp(i)=cc(i)*Z
         ff(i)=(ff(i)-aa(i)*ff(i-1))*z
      enddo
      iupil=iu+il
      do ii=ilp1,iu
         i=iupil-ii
         ff(i)=ff(i)-tmp(i)*ff(i+1)
      enddo
      return
      end