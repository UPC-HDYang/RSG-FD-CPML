!##############################################
!##############################################
!######                                ########
!######       PROGRAM RCPF.f90         ########
!######       Author: Haidi Yang       ########
!######       Time: 04/25/2020         ########
!######       UPC                      ########
!######       Ver 1.1-2                ########
!######                                ########
!##############################################
!##############################################	

program RCPF
implicit none
! RCPF : R - RSG, CP - CPML, F - Free surface
! 2D biot rotated-staggered-grid finite-difference code in velocity-stress formulation
! with convolution PML (CPML) absorbing conditions for an isotropic medium.

! For stress at grid point, and velocity at half grid

! x8t2
! for Free surface, it should be processed with specially attention besides setting the velocity to zero, density to small number close to zero.

! in isotropic medium, averaging the parameter is not needed. But in inhomogeneous medium, should averaging the parameter.

    
	character*77 parfile
	character*77 fname	
    integer*4 nx,nz,nt,xs,zs,np,iw,dst,iwav 
	real*4 dx,dz,dt,vmax,r,chi0,a0,f0,ew,power_display,cutvect 
	integer*4 nr,fxr,fzr,dxr,dzr 
	integer*4 ifs,fslc 
	integer*4 ifa,sl
    real*4 pi
	real*4 c1,c2,c3,c4    
    real*4 dr,maxpar
	real*4 d0,alphamax
    real*4 ts,refdist,Courant_number
    integer*4 ix,iz,it,i
	
	real*4,allocatable,dimension(:,:) :: rho
	real*4,allocatable,dimension(:,:) :: c11,c12,c13,c22,c23,c33
    real*4,allocatable,dimension(:,:) :: rhoxz
    real*4,allocatable,dimension(:,:) :: txx,txz,tzz,vx,vz
	real*4,allocatable,dimension(:,:) :: dxtxx,dztxx,dxtxz,dztxz,dxtzz,dztzz 
    real*4,allocatable,dimension(:,:) :: dxvx,dzvz,dxvz,dzvx
!   parameter for CPML
	real*4,allocatable,dimension(:,:) :: ax,bx,az,bz,chix,chiz,cdx,cdz 
	real*4,allocatable,dimension(:,:) :: ax_half,bx_half,az_half,bz_half,chix_half,chiz_half,cdx_half,cdz_half 
    real*4,allocatable,dimension(:,:) :: psidxtxx,psidztzz,psidxtxz,psidztxz 
    real*4,allocatable,dimension(:,:) :: psidxvx,psidzvx,psidxvz,psidzvz
! output velocity and stress
    real*4,allocatable,dimension(:,:) :: seisvx,seisvz 
    real*4,allocatable,dimension(:,:) :: seistxx,seistzz,seistxz
    real*4,allocatable,dimension(:,:) :: wav 
    real*4,allocatable,dimension(:) :: wavelet,wlet 
	integer :: status=0 

! cpu_time
	real time_begin,time_end 
    
!##########################################################################################################	
!##########################################################################################################
	write(*,*) 'Input the parfile name:'
	read(*,*) parfile 
	print *,parfile 
	
	open(99,file=parfile,form='formatted',access='sequential',iostat=status) 
	write(*,*) status
	read(99,'(64x)')
	
	read(99,'(70x)')
	read(99,*) nx,nz,nt,dx,dz,dt 
	write(*,*) nx,nz,nt,dx,dz,dt
	
	read(99,'(70x)')
	read(99,*) np,vmax,r,chi0,a0 
	write(*,*) np,vmax,r,chi0,a0
	
	read(99,'(70x)')
	read(99,*) f0,iw,dst,iwav,power_display,cutvect 
	write(*,*) f0,iw,dst,iwav,power_display,cutvect
	
	read(99,'(70x)')
	read(99,*) xs,zs 
	write(*,*) xs,zs
	
	read(99,'(70x)')
	read(99,*) nr,fxr,fzr,dxr,dzr 
	write(*,*) nr,fxr,fzr,dxr,dzr
	
	read(99,'(70x)')
	read(99,*) ifs,fslc,ifa
	write(*,*) ifs,fslc,ifa 	

        read(99,'(70x)')
	read(99,*) sl
	write(*,*) sl 
	close(99)		
 
    dr = sqrt(dx**2+dz**2)
    dx = 2.0*dx 
    dz = 2.0*dz 
	!   modify the coefficient of spatial opertor    
    c1 = (1225.0/1024.0)
    c2 = c1/15.0 
    c3 = c1/125.0 
    c4 = c1/1715.0 
    ts = 1.2/f0 
	
	pi = 4.0*atan(1.0) 
    d0 = -3*vmax*log(r)/(2.0*np*dx/2.0) 
    alphamax = pi*a0 
	
!	allocate the arrays
    allocate(rho(nz,nx),c11(nz,nx),c12(nz,nx),c13(nz,nz),c22(nz,nx),c23(nz,nx),c33(nz,nx))
	allocate(txx(nz,nx),txz(nz,nx),tzz(nz,nx),vx(nz,nx),vz(nz,nx))
	allocate(rhoxz(nz,nx))
	allocate(dxtxx(nz,nx),dztxx(nz,nx),dxtxz(nz,nx),dztxz(nz,nx),dxtzz(nz,nx),dztzz(nz,nx))
	allocate(dxvx(nz,nx),dzvz(nz,nx),dxvz(nz,nx),dzvx(nz,nx))
	allocate(ax(nz,nx),bx(nz,nx),az(nz,nx),bz(nz,nx),chix(nz,nx),chiz(nz,nx),cdx(nz,nx),cdz(nz,nx))
	allocate(ax_half(nz,nx),bx_half(nz,nx),az_half(nz,nx),bz_half(nz,nx),chix_half(nz,nx),chiz_half(nz,nx),cdx_half(nz,nx),cdz_half(nz,nx))
	allocate(psidxtxx(nz,nx),psidztzz(nz,nx),psidxtxz(nz,nx),psidztxz(nz,nx))
	allocate(psidxvx(nz,nx),psidzvx(nz,nx),psidxvz(nz,nx),psidzvz(nz,nx))
	allocate(wav(nz,nx))
	allocate(seisvx(nt,nr),seisvz(nt,nr),seistxx(nt,nr),seistzz(nt,nr),seistxz(nt,nr))
	allocate(wavelet(nt))
	allocate(wlet(nt))
	

	open(1,file='rho.dat',access='stream')
    open(2,file='c11.dat',access='stream')
    open(3,file='c12.dat',access='stream')
    open(4,file='c13.dat',access='stream')
    open(5,file='c22.dat',access='stream')
    open(6,file='c23.dat',access='stream')
	open(7,file='c33.dat',access='stream')
    open(8,file='wlet.dat',access='stream')
	read(1)rho
    read(2)c11
	read(3)c12
	read(4)c13
	read(5)c22
	read(6)c23
	read(7)c33
	read(8)wlet
	
	rho=rho*1.0e3
        c11=c11*1.0e9
	c12=c12*1.0e9
	c13=c13*1.0e9
	c22=c22*1.0e9
	c23=c23*1.0e9
	c33=c33*1.0e9
    wlet=wlet
	
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(6)
	close(7)
	close(8)
	if(ifs .eq. 1) then
	    rho(1:fslc-1,:) = rho(1:fslc-1,:)*1.0e-3
        rho(nz:nz,:) = rho(nz:nz,:)*1.0e-3
       c11(1:fslc-1,:) = 0.0
        c12(1:fslc-1,:) = 0.0
        c13(1:fslc-1,:) = 0.0
		c22(1:fslc-1,:) = 0.0
		c23(1:fslc-1,:) = 0.0
		c33(1:fslc-1,:) = 0.0
		c11(nz:nz,:) = 0.0
        c12(nz:nz,:) = 0.0
		c13(nz:nz,:) = 0.0
		c22(nz:nz,:) = 0.0
		c23(nz:nz,:) = 0.0
		c33(nz:nz,:) = 0.0
    endif
	if(ifa .eq. 1) then
	rho(:,1:fslc-1) = rho(:,1:fslc-1)*1.0e-3
    rho(:,nx:nx) = rho(:,nx:nx)*1.0e-3

       c11(:,1:fslc-1)=0.0
       c12(:,1:fslc-1)=0.0
	   c13(:,1:fslc-1)=0.0
       c22(:,1:fslc-1)=0.0
	   c23(:,1:fslc-1)=0.0
	   c33(:,1:fslc-1)=0.0
		c11(:,nx:nx) = 0.0
		c12(:,nx:nx) = 0.0
		c13(:,nx:nx) = 0.0
		c22(:,nx:nx) = 0.0
		c23(:,nx:nx) = 0.0
		c33(:,nx:nx) = 0.0
    endif
 


! interpolate rho and b at half grid
    do iz=1,nz-1 
        do ix=1,nx-1
        rhoxz(iz,ix) = 0.25*(rho(iz  ,ix  )+rho(iz+1,ix  )+rho(iz+1,ix+1)+rho(iz  ,ix+1))
		enddo
    enddo
	rhoxz(nz  ,1:nx-1) = 0.5*(rho(nz  ,1:nx-1)+rho(nz  ,2:nx  ))
    rhoxz(1:nz-1,nx) = 0.5*(rho(1:nz-1  ,nx  )+rho(2:nz  ,nx  ))
    rhoxz(nz,nx) = rho(nz,nx)
    
	if(iw .eq. 1) then
		call dgauss(wavelet,nt,dt,f0,ts)
	elseif(iw .eq. 0) then
		call ricker(wavelet,nt,dt,f0,ts)
	elseif(iw.eq. 2) then
                wavelet=wlet
	else
		wavelet = 0.0
		wavelet(1) = 1.0
    endif

    open(55,file='wavelet.dat',access='stream')
    write(55)wavelet
    close(55)
	if(iwav .eq. 1) then
		do ix = max(1,xs-3),min(nx,xs+3)
			do iz = max(1,zs-3),min(nz,zs+3)
				wav(iz,ix) = exp(-0.2*((real(iz)-zs)**2 + (real(ix)-xs)**2))
			end do
		end do
	else
		wav = 0.0
		wav(zs-2:zs+2,xs-2:xs+2) = 1.0
	endif
	open(56,file='wav.dat',access='stream')
    write(56)wav
    close(56)

! time begin
	call cpu_time(time_begin)

!---------------PML------------------------
! need vmax, so compute vmax_max first. Here assume vmax is known.
    
    ax = 0.0
    az = 0.0
    bx = 1.0
    bz = 1.0
    chix = 1.0
    chiz = 1.0
    cdx = 0.0
    cdz = 0.0

    ax_half = 0.0
    az_half = 0.0
    bx_half = 1.0
    bz_half = 1.0
    chix_half = 1.0
    chiz_half = 1.0
    cdx_half = 0.0
    cdz_half = 0.0

	do ix = 1,nx
	    if(ifs .eq. 0) then
	         do iz = 1,np+5
	              refdist = (np+5.0-(iz))/np
	              if(refdist > 1.0) refdist = 1.0
	              if(refdist < 0.0) refdist = 0.0
	              if(refdist >= 0.0) then
	                  chiz(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	                  cdz(iz,ix) = d0*(refdist)**2
	                  bz(iz,ix) = exp(-(cdz(iz,ix)/chiz(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	                  az(iz,ix) = cdz(iz,ix) / (chiz(iz,ix)*(cdz(iz,ix) + chiz(iz,ix)*alphamax*(1.0-refdist))) * (bz(iz,ix) - 1.0)
	              endif
	              if(ix == 1) print *,'refdist at iz=',iz,'=',refdist

	              refdist = ((np+4.5-(iz))/np)
	              if(refdist > 1.0) refdist = 1.0
	              if(refdist < 0.0) refdist = 0.0
	              if(refdist >= 0.0) then
	                   chiz_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	                   cdz_half(iz,ix) = d0*(refdist)**2
	                   bz_half(iz,ix) = exp(-(cdz_half(iz,ix)/chiz_half(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	                   az_half(iz,ix) = cdz_half(iz,ix) / (chiz_half(iz,ix)*(cdz_half(iz,ix) + chiz_half(iz,ix)*alphamax*(1.0-refdist))) * (bz_half(iz,ix) - 1.0)
	               endif
	               if(ix == 1) print *,'refdist half at iz=',iz,'=',refdist
	         enddo
        endif

     do iz = (nz-4-np),nz
          refdist = ((iz-0.5)-(nz-4.0-np))/np
         if(refdist > 1.0) refdist = 1.0
          if(refdist < 0.0) refdist = 0.0
          if(refdist >= 0.0) then
              chiz(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
              cdz(iz,ix) = d0*(refdist)**2
              bz(iz,ix) = exp(-(cdz(iz,ix)/chiz(iz,ix) + alphamax*(1.0-refdist)) * dt)
              az(iz,ix) = cdz(iz,ix) / (chiz(iz,ix)*(cdz(iz,ix) + chiz(iz,ix)*alphamax*(1.0-refdist))) * (bz(iz,ix) - 1.0)
         endif
         if(ix == 1) print *,'refdist at iz=',iz,'=',refdist

          refdist = (((iz-0.0)-(nz-4.0-np))/np)
         if(refdist > 1.0) refdist = 1.0
          if(refdist < 0.0) refdist = 0.0
          if(refdist >= 0.0) then
               chiz_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
               cdz_half(iz,ix) = d0*(refdist)**2
               bz_half(iz,ix) = exp(-(cdz_half(iz,ix)/chiz_half(iz,ix) + alphamax*(1.0-refdist)) * dt)
               az_half(iz,ix) = cdz_half(iz,ix) / (chiz_half(iz,ix)*(cdz_half(iz,ix) + chiz_half(iz,ix)*alphamax*(1.0-refdist))) * (bz_half(iz,ix) - 1.0)
          endif
          if(ix == 1) print *,'refdist half at iz=',iz,'=',refdist
     enddo
    enddo

	do iz = 1,nz
    if(ifa .eq. 0) then
	 do ix = 1,np+5
	  refdist = (np+5.0-(ix))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	  chix(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	  cdx(iz,ix) = d0*(refdist)**2
	  bx(iz,ix) = exp(-(cdx(iz,ix)/chix(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	  ax(iz,ix) = cdx(iz,ix) / (chix(iz,ix)*(cdx(iz,ix) + chix(iz,ix)*alphamax*(1.0-refdist))) * (bx(iz,ix) - 1.0)
	 endif
	 if(iz == 1) print *,'refdist at ix=',ix,'=',refdist

	  refdist = (np+4.5-(ix))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	   chix_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	   cdx_half(iz,ix) = d0*(refdist)**2
	   bx_half(iz,ix) = exp(-(cdx_half(iz,ix)/chix_half(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	   ax_half(iz,ix) = cdx_half(iz,ix) / (chix_half(iz,ix)*(cdx_half(iz,ix) + chix_half(iz,ix)*alphamax*(1.0-refdist))) * (bx_half(iz,ix) - 1.0)
	  endif
	  if(iz == 1) print *,'refdist half at ix=',ix,'=',refdist
	 enddo

	 do ix = (nx-4-np),nx
	  refdist = ((ix-0.5)-(nx-4-np))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	  chix(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	  cdx(iz,ix) = d0*(refdist)**2
	  bx(iz,ix) = exp(-(cdx(iz,ix)/chix(iz,ix) + alphamax*(1.0-refdist)) * dt)
	  ax(iz,ix) = cdx(iz,ix) / (chix(iz,ix)*(cdx(iz,ix) + chix(iz,ix)*alphamax*(1.0-refdist))) * (bx(iz,ix) - 1.0)
	 endif
	 if(iz == 1) print *,'refdist at ix=',ix,'=',refdist

	  refdist = ((ix+0.0)-(nx-4-np))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	   chix_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	   cdx_half(iz,ix) = d0*(refdist)**2
	   bx_half(iz,ix) = exp(-(cdx_half(iz,ix)/chix_half(iz,ix) + alphamax*(1.0-refdist)) * dt)
	   ax_half(iz,ix) = cdx_half(iz,ix) / (chix_half(iz,ix)*(cdx_half(iz,ix) + chix_half(iz,ix)*alphamax*(1.0-refdist))) * (bx_half(iz,ix) - 1.0)
	  endif
	  if(iz == 1) print *,'refdist half at ix=',ix,'=',refdist

	 enddo
	 end if
	enddo

! print position of the source
	print *,'Position of the source:'
	print *
	print *,'x = ',xs
	print *,'z = ',zs
	print *

! define location of receivers
	print *,'There are ',nr,' receivers'
	print *

	print *,'the 1st receiver: (x_target,z_target) = ',fxr,fzr
	print *,'the last receiver: (x_target,z_target) = ',fxr+real(nr-1)*dxr,fzr+real(nr-1)*dzr
	print *,'distance between receivers: (dx,dz) = ',dxr,dzr

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
	Courant_number = vmax * dt / dr
	print *,'Courant number is ',Courant_number
	print *
	if(Courant_number > 1.0) then
		pause 'time step is too large, simulation will be unstable'
	endif

!-----------------------------------------------------------------------
!                   compute code
!-----------------------------------------------------------------------
! initial the initial values
	psidxtxx = 0.0
    psidztzz = 0.0
    psidxtxz = 0.0
    psidztxz = 0.0
	
    psidxvx = 0.0
    psidzvx = 0.0
    psidxvz = 0.0
    psidzvz = 0.0
	
	
    vx = 0.0
    vz = 0.0
	
	txx = 0.0
	txz = 0.0
	tzz = 0.0


	write(*,*) 'Begin time step computing ...'
	

	do it = 1,nt
        
        if(mod(it,100) .eq. 0) print *,'No. ',it,' of ',nt


!   compute the velocity of nonstiff equations
!   compute the spatial difference of stress and pressure
        dxtxx = 0.0
        dztxx = 0.0
        dxtxz = 0.0
        dztxz = 0.0
        dxtzz = 0.0
        dztzz = 0.0

        
        do ix=4,nx-4
            do iz=4,nz-4
                dxtxx(iz,ix) = (c1*(txx(iz  ,ix+1)-txx(iz+1,ix  )) - c2*(txx(iz-1,ix+2)-txx(iz+2,ix-1)) + c3*(txx(iz-2,ix+3)-txx(iz+3,ix-2)) -c4*(txx(iz-3,ix+4)-txx(iz+4,ix-3)))
                dztxx(iz,ix) = (c1*(txx(iz+1,ix+1)-txx(iz  ,ix  )) - c2*(txx(iz+2,ix+2)-txx(iz-1,ix-1)) + c3*(txx(iz+3,ix+3)-txx(iz-2,ix-2)) -c4*(txx(iz+4,ix+4)-txx(iz-3,ix-3)))
                
                dxtzz(iz,ix) = (c1*(tzz(iz  ,ix+1)-tzz(iz+1,ix  )) - c2*(tzz(iz-1,ix+2)-tzz(iz+2,ix-1)) + c3*(tzz(iz-2,ix+3)-tzz(iz+3,ix-2)) -c4*(tzz(iz-3,ix+4)-tzz(iz+4,ix-3)))
                dztzz(iz,ix) = (c1*(tzz(iz+1,ix+1)-tzz(iz  ,ix  )) - c2*(tzz(iz+2,ix+2)-tzz(iz-1,ix-1)) + c3*(tzz(iz+3,ix+3)-tzz(iz-2,ix-2)) -c4*(tzz(iz+4,ix+4)-tzz(iz-3,ix-3)))
                
                dxtxz(iz,ix) = (c1*(txz(iz  ,ix+1)-txz(iz+1,ix  )) - c2*(txz(iz-1,ix+2)-txz(iz+2,ix-1)) + c3*(txz(iz-2,ix+3)-txz(iz+3,ix-2)) -c4*(txz(iz-3,ix+4)-txz(iz+4,ix-3)))
                dztxz(iz,ix) = (c1*(txz(iz+1,ix+1)-txz(iz  ,ix  )) - c2*(txz(iz+2,ix+2)-txz(iz-1,ix-1)) + c3*(txz(iz+3,ix+3)-txz(iz-2,ix-2)) -c4*(txz(iz+4,ix+4)-txz(iz-3,ix-3)))
                

            enddo
        enddo

		 dxtxx = (dztxx+dxtxx) / dx
		 dztzz = (dztzz-dxtzz) / dz
		 dxtxz = (dztxz+dxtxz) / dx
		 dztxz = (2.0*dztxz-dx*dxtxz) / dz

        
! update the psi_function of stress in the CPML
		psidxtxx = bx_half*psidxtxx + ax_half*dxtxx
		psidztzz = bz_half*psidztzz + az_half*dztzz
		psidxtxz = bx_half*psidxtxz + ax_half*dxtxz
		psidztxz = bz_half*psidztxz + az_half*dztxz

		
! take dxtxx as dxtxx+psidxtxx
       dxtxx = dxtxx/chix_half + psidxtxx
        dztzz = dztzz/chiz_half + psidztzz
       dxtxz = dxtxz/chix_half + psidxtxz
        dztxz = dztxz/chiz_half + psidztxz

		   
! update the x and z component of solid and fluid velocity 
    	vx = vx + dt/rhoxz*((dxtxx+dztxz))
    	vz = vz + dt/rhoxz*((dxtxz+dztzz))
 	


!   compute the stress and pressure of nonstiff equations
!   compute the spatial difference of velocity
        dxvx = 0.0
        dzvx = 0.0
        dzvz = 0.0
        dxvz = 0.0
		do ix = 5,nx-3
			do iz = 5,nz-3
			    dzvx(iz,ix) = (c1*(vx(iz  ,ix  )-vx(iz-1,ix-1)) - c2*(vx(iz+1,ix+1)-vx(iz-2,ix-2)) +c3*(vx(iz+2,ix+2)-vx(iz-3,ix-3)) -c4*(vx(iz+3,ix+3)-vx(iz-4,ix-4)))
				dxvx(iz,ix) = (c1*(vx(iz-1,ix  )-vx(iz  ,ix-1)) - c2*(vx(iz-2,ix+1)-vx(iz+1,ix-2)) +c3*(vx(iz-3,ix+2)-vx(iz+2,ix-3)) -c4*(vx(iz-4,ix+3)-vx(iz+3,ix-4)))
				
				dzvz(iz,ix) = (c1*(vz(iz  ,ix  )-vz(iz-1,ix-1)) - c2*(vz(iz+1,ix+1)-vz(iz-2,ix-2)) +c3*(vz(iz+2,ix+2)-vz(iz-3,ix-3)) -c4*(vz(iz+3,ix+3)-vz(iz-4,ix-4)))
				dxvz(iz,ix) = (c1*(vz(iz-1,ix  )-vz(iz  ,ix-1)) - c2*(vz(iz-2,ix+1)-vz(iz+1,ix-2)) +c3*(vz(iz-3,ix+2)-vz(iz+2,ix-3)) -c4*(vz(iz-4,ix+3)-vz(iz+3,ix-4)))
				
			enddo
		enddo	

   dxvx = (dzvx+dxvx) / dx
   dzvx = (2.0*dzvx-dx*dxvx) / dz
   dxvz = (dzvz+dxvz) / dx
   dzvz = (2.0*dzvz-dx*dxvz) / dz

		
! update the psi_function of the CPML
!	the CPML parameter here is at half grid
		psidxvx = bx*psidxvx + ax*dxvx
		psidzvx = bz*psidzvx + az*dzvx
		psidxvz = bx*psidxvz + ax*dxvz
		psidzvz = bz*psidzvz + az*dzvz
		
		
! take dxvx as dxvx+psidxvx	    
        dxvx = dxvx/chix + psidxvx
        dzvx = dzvx/chiz + psidzvx
        dxvz = dxvz/chix + psidxvz
        dzvz = dzvz/chiz + psidzvz
                
			
		
! update the stress and pressure 
       txx = txx + ((c11)*(dxvx) + c12*(dzvz)+c13*(dxvz+dzvx))*dt
        tzz = tzz + ((c22)*(dzvz) + c12*(dxvx)+c23*(dxvz+dzvx))*dt
        txz = txz + ((c13)*(dxvx) + c23*(dzvz)+c33*(dxvz+dzvx))*dt
 
! store seismograms
		do i = 1,nr
			seisvx(it,i) = vx(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
			seisvz(it,i) = vz(fzr+(i-1)*dzr,fxr+(i-1)*dxr)

			
			seistxx(it,i) = txx(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
			seistxz(it,i) = txz(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
			seistzz(it,i) = tzz(fzr+(i-1)*dzr,fxr+(i-1)*dxr)

		end do
		if(sl.eq.1) then
		txx = txx + wavelet(it)*wav
!		tzz = tzz + wavelet(it)*wav
!		txz = txz + wavelet(it)*(1.0-ew)*wav/2.0
   
        endif
	    if(sl.eq.2) then
!		txx = txx + wavelet(it)*wav
		tzz = tzz+wavelet(it)*wav
!		txz = txz + wavelet(it)*(1.0-ew)*wav/2.0
  
        endif
        if(sl.eq.3) then
!		txx = txx + wavelet(it)*wav
!		tzz = tzz + wavelet(it)*wav
		txz = txz + wavelet(it)*(1.0-ew)*wav/2.0
!        p = p + wavelet(it)*wav
        endif
		
		
! output the snapshot
		if(mod(it,dst) .eq. 0) then
			print *,'Time step # ',it
			print *,'Time: ',sngl((it-1)*dt),' seconds'
			print *,'Max norm velocity vector V (m/s) = ',max(maxval(sqrt(vx**2 + vz**2)),maxval(sqrt(vz**2)))
			print *	
			
			write(fname,"('snap_',i6.6,'_Vx.dat')") it
			maxpar = 1.0
			!maxpar = maxval(abs(vx))
			!maxpar = abs(vz(zs,xs))
!			tmp = vx(np+5:nz-4-np,np+5:nx-4-np)
			open(11,file=fname,access='stream')
			write(11) vx
			close(11)
!			call creat_2D_image(vx,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,0,power_display,cutvect)
            write(fname,"('snap_',i6.6,'_Vz.dat')") it
			!maxpar = maxval(abs(vz))
			!maxpar = 1.0
!			tmp = vz(np+5:nz-4-np,np+5:nx-4-np)
			open(12,file=fname,access='stream')
			write(12) vz
			close(12)


            write(fname,"('snap_',i6.6,'_txx.dat')") it
			!maxpar = maxval(txx)
!			tmp = txx(np+5:nz-4-np,np+5:nx-4-np)
			open(15,file=fname,access='stream')
			write(15) txx
			close(15)
!			call creat_2D_image(txx,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,4,power_display,cutvect)

            write(fname,"('snap_',i6.6,'_tzz.dat')") it
			!maxpar = maxval(tzz)
!			tmp = tzz(np+5:nz-4-np,np+5:nx-4-np)
			open(16,file=fname,access='stream')
			write(16) tzz
			close(16)
!			call creat_2D_image(tzz,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,5,power_display,cutvect)

            
            write(fname,"('snap_',i6.6,'_txz.dat')") it
!			!maxpar = maxval(txz)
!			tmp = txz(np+5:nz-4-np,np+5:nx-4-np)
			open(17,file=fname,access='stream')
			write(17) txz
			close(17)
		end if
		
	end do
	
	call cpu_time(time_end)
	write(*,*) 'Time of operation was ',time_end-time_begin,' seconds'


!   write the seismgram
    open(109,file='seisvz.dat',access='stream')
    write(109) seisvz
    close(109)


	open(111,file='seisvx.dat',access='stream')
    write(111) seisvx
    close(111)


	open(114,file='seistxx.dat',access='stream')
    write(114) seistxx
    close(114)


	open(116,file='seistzz.dat',access='stream')
    write(116) seistzz
    close(116)


    deallocate(rhoxz)	
	deallocate(txx,txz,tzz,vx,vz,seisvx,seisvz,seistxx,seistzz,seistxz)
	deallocate(dxtxx,dztxx,dxtxz,dztxz,dxtzz,dztzz,dxvx,dzvz,dxvz,dzvx)
	deallocate(psidxtxx,psidztzz,psidxtxz,psidztxz,psidzvx,psidxvz,psidzvz)
	deallocate(ax,bx,az,bz,chix,chiz,cdx,cdz,ax_half,bx_half,az_half,bz_half,chix_half,chiz_half,cdx_half,cdz_half)
	deallocate(wavelet,wav)

    
	print *
	print *,'End of the simulation'
	pause
	print *
    
end program    
!*********************************
!*     END
!*********************************   



    
!*********************************
!*     Ricker Wavelet
!*********************************
      subroutine ricker(wv,nt,dt,f0,ts)
	  integer*4 nt
      real*4 wv(nt)
      real*4 aa,dt,f0,ts,pi
      integer ii

!cccccccccccccccccccccccccccccccccccccccc
    pi=4.*atan(1.0)
    do ii=1,nt
        aa=pi*f0*((ii-1.)*dt-ts)
        aa=aa*aa
        wv(ii)=(1.-2.*aa)*exp(-aa)
    enddo    
!c---------------------------------------
      return
      end subroutine

!*********************************
!*     dgauss wavelet
!*********************************
    subroutine dgauss(wv,nt,dt,f0,ts)
	integer*4 nt
    real*4 wv(nt)
    real*4 aa,dt,f0,ts,pi
    integer ii
!cccccccccccccccccccccccccccccccccccccccc
    pi = 4. * atan(1.0)
    do ii=1,nt
        aa = pi*f0*pi*f0
        
        wv(ii) = -2.0*aa*((ii-1)*dt-ts)*exp(-aa*((ii-1)*dt-ts)**2)
    end do
!cccccccccccccccccccccccccccccccccccccccc
    return
    end subroutine



!*********************************
!*     draw the snapshot
!*********************************
	subroutine creat_2D_image(sndat,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np,ifield,power_display,cutvect)
	implicit none
!	real,parameter :: power_display=0.2
!	real,parameter :: cutvect=0.05
	logical,parameter :: white_background=.true.
	integer,parameter :: width_cross=5,thickness_cross=1,size_square=3

		integer nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,np,nr
		real,dimension(nz,nx) :: sndat
		!integer,dimension(nr) :: xr,zr
		integer :: ix,iz,irec
		character(len=100) :: fname,system_command
		integer :: R,G,B
		real :: normalized,maxamp
		integer :: ifield
		real :: power_display,cutvect		

		if(ifield .eq. 0) then
			write(fname,"('image',i6.6,'_Vz.pnm')") it
			write(system_command,"('convert image',i6.6,'_Vz.pnm image',i6.6, &
				'_Vz.tif ; rm image',i6.6,'_Vz.pnm')") it,it,it
		else
			write(fname,"('image',i6.6,'_p.pnm')") it
			write(system_command,"('convert image',i6.6,'_p.pnm image',i6.6, &
				'_p.tif ; rm image',i6.6,'_p.pnm')") it,it,it
		endif

		open(unit=27, file=fname, status='unknown')

		write(27,"('P3')") ! write image in PNM P3 format

		write(27,*) nx,nz ! write image size
		write(27,*) '255' ! maximum value of each pixel color

! compute maximum amplitude
		maxamp = maxval(abs(sndat))

! image starts in upper-left corner in PNM format
		  do iz=nz,1,-1
			do ix=1,nx

! define data as vector component normalized to [-1:1] and rounded to nearest integer
! keeping in mind that amplitude can be negative
				normalized = sndat(iz,ix) / maxamp

! suppress values that are outside [-1:+1] to avoid small edge effects
				if(normalized < -1.0) normalized = -1.0
				if(normalized > 1.0) normalized = 1.0

! draw an orange cross to represent the source
				if((ix >= xs - width_cross .and. ix <= xs + width_cross .and. &
					iz >= zs - thickness_cross .and. iz <= zs + thickness_cross) .or. &
					(ix >= xs - thickness_cross .and. ix <= xs + thickness_cross .and. &
					iz >= zs - width_cross .and. iz <= zs + width_cross)) then
						 R = 255
						 G = 157
						 B = 0

! display two-pixel-thick black frame around the image
				else if(ix <= 2 .or. ix >= nx-1 .or. iz <= 2 .or. iz >= nz-1) then
						R = 0
						G = 0
						B = 0

! display edges of the PML layers
				else if((ix == np) .or. &
					(ix == nx - np) .or. &
					(iz == np) .or. &
					(iz == nz - np)) then
						R = 255
						G = 150
						B = 0

! suppress all the values that are below the threshold
				else if(abs(normalized) <= cutvect) then

! use a black or white background for points that are below the threshold
					if(white_background) then
						R = 255
						G = 255
						B = 255
					else
						R = 0
						G = 0
						B = 0
					endif

! represent regular image points using red if value is positive, blue if negative
				else if(normalized >= 0.0) then
					R = nint(255.0*normalized**power_display)
					G = 0
					B = 0
				else
					R = 0
					G = 0
					B = nint(255.0*abs(normalized)**power_display)
				endif

! draw a green square to represent the receivers
				do irec = 1,nr
					if((ix >= fxr+(irec-1)*dxr - size_square .and. ix <= fxr+(irec-1)*dxr  + size_square .and. &
						iz >= fzr+(irec-1)*dzr - size_square .and. iz <= fzr+(irec-1)*dzr + size_square) .or. &
						(ix >= fxr+(irec-1)*dxr  - size_square .and. ix <= fxr+(irec-1)*dxr  + size_square .and. &
						iz >= fzr+(irec-1)*dzr  - size_square .and. iz <= fzr+(irec-1)*dzr + size_square)) then
! use dark green color
							R = 30
							G = 180
							B = 60
					endif
				enddo

! write color pixel
				write(27,"(i3,' ',i3,' ',i3)") R,G,B

			enddo
		enddo

! close file
		close(27)

! call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
		call system(system_command)
		print *,'convert over!'

		
		
	end subroutine 


