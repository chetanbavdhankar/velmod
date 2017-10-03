 
      PROGRAM VELMOD
      implicit none
      integer i,j,io,il,k,ii, nctot,nptot,iter,iseed, iw
      real epsilon, unifd, likmin
      
c     Latest modified by Chetan Bavdhankar, 1/10/2017
c     Subroutines can be found in files in VELMOD folder

**    Internal variables - common
      include     'parm.inc'
      include     'arrays.inc'
      include     'tf.inc'       ! external TF param

**    Data    
      pi = acos(-1.)
      epsilon=1e-10
      write(*,*)'nput random seed'
      read(*,*)iseed
      call getparams 

      nctot=0
      do io=1,3                 ! cicle over the catalogs
         nctot=nctot+ic(io)
      enddo

      write(*,*)'Total number of catalogues',nctot
      if(nctot.ne.icatnum) then
         write(*,*)'Mismatch in the number of catalogues'
         stop
      endif

c     number of free parameter per calalogue = 3 (tf parameters)

      nptot=0
      do io=1,3                 ! cicle over the catalogs
         if(ic(io).eq.1)then
c            write(*,*)io
            write(*,*)'starting TF par in cat',io,'=',A(io),b(io),stf(io)
            write(*,*)'Range of the radial shell',io,'=',zinf(io),zsup(io)
            if(iconst(io).eq.0)nptot=nptot+1
         endif
      enddo
c      write(*,*)'io',io
      nptot=nptot*3
      if(iconst(io).eq.0)nptot=nptot+1  ! the extra free parameter is sigma_v

c     display the number of radial bins used
      write(*,*)'The number of radial bins is........',nbin_int
      write(*,*)'The initial velocity dispersion is..',sigmav
      write(*,*)'The number of free parameters is.....',nptot

c     reads in the various dataset

      call getdata
c      call getdata(dd,vx,vy,vz)

      write(*,*)'Number of catalogs to be analyzed.',icatnum
      do ii=1,3 ! cicle over catalogs
         if(ic(ii).eq.1)then
            write(4,*)'Catalog n. ',ii,' included'    
            write(4,*)'containing ',ngals(ii),' galaxies'
          else
            write(4,*)'Catalog n.',ii,'not included'    
          endif
      enddo      

c     definyig the array of the free parameters to be used as  
c     starting point for the minimization routine.

      il=0
      do ii=1,3  ! cicle over catalogs
         if(ic(ii).eq.1)then
            if(iconst(ii).eq.0)then
               startp(il+1)=A(ii)
               startp(il+2)=b(ii)
               startp(il+3)=stf(ii)
               il=il+3
             endif
         endif
      enddo
c      write(*,*)'ii=',ii
      if(iconst(ii).eq.0)then
         il=il+1
         startp(il)=sigmav
      endif
      if(il.ne.nptot)then
         write(*,*)'Mismatch in the number of free parameters'
         stop
      endif
      do i=1,nptot
         startp0(i)=startp(i)
         startt(i)=0
      enddo

c     here is the minimization routine
c     on input:
c               startp is the starting point for the search
c               xi is the matrix that defines the initial direction
c               nptot is the number of free parameters
c               npartot fixes the physical dimension for startp and xi
c               ftol is the required fractional tolerance in the search
c               iter is the number of required iteration
c     on output:
c               startp is the best set of free parameters            
c               fret is the correponding likelihood value 
c               xi is the current direction
c               iter is the number of iterations taken
c

c     set tolerance

      ftol=0.03
      likmin=100000
      do k=1,1
c     definyif the matrix  used to determine the  
c     starting direction in the minimization routine.
c     (it's a diagonal matrix whch represents the unit 
c     vectors in the space of free parameters)
         
         do i=1,nptot
            do j=1,nptot
               if(i.eq.j)then
                  xi(i,j)=1.0
               else
                  xi(i,j)=0.0
               end if
            enddo
         enddo
 88      write(*,*)'Minimization run',k
         do j=1,nptot
            startp(j)=startp0(j)+(unifd(iseed)-0.5)*startp0(j)/2.5
            write(*,*)'Parameter ',j,'=',startp(j)
         enddo
         iw=0
         call powell(startp,xi,nptot,npartot,ftol,iter,fret,iw)               
         if(iw.eq.1)goto 88
         if(fret.le.likmin)likmin=fret
         il=0
         do ii=1,3        ! cicle over catalogs
            if(ic(ii).eq.1.and.iconst(ii).eq.0)then
               write(*,*)'Catalogue number',ii
               write(*,*)'Zeropoint for the TF...',startp(il+1)
               write(*,*)'Slope for the TF.......',startp(il+2)
               write(*,*)'dispersion for the TF..',startp(il+3)
               if(fret.le.likmin)then
                  startt(il+1)=startp(il+1)
                  startt(il+2)=startp(il+2)
                  startt(il+3)=startp(il+3)
                endif
               il=il+3
            endif
         enddo
         il=il+1
         if(iconst(il).eq.0)then
            write(*,*)'Velocity dispersion....',startp(il)
            if(fret.le.likmin)startt(il)=startp(il)
         endif
         
         write(*,*)'Value of the likelihood.',fret
         
      enddo
      

      write(*,*)'BEST VALUES'
      write(*,*)'Value of the likelihood.',likmin
      il=0
      do ii=1,3           ! cicle over catalogs
         if(ic(ii).eq.1.and.iconst(ii).eq.0)then
            write(*,*)'Catalogue number',ii
            write(*,*)'Zeropoint for the TF...',startt(il+1)
            write(*,*)'Slope for the TF.......',startt(il+2)
            write(*,*)'dispersion for the TF..',startt(il+3)            
            write(12,44)ii,beta,startt(il+1),startt(il+2),startt(il+3)
     1            ,startt(il+4),likmin
            il=il+3
         endif
      enddo
      il=il+1
      if(iconst(il).eq.0)then
         write(*,*)'Velocity dispersion....',startt(il)
         write(13,*)beta,startt(il),likmin,ngals(1)
      endif
      
 44   format(1x,i2,6(1x,f9.3))
      
      stop
      end


