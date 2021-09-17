



      program dat2paraview
	  
      implicit double precision (a-h,o-z)
	  
      integer npoin,nele,ipoin,iele,i,j,IFAIL,icelnod,n
      character*255 fname,execmd,pwd,fwork
      integer nva
      parameter (nva=999000)             
      LOGICAL eulfs,neo 	  
      double precision dumt,Ma,k

      dimension x(nva),y(nva),rho(nva),p(nva),Ma(nva),S(nva),T(nva)
      dimension h(nva),k(nva),u(nva),v(nva),z(nva,4)
      dimension icelnod(3,nva)



c-------------------------------------------------
c     INPUT: number of points and elements.
c------------------------------------------------ 

 

      EXECMD='mkdir tmp'
      IFAIL = SYSTEM(EXECMD)	
      EXECMD='cp *.1.node tmp/.'
      IFAIL = SYSTEM(EXECMD)




!      EXECMD='cd tmp'
!      IFAIL = SYSTEM(EXECMD)
      EXECMD='cd tmp; ls -t | head -n 1 > ../scratch.log'
      IFAIL = SYSTEM(EXECMD)
      EXECMD=' cd ..'
      IFAIL = SYSTEM(EXECMD)	  
      open(17,file='scratch.log')
        read(17,*) fname
      close(17)
      open (9,file=fname)
      read(9,*) npoin,dum,dum,dum  
      close(9)
      EXECMD='cp *.1.ele tmp/.'
      IFAIL = SYSTEM(EXECMD)
      EXECMD='cd tmp; ls -t | head -n 1 > ../scratch.log; cd ..'
      IFAIL = SYSTEM(EXECMD)
	  
      open(17,file='scratch.log')
        read(17,*) fname
      close(17)
	    
	  
      open (9,file=fname)
      read(9,*) nele,dum,dum,dum	  
      close(9)
      
      EXECMD='cd ..'
      IFAIL = SYSTEM(EXECMD)
      
      EXECMD='rm -r tmp'
      IFAIL = SYSTEM(EXECMD)

c/////////////////////////////////////////////////////////////


       gam=1.4d0
       gam1=0.4d0
       
        eulfs=.false.
        neo=.false.
       INQUIRE( FILE='file003.dat', EXIST=eulfs ) 
       INQUIRE( FILE='vvvv.dat', EXIST=neo ) 	   
       
        IF (eulfs.eqv..true.)  then

          open(9,file='file012.dat')
		  do i = 1,4
          read(9,FMT=*)
          end do	  
 
          do i=1, npoin
           read(9,*) x(i),y(i),z(i,1:4)
		   
          rho(i)=z(i,1)**2		   
          h(i)=z(i,2)/z(i,1)			   
          k(i)=0.5*(z(i,3)**2+z(i,4)**2)/rho(i)
          Ma(i)=(k(i)*2d0/(gam1*(h(i)-k(i))))**0.5	
		  p(i)=gam1/gam*rho(i)*(h(i) - k(i) )
		  S(i)=p(i)/(rho(i)**gam)
		  u(i)=z(i,3)/z(i,1)
		  v(i)=z(i,4)/z(i,1)
		  T(i)=p(i)/rho(i)
		  end do

          do i=1,nele
            read(9,*) icelnod(1:3,i)
         end do			

              close(9)

       open(10,file='paraviewE.dat')

        write(10,*) 'TITLE = "Full grid"'
      write(10,*) 'VARIABLES = X Y rho u  v p H  Ma  s T'
      write(10,*) 'ZONE T="sampletext", DATAPACKING=POINT, NODES='
     &,npoin,', ELEMENTS=',nele,', ZONETYPE=FETRIANGLE'
      write(10,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE
     & DOUBLE DOUBLE DOUBLE )'
	 
	 
          do i=1, npoin
           WRITE(10,*) x(i),y(i),rho(i), u(i), v(i), p(i), H(i),  
     +              		   Ma(i),  s(i), T(i)
           end do	 
	 
          do i=1,nele
            write(10,*) icelnod(1:3,i)
         end do		 
         close(10)

       END IF

      IF (neo.eqv..true.)  then
	  
            open(9,file='vvvv.dat')
		  do i = 1,3
          read(9,FMT=*)
		  write(*,*)
          end do	  

          do i=1, npoin
           read(9,*) x(i),y(i),rho(i), u(i), v(i), p(i), H(i),  
     +              		   Ma(i), s(i), T(i),dumt,dumt
           end do
		   
          do i=1,nele
            read(9,*) icelnod(1:3,i)
         end do			

              close(9)

       open(10,file='paraviewN.dat')

        write(10,*) 'TITLE = "Full grid"'
      write(10,*) 'VARIABLES = X Y rho u  v p H  Ma  s T'
      write(10,*) 'ZONE T="sampletext", DATAPACKING=POINT, NODES='
     &,npoin,', ELEMENTS=',nele,', ZONETYPE=FETRIANGLE'
      write(10,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE
     & DOUBLE DOUBLE DOUBLE )'
	 
	 
          do i=1, npoin
           WRITE(10,*) x(i),y(i),rho(i), u(i), v(i), p(i), H(i),  
     +              		   Ma(i),  s(i), T(i)
           end do	 
	 
          do i=1,nele
            write(10,*) icelnod(1:3,i)
         end do		 
         close(10)


         END IF
		 
        


		 
        end program
	 
