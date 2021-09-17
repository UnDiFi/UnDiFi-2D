    subroutine massaction(rhoN,kf,kb,S,is)

!   Compute source term for each chemical species

    implicit none
    include 'plasma.h'
    include 'chem.h'
    include 'paramchem.h'
    include 'commonv.inc'
    include 'commonchem.inc'

    integer i,j,is
    double precision kf(NRMAX),kb(NRMAX),rhoN(NSP),S,Sn(NRMAX)
    double precision qr,qp

    do i = 1,NREAC
        
        qr = 1.d0
        qp = 1.d0
        
        do j = 1,NSP
            qr = qr * rhoN(j)**scr(j,i)
            qp = qp * rhoN(j)**scp(j,i)
        enddo

        Sn(i) = kf(i) * qr - kb(i) * qp
        Sn(i) = Sn(i) * (scp(is,i) - scr(is,i)) 

    enddo 

!   1/(cm^3*s)    
    S = 0.d0
    do i = 1,NREAC
        S = S + Sn(i)
    enddo    

!   g/(mol*cm^3*s)
    S = S * Ms(is)     

    return
    end
