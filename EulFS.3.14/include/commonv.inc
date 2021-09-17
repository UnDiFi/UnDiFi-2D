        integer Nfit,Ncol,Nbac
        parameter(Nfit=5,Ncol=3,Nbac=3)
        integer scr(NSP,NRMAX),scp(NSP,NRMAX)
        real :: Ms(NSP),HF0(NSP),RGASS(NSP)
        real :: mc(Nfit,2*Ncol+1),mbv(2,Nbac)
 	common/commonv/scr, scp, mc, mbv, Ms, HF0, RGASS


