!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       NAHOON_1D.F90                                                                          
!       Wakelam V. and Hersant F. Jan. 2011 
!
! 	Chemical eovlution is computed as a function of time. 
!       This version of the code can compute the chemistry in more than 
!       one shell of gas but does not include any diffusion terme.
!       The jacobian is explicitly given in this version of the code.
!	Output times are read in an input file (timeres.dat). 
! 	In 0D, the temperature and density are read in the input_parameter.dat 
!	file.
!	In 1D, the temperature and density can be read in an external file but the 
!	user must change the format in the code.                         
!                                                                                                                                                                                 
!       DLSODES integretor is used.    
!                                             
! 	Output files:    
!       - Kout.dat contain the values of the rate coefficients          
!       - plot.dat the abundances as a function of time                 
!       - verif.dat contains some balance verifications and the         
!       rates of formation and destruction for each species at          
!       specific times                                                  
!       - output.dat is the chemical composition at one time (by default
!       in the same format as cond_initial.dat  
!	- photodiss_XX.dat files give the column densities of H2 and (H2 and CO) 
!       with the photorates  
!
!	April 2011 VW  Modifications: 
!		- read new format of KIDA 
!		- including a test on the T range of validity of rate coefficients   
!		- including a test on duplicated reactions
!		- changing the H2 formation rate              
!       August 2013 VW Bug fixing
!       STRONG WARNING: FOR TEMPERATURES OUTSIDE THE T RANGE, RATE COEFFICIENTS ARE 
!       NOT EXTRAPOLATED                              
!              - modification of the way rate coefficients are calculated accoding to T range
!                  in the case of multiple reactions                      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
                                                                        
       PROGRAM NAHOON_1D
                                                                        
       IMPLICIT NONE 
                                                                        
       include 'header.f90' 
                                                                        
       INTEGER ITYPE,IPROD,                                        &
     &              ITEST,ITEST2,ELEMENT,I,J,ISP,JSTEP,L,NPLOT,IH,      &
     &              IGRAIN0,IGRAINN,JSPACE,IE,ICO,IH2,NUMBER,FORMULA                      
       PARAMETER (NPLOT=NS) 
       DOUBLE PRECISION SN,Y,AB,PLOTAB,TPLOT,TD,TIMERES,TYR,            &
     &              TIMESELEC,XK1,YGRAIN,XCO,XH2,A,B,C,RANDOM1,UNC,ZETA,&
     &              GTODN,TAU,NCO,T_CO,NH2_1,T_H2_1,AV,T_AV,NH2_2,T_H2_2,&
     &          SNGRAIN,ABBIS,NH2TOT,PLUS,XELEC,ALPHA,Tmin,Tmax                        
                                                                        
       CHARACTER*10 SPEC,REACTANT 
       DIMENSION SN(NS,nx),Y(NS,nx),AB(NS,nx),PLOTAB(NTIME,NPLOT,nx),&
     &              A(NRTOT),B(NRTOT),C(NRTOT),TD(NX),ABBIS(NS),XK1(NX),&
     &              TAU(NX),NH2TOT(NX),TIMERES(NTIME),SPEC(NS+1),ITEST(6),&
     &              ELEMENT(NELEM,NS+1),REACTANT(3,NRTOT),ITYPE(NRTOT),   &
     &              NCO(52),NH2_2(105),Tmin(NRTOT),Tmax(NRTOT),       &
     &              T_H2_2(105),T_CO(52),NH2_1(43),T_H2_1(43),AV(43),       &
     &              T_AV(43),PLUS(NX),XELEC(NX),XCO(NX),XH2(NX),ALPHA(NRTOT), &
     &              NUMBER(NRTOT),FORMULA(NRTOT)                                            
       CHARACTER AA*78,RANDMODE 
                                                                        
       INTEGER REACT(NRTOT,7) 
       DOUBLE PRECISION K(NRTOT,NX),NHTOT(NX),RADIUS(2,NX) 
       COMMON/EQUA/K,NHTOT,RADIUS,REACT 
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        Variable declaration FOR THE DLSODES SUBROUTINE                
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       EXTERNAL F, JAC, JACFH
       INTEGER NEQ,ITOL,ITASK,ISTATE,IOPT,LRW,IWORK,LIW,JEX,MF 
       DOUBLE PRECISION T, TOUT, RTOL, ATOL, RWORK 
       PARAMETER (NEQ=NS*NX,LIW=31+NEQ+2000*NEQ,LRW=20+9*NEQ+2000*NEQ) 
       DIMENSION IWORK(LIW), RWORK(LRW)                                 &
     &       , ATOL(NS,NX)                                              
                                                                        
       DATA ITOL, ITASK, ISTATE, IOPT, MF, RTOL                         &
     &              /2, 1, 1, 1, 21, 1.D-6/                            
! Specific variables for the jacobian geometry
       integer, dimension(NS*NX+1) :: IA, IAN
       integer, dimension(liw) :: JA, JAN
       integer :: NNZ
       integer, dimension(NEQ*NEQ+20) :: IWK
       equivalence(iwk(1),rwork(21))

                                                                        
       IWORK(6)=1000
       RWORK(6)=3.154E+15 
       IWORK(5)=5 
       IWORK(7)=10 
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Files declaration                                                                
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                        
       OPEN (UNIT=1,FILE='kida.uva.2014',STATUS='OLD') 
       OPEN (UNIT=2,FILE='cond_initial_kida.uva.2014.dat',STATUS='OLD') 
       OPEN (UNIT=3,FILE='input_parameter.dat',STATUS='OLD') 
       OPEN (UNIT=4,FILE='timeres.dat',STATUS='OLD') 
       OPEN (UNIT=7,FILE='Self_Shielding_data',STATUS='OLD') 
                                                                        
       OPEN (UNIT=9,FILE='Kout.dat',STATUS='UNKNOWN') 
       OPEN (UNIT=10,FILE='plot.dat',STATUS='UNKNOWN') 
       OPEN (UNIT=11,FILE='verif.dat',STATUS='UNKNOWN') 
       OPEN (UNIT=12,FILE='output.dat',STATUS='UNKNOWN')
       OPEN (UNIT=13,FILE='photodiss_H2.dat',STATUS='UNKNOWN')
       OPEN (UNIT=14,FILE='photodiss_CO.dat',STATUS='UNKNOWN')
      
                                                                        
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
!       read data files                                                 
!                                                                       
!       AB is the abundance species with respect to total               
!       SN is in number density (cm-3)                                  
                                                                        
!       First put some numbers to zero                                  
       DO I=1,NS 
              DO JSPACE=1,NX 
                     SN(I,JSPACE)=0.D0 
                     Y(I,JSPACE)=0.D0 
                     AB(I,JSPACE)=0.D0 
                     PLUS(JSPACE)=0.D0
              ENDDO 
              DO ISP=1,NELEM 
                     ELEMENT(ISP,I)=0 
              ENDDO 
       ENDDO 
       DO ISP=1,NELEM 
              ELEMENT(ISP,NS+1)=0 
       ENDDO 
                                                                       
       DO I=1,NRTOT 
              DO JSPACE=1,NX 
                  K(I,JSPACE) = 0.0D0 
              enddo 
       ENDDO 
                                                                        
       DO I=1,NX 
              TD(I)=0.D0 
              NHTOT(I)=0.D0 
              TAU(I)=0.D0 
       ENDDO 
                                                                        
       T=0.D0 
       JSTEP=0 
                                                                        
!       read output times                                               
       READ(4,*) TIMERES 
                                                                        
!       read initial conditions from cond_initial.dat (abundance /H)                          
!       WHEN THE COND_INITIAL.DAT IS 0D    
!       if you want to read 1D initial conditions, the format has to 
!       be changed here                        
       DO I=1,1
              READ (2,*) AA 
       ENDDO 
       DO I=1,NS 
              READ (2,1) SPEC(I),(ELEMENT(ISP,I),ISP=1,NELEM),          &
     &              ABBIS(I)                                            
               DO JSPACE=1,NX 
                     AB(I,JSPACE)=ABBIS(I) 
                     IF (ELEMENT(1,I).EQ.1) PLUS(JSPACE)=PLUS(JSPACE)+AB(I,JSPACE) 
              ENDDO 
              IF (SPEC(I).EQ.'H         ') IH=I 
              IF (SPEC(I).EQ.'GRAIN0    ') IGRAIN0=I 
              IF (SPEC(I).EQ.'GRAIN-    ') IGRAINN=I 
	      IF (SPEC(I).EQ.'e-        ') IE=I                                                             
              IF (SPEC(I).EQ.'CO        ') ICO=I
              IF (SPEC(I).EQ.'H2        ') IH2=I
       ENDDO 
    1   FORMAT(5X,A10,14(I3),2X,D14.8) 

                                                                        
!       for the additional blank species                                
        SPEC(NS+1)='          ' 

!       read the parameters from the input_parameter.dat file
! 	if the number of physical points are larger than one, T and n are 
!	read in another file                                            

	if (nx.gt.1) then 
       DO I=1,9 
              READ (3,14) AA 
       ENDDO 
       READ(3,2) RANDMODE 
       READ (3,14) AA 
       READ (3,14) AA 
       READ (3,14) AA 
       READ(3,3) TIMESELEC 
          endif
	 
! 	if the number of physical points is one, T and n are 
!	read in the input_parameter.dat file                                            
	 
	 if (nx.eq.1) then 
	 DO I=1,9 
              READ (3,14) AA 
       ENDDO 
       READ(3,2) RANDMODE 
       READ(3,3) NHTOT(1) 
       READ(3,3) TD(1) 
       READ(3,3) TAU(1)                                                  
       READ(3,3) TIMESELEC 
	endif
 14    FORMAT(A78) 
 2   FORMAT(A1) 
 3   FORMAT(D9.3) 


!-----------------------------------------------------------------                                                                        
!        READ THE TEMPERATURE, DENSITY AND VISUAL EXTINCTION IN AN      
!       EXTERNAL FILE, current format of ratran code output             
!       THIS HAS TO BE CHANGED ACCORDING TO YOUR CASE !!!!!             
       
!       if (nx.gt.1) then                                                                 
!       DO I=1,7 
!              READ (8,14) AA 
!       ENDDO 
!       DO I=1,NX               
!              READ (8,*) RADIUS(1,I),RADIUS(2,I),NHTOT(I),   &
!     &        TD(I),TAU(I)                                          
!!              NHTOT(I)=2.D0*NH2TOT(I) 
!       ENDDO 
!  100     FORMAT(3X,3E13.6,13X,F7.3)                                                                    
!	endif
!-----------------------------------------------------------------                                                                        
	                                                          
!       compute the species densities from the initial abundances         
!       SNGRAIN is the total abundance of grains                        
       DO I=1,NS 
              DO JSPACE=1,NX 
!                     IF (AB(I,JSPACE).EQ.0.D0) AB(I,JSPACE)=1.D-50     
                     SN(I,JSPACE)=AB(I,JSPACE)*NHTOT(JSPACE) 
                     PLOTAB(1,I,JSPACE)=AB(I,JSPACE) 
              ENDDO 
        ENDDO 
	
!      COMPUTE THE TOTAL ABUNDANCE OF GRAINS	
       SNGRAIN=SN(IGRAIN0,1)+SN(IGRAINN,1) 

!        COMPUTE THE INITIAL ABUNDANCE OF ELECTRONS, THIS IS TO BE SURE 
!	THAT WE HAVE NEUTRALITY
        
	DO JSPACE=1,NX
	   XELEC(JSPACE)=AB(IE,JSPACE)
	   AB(IE,JSPACE)=PLUS(JSPACE)
	ENDDO
	   
        DO I=1,NS
	   DO JSPACE=1,NX 
	      IF ((ELEMENT(1,I).EQ.-1).AND.(I.NE.IE)) AB(IE,JSPACE)=AB(IE,JSPACE)-AB(I,JSPACE)
           ENDDO
        ENDDO 
        IF (AB(IE,1).NE.XELEC(1)) then 
	
	WRITE(*,*) 'THE INITIAL ABUNDANCE OF ELECTRONS HAVE BEEN MODIFIED TO HAVE NEUTRALITY'
        write(*,*) 'previous electron abundance ', XELEC
        write(*,*) 'new electron abundance ', (ab(ie,jspace), jspace=1,nx)
	  endif 
	    
!        WE SET THE EXTERNAL CONDITIONS FOR THE H2 AND CO SELF-SHIELDING. POINT 1 OF 
!	JSPACE HAS TO BE THE EXTERNAL POINT WHERE AV IS THE SMALLER ONE.       
      XCO(1)=XCO_0
      XH2(1)=XH2_0                                               

	    
	                                                           
!       read the chemical database                                      
       CALL READ_REACT(ELEMENT,A,B,C,REACTANT,         &
     &       SPEC,ITYPE,ZETA,GTODN,TAU,NCO,T_CO,NH2_1,      &
     &       T_H2_1,AV,T_AV,NH2_2,T_H2_2,Tmin,Tmax,NUMBER,FORMULA)                               
                                                                        
                                                                        
!
!       compute the rate coefficients                                   
!                                                                                                                                              
       CALL RATE_COEFF(A,B,C,GTODN,TAU,ZETA,TD,XH2,XCO,                 &
     &              NCO,T_CO,NH2_1,T_H2_1,AV,T_AV,NH2_2,T_H2_2, &
     &                   JSTEP,ITYPE,RANDMODE,REACTANT,SNGRAIN,Tmin,Tmax,NUMBER,FORMULA)     
                                                                        

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
!       compte the chemical evolution                                   
!                                                                       
                                                                        
!       definition of the times at which the rates of formation/destruct
!       are writen in verif.dat                                         
       ITEST(1) = 2 
       ITEST(2) = 11 
       ITEST(3) = 34 
       ITEST(4) = 56 
       ITEST(5) = 79 
       ITEST(6) = 101 
       ITEST2=1 
                                                                               
                                                                        
       DO ISP=1,NS 
              DO JSPACE=1,NX 
                     Y(ISP,JSPACE) = AB(ISP,JSPACE)                                                   
              ENDDO 
       ENDDO 
                                                                        
       TOUT=0.D0 
                                                                                                                                                
       ISTATE=1 
                                                                        
!       starting the loop                                               
                                                                        
       WRITE (*,*) 'It is computing.... wait !' 
                                                                        
       DO JSTEP=2,NTIME 
                                                                        
!       TOUT is the output time of the dlsode subroutine and is set to b
!       time in sec from the table timeres                              
!       TYR is the time in year writen in the output file               
       TOUT=TIMERES(JSTEP)*3.154D7 
       TYR=TIMERES(JSTEP) 
                                                                        
!       we transfer the values of the species densities into the        
!       variable Y used by dlsode                                       
       DO ISP=1,NS 
              DO JSPACE=1,NX 
                     Y(ISP,JSPACE) = AB(ISP,JSPACE) 
!       defining atol                                                   
              ATOL(ISP,jspace)=1.D-30                                  
                     if ((Y(ISP,JSPACE)*1.D-12).ge.1.D-24)              &
     &                     ATOL(ISP,jspace)=Y(ISP,JSPACE)*1.D-12        
              ENDDO 
       ENDDO 
                                                                               
                                                                        
! Fancy "do while" loop
                  
       do while (T.lt.TOUT)
                                                      
!------------------------------------------------------------------
! Determine the Jacobian geometry

       IWORK(31:liw)=0

       call computeIAJA(IA,JA,liw,SPEC,Y)
       NNZ=IA(NEQ+1)-1
       iwork(30+1:30+NEQ+1)=IA(1:NEQ+1)
       iwork(31+NEQ+1:31+NEQ+NNZ)=JA(1:NNZ)

!------------------------------------------------------------------

       ISTATE=1                                                        
!       we call the subroutine to solve the differential equations      
       CALL DLSODES (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,       &
     &                ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JACFH, MF)    


      enddo ! fancy while loop


       print *,'TIME = ',T/3.154d7, 'YR'
                                                                        
!       put the densities as abundances compared to total H             
       DO JSPACE=1,NX 
              DO ISP=1,NS 
              IF (SN(ISP,JSPACE).LT.0.D0) SN(ISP,JSPACE)=0.D0 
              SN(ISP,JSPACE) = Y(ISP,JSPACE)*NHTOT(JSPACE) 
              AB(ISP,JSPACE) = Y(ISP,JSPACE) 
              ENDDO 
                                                                        
!       writing the abundances in the PLOTAB for the plot.dat file      
              DO ISP=1,NS 
                         PLOTAB(JSTEP,ISP,JSPACE)=AB(ISP,JSPACE) 
              ENDDO 
       ENDDO 
                                                                        
                                                                        
!       write the verif.dat file, the output.dat file                   
       if (itest2.ne.7) then
          IF (JSTEP.EQ.ITEST(ITEST2)) THEN 
              WRITE(*,*) 'T(YR)=',TYR 
              WRITE(11,*) 'T(YR)=',TYR 
                                                                        
              CALL CHECKING(ELEMENT,SN,IGRAIN0,IGRAINN) 
              CALL COMP_RATES(SN,SPEC,NUMBER) 
              ITEST2=ITEST2+1 
          ENDIF 
       endif 
                                                                        
       IF (TYR.EQ.TIMESELEC) THEN 
       WRITE(12,5) TIMESELEC,' +-  H  He C  O  S  G    ' 
    5   FORMAT(D8.2,A70) 
       DO JSPACE=1,NX 
       write(12,*) 'JSPACE = ',JSPACE 
              DO I=1,NS 
              WRITE(12,6) I,SPEC(I),(ELEMENT(ISP,I),ISP=1,NELEM),       &
     &              AB(I,JSPACE)                                        
              ENDDO 
       ENDDO 
       ENDIF 
    6   FORMAT (I4,1X,A10,14(I3),2X,D14.8) 
              
	                                                                
!	HERE WE DEFINE THE ABUNDANCES OF CO AND H2 THAT WE HAVE IN FRONT OF 
!	EACH POINT       
        IF (NX.GT.1) THEN
	DO JSPACE=2,NX 
       	      XCO(JSPACE)=XCO(JSPACE-1)+Y(ICO,JSPACE-1)
      	      XH2(JSPACE)=XH2(JSPACE-1)+Y(IH2,JSPACE-1)
	ENDDO      
        ENDIF
	                                                      
!       re-compute the rate coefficients.                               
!       This is useful if one wants the temperature to evolve and       
!       to compute the H2 and CO self-shielding which depend on the H2 a
        CALL RATE_COEFF(A,B,C,GTODN,TAU,ZETA,TD,XH2,XCO,                 &
     &              NCO,T_CO,NH2_1,T_H2_1,AV,T_AV,NH2_2,T_H2_2, &
     &                   JSTEP,ITYPE,RANDMODE,REACTANT,SNGRAIN,Tmin,Tmax,NUMBER,FORMULA)     
                                                                        
       ENDDO 
                                                                        
!       end of big loop                                                 
                                                                        
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                 
!       write the plot.dat file for all times                           
                                                                        
       DO JSPACE=1,NX 
              WRITE (10,*) 'Jspace = ',JSPACE 
              WRITE (10,7) 'Gas temperature (K)', TD(JSPACE) 
              WRITE (10,7) 'Total H density (cm-3)',NHTOT(JSPACE) 
    7          FORMAT(A23,E10.2) 
              WRITE(10,8) TIMERES 
    8          FORMAT(10X,124(2X,D14.8)) 
              DO ISP=1,NS 
              WRITE(10,9) SPEC(ISP),(PLOTAB(I,ISP,JSPACE),I=1,NTIME) 
              ENDDO 
    9          FORMAT(A10,124(2X,D14.8)) 
       ENDDO 
                                                                        
       STOP 
      END                                           
                                                                        
!       end of main program                                             
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc                  
                                                                        
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                   
!       this subroutine reads the reactions in the database             
                                                                        
       SUBROUTINE READ_REACT(ELEMENT,A,B,C,REACTANT,   &
     &              SPEC,ITYPE,ZETA,GTODN,TAU,NCO,T_CO,NH2_1,&
     &              T_H2_1,AV,T_AV,NH2_2,T_H2_2,Tmin,Tmax,NUMBER,FORMULA)                        
                                                                        
       IMPLICIT NONE 
                                                                        
       include 'header.f90' 
                                                                        
       DOUBLE PRECISION A,B,C,ZETA,DTOGM,RD,RHOD,TAU,       &
     &              NCO,T_CO,NH2_1,T_H2_1,AV,T_AV,NH2_2,T_H2_2,GTODN,   &
     &          SUM1,SUM2,AMH,BOLTZ,PI,MOLE,Tmin, Tmax                           
       INTEGER ITYPE,ELEMENT,IR1,IR2,IR3, Tmini,Tmaxi,     &
     &              DES1,DES2,IPROD1,IPROD2,IPROD3,IPROD4,L,I,J, NUMBER, FORMULA         
                                                                        
       CHARACTER REACTANT*10,PRODUIT*10,SPEC*10,AA*78
       DIMENSION ITYPE(NRTOT),A(NRTOT),B(NRTOT),C(NRTOT), TAU(NX),NCO(52),    &
     &           T_CO(52),NH2_1(43),T_H2_1(43),AV(43),T_AV(43),        &
     &           NH2_2(105),T_H2_2(105),ELEMENT(NELEM,NS+1),         &
     &           REACTANT(3,NRTOT),PRODUIT(4,NRTOT),SPEC(NS+1),Tmin(NRTOT),&
     &           Tmax(NRTOT),NUMBER(NRTOT), FORMULA(NRTOT)       
                                                                        
       INTEGER REACT(NRTOT,7) 
       DOUBLE PRECISION K(NRTOT,NX),NHTOT(NX),RADIUS(2,NX) 
       COMMON/EQUA/K,NHTOT,RADIUS,REACT 
                                                                        
       DATA AMH,BOLTZ,PI,MOLE                                           &
     &     /1.66043D-24, 1.38054D-16, 3.14159265, 6.023D23/             
                                                                        
                                                                        
!       reading other parameters                                        
       READ(3,10) ZETA 
       READ(3,10) DTOGM 
       READ(3,10) RD 
       READ(3,10) RHOD 
                                                                        
   10    FORMAT (D9.3) 
                                                                        
!       reading the CO and H2 self-shielding parameters                 
       READ(7,14) AA                                                   
       DO J=1,43                                                       
              READ(7,12) NCO(J),T_CO(J),NH2_1(J),T_H2_1(J),AV(J),T_AV(J)
       ENDDO                                                           
       DO J=44,52                                                      
              READ(7,13) NCO(J),T_CO(J)                                
       ENDDO                                                           
                                                                       
       READ(7,15) AA                                                   
14       FORMAT(A52)                                                   
15       FORMAT(A17)                                                   
                                                                       
       DO J=1,105                                                      
             READ(7,13) NH2_2(J),T_H2_2(J)                            
       ENDDO                                                           
12       FORMAT (6(E9.3,3X))                                           
13       FORMAT (E9.3,3X,E9.3)                                         
                                                                        
!       reading the chemical database                                   
       DO J=1,NRTOT 
              READ (1,11) REACTANT(1,J),REACTANT(2,J),REACTANT(3,J),    &
     &         PRODUIT(1,J),PRODUIT(2,J),PRODUIT(3,J),PRODUIT(4,J),     &
     &          A(J),B(J),C(J), ITYPE(J),   Tmini,Tmaxi, FORMULA(J)  , NUMBER(J)                      
 
 
        Tmin(J)=real(Tmini)
        Tmax(J)=real(Tmaxi)
                                                             
                                                                   
!       replace the species names by blanks for non chemical species                                                                        
      		if ((REACTANT(1,J).eq.'CR        ').or.(REACTANT(1,J).eq.'CRP       ').or.&
     &		(REACTANT(1,J).eq.'Photon    ')) reactant(1,J)='          '
     		if ((REACTANT(2,J).eq.'CR        ').or.(REACTANT(2,J).eq.'CRP       ').or.&
     &		(REACTANT(2,J).eq.'Photon    ')) reactant(2,J)='          '
     		if ((REACTANT(3,J).eq.'CR        ').or.(REACTANT(3,J).eq.'CRP       ').or.&
     &		(REACTANT(3,J).eq.'Photon    ')) reactant(3,J)='          '
      		if ((Produit(1,J).eq.'CR        ').or.(PRODUIT(1,J).eq.'CRP       ').or.&
     &		(PRODUIT(1,J).eq.'Photon    ')) PRODUIT(1,J)='          '
      		if ((Produit(2,J).eq.'CR        ').or.(PRODUIT(2,J).eq.'CRP       ').or.&
     &		(PRODUIT(2,J).eq.'Photon    ')) PRODUIT(2,J)='          '
      		if ((Produit(3,J).eq.'CR        ').or.(PRODUIT(3,J).eq.'CRP       ').or.&
     &		(PRODUIT(3,J).eq.'Photon    ')) PRODUIT(3,J)='          '
      		if ((Produit(4,J).eq.'CR        ').or.(PRODUIT(4,J).eq.'CRP       ').or.&
     &		(PRODUIT(4,J).eq.'Photon    ')) PRODUIT(4,J)='          '
                                                                        
       ENDDO 
     
  11    FORMAT (3(a10,1x),1x,4(a10,1x),11x,3(e11.3),24x,i2,2(i7),i3,I6)
                                                                                                                                         
                                                                        
!       GTODN is the ratio between the density of gas and the density of
       GTODN=(4.d0*PI*RHOD*RD*RD*RD)/(3.d0*DTOGM*AMH) 
       write(*,*) GTODN                                                     
							               
!       construct the table REACT which is used to construct the differe
        DO I=1,NRTOT 
              DO J=1,NS+1 
              DO L=1,3 
              IF (REACTANT(L,I).EQ.SPEC(J)) REACT(I,L)=J 
              ENDDO 
              DO L=1,4 
              IF (PRODUIT(L,I).EQ.SPEC(J)) REACT(I,L+3)=J 
              ENDDO 
              ENDDO 
       ENDDO 
                                                                        
!       check the balance of the reactions (elements and charges)       
      DO I=1,NRTOT 
              IR1=REACT(I,1) 
              IR2=REACT(I,2) 
              IR3=REACT(I,3) 
              IPROD1=REACT(I,4) 
              IPROD2=REACT(I,5) 
              IPROD3=REACT(I,6) 
              IPROD4=REACT(I,7) 
              DO J=1,NELEM 
              SUM1=ELEMENT(J,IR1)+ELEMENT(J,IR2)+ELEMENT(J,IR3) 
              SUM2=ELEMENT(J,IPROD1)+ELEMENT(J,IPROD2)+                 &
     &                     ELEMENT(J,IPROD3)+ELEMENT(J,IPROD4)          
                   IF (SUM1.NE.SUM2)                                    &
     &              WRITE(*,*) 'BALANCE PROBLEM AT THE REACTION',I      
              ENDDO 
	      
       ENDDO   
       
!       check for the reactions that are present more than once in the network



       
                    
                                                                        
       RETURN 
      END                                           
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
                                                                        
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
!       this subroutine computes the rate coefficients of the reactions 
!                                                                       
!       the formula to compute the rates depends on the type of reaction                                  
                                                                        
!       if you add a new kind of reaction, you have to add here the way 
!       to compute the rate coefficient                                 
                                                                        
                                                                        
       SUBROUTINE RATE_COEFF(A,B,C,GTODN,TAU,ZETA,TD,XH2,XCO,           &
     &              NCO,T_CO,NH2_1,T_H2_1,AV,T_AV,NH2_2,T_H2_2, &
     &                   JSTEP,ITYPE,RANDMODE,REACTANT,SNGRAIN,Tmin,Tmax,NUMBER,FORMULA)     
                                                                        
       IMPLICIT NONE 
                                                                        
       INCLUDE 'header.f90' 
                                                                        
       DOUBLE PRECISION MOLE,GTODN,                                     &
     &         PI,AMH,BOLTZ,TT,TD,ZETA,A,B,C,TAU,ALPHA,          &
     &         NCO,T_CO,NH2_1,T_H2_1,AV,T_AV,NH2_2,T_H2_2,         &
     &         SNGRAIN,NCOLLH2,XH2,XCO,TETA,NCOLLCO,TETA1,TETA2,TETA3
        DOUBLE PRECISION U,V,CHD,MG,RD,RHOD,EV,MU0,DS,TCOLL,SIGD ,Tmin,Tmax
        real(kind=8), dimension(10) :: distmin, distmax
                                                                
       INTEGER ITYPE,L,M,N,W,JSTEP,I,J,ISP,JSPACE,NUMBER,INDICE ,FORMULA
       CHARACTER REACTANT*10,RANDMODE,AA*78 
       DIMENSION ITYPE(NRTOT),A(NRTOT),B(NRTOT),C(NRTOT), TAU(NX),      &
     &              REACTANT(3,NRTOT),ALPHA(NRTOT),Tmin(NRTOT),Tmax(NRTOT), &
     &              NCO(52),T_CO(52),NH2_1(43),T_H2_1(43),AV(43),T_AV(43),&
     &              NH2_2(105),T_H2_2(105),TD(NX),TT(NX),XCO(NX),XH2(NX),  &
     &              INDICE(10),NUMBER(NRTOT),FORMULA(NRTOT)                
                                                                        
       INTEGER REACT(NRTOT,7) 
       DOUBLE PRECISION K(NRTOT,NX),NHTOT(NX),RADIUS(2,NX) 
       COMMON/EQUA/K,NHTOT,RADIUS,REACT 
                                                                        
       DATA AMH,BOLTZ,PI,MOLE                                        &
     &     /1.66043D-24, 1.38054D-16, 3.14159265, 6.023D23/   
                                                                        

!      READ THE ALPHAS IF RANDMODE EQ 'Y'
       
       IF (RANDMODE.EQ.'Y') THEN                                
          OPEN (UNIT=5,FILE='alpha',STATUS='OLD') 
	  DO J=1,NRTOT
	  	READ(5,*) ALPHA(J)
!		REPLACE THE A VALUES BY THE ONES COMPUTED BY THE 
!    		MONTE CARLOT SIMULATION		
		A(J)=ALPHA(J)
	  ENDDO
	  close(unit=5)
       ENDIF	  			       

	                                                                    
       DO JSPACE=1,NX 
              TT(JSPACE)=TD(JSPACE)/300.D0 
       ENDDO 
                                                                        
       DO JSPACE=1,NX 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	reactions for kida networks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	W=1
	N=0
	INDICE(:)=0

!       start computing the rate coefficients depending on the formula
       DO J=1,NRTOT                                                             


!	REACTIONS WITH GRAINS

         IF (FORMULA(J).EQ.0) THEN 
              K(J,jspace)=A(J)*(TT(jspace)**B(J))
         ENDIF 
                  	
          IF  	(FORMULA(J).EQ.10) K(J,jspace)=A(J)*1.186D7*exp(225.D0/TD(jspace))**(-1)*GTODN/NHTOT(jspace)
          IF  	(FORMULA(J).EQ.11) K(J,jspace)=A(J)*(TT(jspace)**B(J))*NHTOT(jspace)/GTODN
                                             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	REACTIONS WITH COSMIC-RAYS

	   IF (FORMULA(J).EQ.1) K(J,jspace)=A(J)*ZETA
	   
	   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	PHOTOREACTIONS

	   IF (FORMULA(J).EQ.2) THEN
	       K(J,jspace)=A(J)*EXP(-C(J)*TAU(jspace))
		                                                                        
!             computation of the H2 and CO self-shielding by linear extrapola

              IF (REACTANT(1,J).EQ.'H2        ') THEN                   
                  NCOLLH2=(TAU(jspace)/5.34D-22)*XH2(JSPACE)   
                  TETA=1                                            
                  DO L=1,104                                        
                     IF ((NH2_2(L).LE.NCOLLH2).AND.(NH2_2(L+1).GE.NCOLLH2))         &                
     &                            TETA=T_H2_2(L)+(NCOLLH2-NH2_2(L))*                &
     &                            (T_H2_2(L+1)-T_H2_2(L))/(NH2_2(L+1)-NH2_2(L))
                  ENDDO                                             
 		                                                                                
                                                                       
!                     IF (NCOLLH2.GT.NH2_2(105)) K(J,JSPACE)=A(J)*EXP(-C(J)*TAU(jspace))
!
!     VW Sept 08, we change the way the photodissociation is computed when the H2
!     column density exceeds the largest value in Self-Shielding_data file, we use the smallest value of teta. 
!     If H2 column density is smaller 
!     than the smallest value in the file, we use a photorate of 2.54e-11 (teta=1).

                  IF (NCOLLH2.GT.NH2_2(105)) TETA=T_H2_2(105)

		  K(J,JSPACE)=2.54D-11*TETA

			write(13,*)  NCOLLH2, K(J,JSPACE)

              ENDIF                                               
	           
                  IF (REACTANT(1,J).EQ.'CO        ') THEN           
                     NCOLLH2=(TAU(jspace)/5.34D-22)*XH2(JSPACE)                     
                     NCOLLCO=(TAU(jspace)/5.34D-22)*XCO(JSPACE)                        
                     TETA1=1                                           
                     TETA2=1                                           
                     TETA3=1                                           
                     DO L=1,51                                         
                            IF ((NCO(L).LE.NCOLLCO).AND.(NCO(L+1).GE.NCOLLH2))             &
     &                            TETA2=T_CO(L)+(NCOLLCO-NCO(L))*(T_CO(L+1)-T_CO(L))        &
     &                            /(NCO(L+1)-NCO(L))                                
                     ENDDO                                             
                     DO L=1,42                                         
                            IF ((NH2_1(L).LE.NCOLLH2).AND.(NH2_1(L+1).GE.NCOLLH2))             &
     &                            TETA1=T_H2_1(L)+(NCOLLH2-NH2_1(L))*(T_H2_1(L+1)-T_H2_1(L))   &
     &                            /(NH2_1(L+1)-NH2_1(L))               
                            IF ((AV(L).LE.TAU(jspace)).AND.(AV(L+1).GE.TAU(jspace)))                &
     &                            TETA3=T_AV(L)+(TAU(jspace)-AV(L))*(T_AV(L+1)-T_AV(L))     &
     &                            /(AV(L+1)-AV(L))                     
                     ENDDO                                             
                     
                                                 
                                   
                                                                      
!     VW Sept 08, we change the way the photodissociation is computed when the H2 and/or CO
!     column density and tau exceeds the largest value in Self-Shielding_data file, we use the smallest values of tetas. 
!     If H2 column density is smaller 
!     than the smallest value in the file, we use a photorate of 2.54e-11 (teta=1).
								      
                     IF (NCOLLH2.GT.NH2_1(43)) TETA1=T_H2_1(43)
                                       
                     IF (NCOLLCO.GT.NCO(52)) TETA2=T_CO(52)
                                       
                     IF (TAU(jspace).GT.AV(43)) TETA3=T_AV(43)
                                   
		     K(J,JSPACE)=1.03D-10*TETA1*TETA2*TETA3    		      

                      write(14,*)  NCOLLH2, NCOLLCO, K(J,JSPACE)
                                                 
                   ENDIF                                               
          ENDIF                                                        
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!	KOOIJ FORMULA	

	IF (FORMULA(J).EQ.3) 	then		                            
     	       K(J,jspace)=A(J)*(TT(jspace)**B(J))*EXP(-C(J)/TD(jspace))

!                  Check for temperature bounderies !!! Temperatures outside the T range are not extrapolated
		if (TD(jspace).LT.Tmin(J)) K(J,jspace)=A(J)*((Tmin(J)/300.D0)**B(J))*EXP(-C(J)/Tmin(J))
		if (TD(jspace).GT.Tmax(J)) K(J,jspace)=A(J)*((Tmax(J)/300.D0)**B(J))*EXP(-C(J)/Tmax(J))
		   
!                  Check for the presence of several rate coefficients present in the network for the 
!                  the same reaction
                IF (NUMBER(J+1).EQ.NUMBER(J)) THEN
			INDICE(W)=J
			W=W+1
			distmin(w)=tmin(j)-TD(jspace)
              		distmax(w)=TD(jspace)-tmax(j)
         		if (jstep.eq.0) write(*,*) 'Reaction ',NUMBER(J),' is present more than once in the network'
         	ENDIF

         	IF ((NUMBER(J+1).NE.NUMBER(J)).AND.(W.NE.1)) THEN
              		INDICE(W)=J
              		distmin(w)=tmin(j)-TD(jspace)
              		distmax(w)=TD(jspace)-tmax(j)

              		DO M=1,W
                		N=INDICE(M)
                		IF (TD(jspace).LT.Tmin(N)) K(N,jspace)=0.d0
                		IF (TD(jspace).GT.Tmax(N)) K(N,jspace)=0.d0
              		ENDDO

              		if (maxval(K(indice(1:w),jspace)).lt.1.d-99) then
              			if (minval(abs(distmin)).lt.minval(abs(distmax))) then
              				N=indice(minloc(abs(distmin),dim=1))
              				K(N,jspace)=A(N)*((Tmin(N)/300.D0)**B(N))*EXP(-C(N)/Tmin(N)) 
              			else
              				N=indice(minloc(abs(distmax),dim=1))
              				K(N,jspace)=A(N)*((Tmax(N)/300.D0)**B(N))*EXP(-C(N)/Tmax(N)) 
              			endif
              		endif

              		W=1
              		INDICE(:)=0
              		distmin(:)=9999.
              		distmax(:)=9999.
         	ENDIF       
         endif
 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
!	IONPOL1 and IONPOL2 FORMULA
	
 	IF ((FORMULA(J).EQ.4).OR.(FORMULA(J).EQ.5))	then		                            
     	        if (FORMULA(J).EQ.4)   K(J,jspace)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/TD(JSPACE))**0.5))
		if (FORMULA(J).EQ.5)   K(J,jspace)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TD(JSPACE))**0.5)+(C(J)**2*300.d0/(10.526*TD(JSPACE))))  

		
!                  Check for temperature bounderies
		if ((TD(jspace).LT.Tmin(J)).AND.(FORMULA(J).EQ.4)) K(J,jspace)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/Tmin(J))**0.5)) 
		if ((TD(jspace).GT.Tmax(J)).AND.(FORMULA(J).EQ.4)) K(J,jspace)=A(J)*B(J)*(0.62d0+0.4767d0*C(J)*((300.D0/TMAX(J))**0.5)) 
		if ((TD(jspace).LT.Tmin(J)).AND.(FORMULA(J).EQ.5)) &
     &	K(J,jspace)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMIN(J))**0.5)+(C(J)**2*300.d0/(10.526*TMIN(J)))) 
		if ((TD(jspace).GT.Tmax(J)).AND.(FORMULA(J).EQ.5)) &
     &	K(J,jspace)=A(J)*B(J)*((1.d0+0.0967*C(J)*(300.D0/TMAX(J))**0.5)+(C(J)**2*300.d0/(10.526*TMAX(J)))) 
   
		   
		   
!                  Check for the presence of several rate coefficients present in the network for the 
!                  the same reaction
                IF (NUMBER(J+1).EQ.NUMBER(J)) THEN
               			INDICE(W)=J
               			distmin(w)=tmin(j)-TD(jspace)
               			distmax(w)=TD(jspace)-tmax(j)
               			W=W+1
               			if (jstep.eq.0) write(*,*) 'Reaction ',NUMBER(J),' is present more than once in the network'
             	ENDIF
             	IF ((NUMBER(J+1).NE.NUMBER(J)).AND.(W.NE.1)) THEN
               		INDICE(W)=J
               		distmin(w)=tmin(j)-TD(jspace)
               		distmax(w)=TD(jspace)-tmax(j)

               		DO M=1,W
               			N=INDICE(M)
               			IF (TD(jspace).LT.Tmin(N)) K(N,jspace)=       0.d0
               			IF (TD(jspace).GT.Tmax(N)) K(N,jspace)=       0.d0
               		ENDDO
              
               		if (maxval(K(indice(1:w),jspace)).lt.1.d-99) then
               			if (minval(abs(distmin)).lt.minval(abs(distmax))) then
               				N=indice(minloc(abs(distmin),dim=1))
					if (FORMULA(N).EQ.4)   K(N,jspace)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/Tmin(N))**0.5))
					if (FORMULA(N).EQ.5)   &
     &					K(N,jspace)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMIN(N))**0.5)+(C(N)**2*300.d0/(10.526*TMIN(N))))
               			else
               				N=indice(minloc(abs(distmax),dim=1))
               				if (FORMULA(N).EQ.4)   K(N,jspace)=A(N)*B(N)*(0.62d0+0.4767d0*C(N)*((300.D0/TMAX(N))**0.5))
					if (FORMULA(N).EQ.5)   &
     &					K(N,jspace)=A(N)*B(N)*((1.d0+0.0967*C(N)*(300.D0/TMAX(N))**0.5)+(C(N)**2*300.d0/(10.526*TMAX(N))))
               			endif
               		endif

              		W=1
              		INDICE(:)=0
               		distmin(:)=9999.
               		distmax(:)=9999.
             	ENDIF   
			    
         ENDIF
	
        enddo                                                                                                                                      
!       remove the very small rate coefficients                         
       DO J=1,NRTOT 
              IF (K(J,JSPACE).LT.1.D-50) K(J,JSPACE)=0.D0 
       ENDDO 
                                                                        
!       writing the rate coefficients in the Kout.dat file              
       IF (JSTEP.EQ.0) THEN 
              WRITE (9,*) 'Rate coefficients of the reactions' 
              WRITE (9,*) 'Tk =    ',TD(JSPACE) 
              DO J=1,NRTOT 
                        WRITE (9,*) NUMBER(J), K(J,JSPACE) 
              ENDDO 
       ENDIF 
       enddo 
                                                                        
       RETURN 
      END                                           
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!       this subroutine contains the differential equations             
!                                                                       
       SUBROUTINE F (NEQ, T, Y, YDOT) 
                                                                        
       IMPLICIT NONE 
                                                                        
       include 'header.f90' 
                                                                        
       INTEGER NEQ,jspace 
                                                                        
       INTEGER IR1,IR2,IR3,DES1,DES2,IPROD1,IPROD2,IPROD3,IPROD4,I,J 
                                                                        
       DOUBLE PRECISION Y, YDOT, UP, DOWN, RATE,T,DY,DD,NHTOT2,      &
     &              RADIUS2                                            
       DIMENSION Y(NS,NX), YDOT(NS,NX), UP(NS,NX), DOWN(NS,NX),DD(NX),  &
     &              NHTOT2(NX),RADIUS2(NX),DY(NX)                       
                                                                        
       INTEGER REACT(NRTOT,7) 
       DOUBLE PRECISION K(NRTOT,NX),NHTOT(NX), RADIUS(2,NX) 
       COMMON/EQUA/K,NHTOT,RADIUS,REACT 
                                                                        
                                                                                                                                               
       DO JSPACE=1,NX 
              DO I=1,NS 
                     UP(I,JSPACE)=0.d0 
                     DOWN(I,JSPACE)=0.d0 
              ENDDO 
       ENDDO 
                                                                        
                                                                        
!       the differential equations are calcultaed in a loop here        
       DO I=1,NRTOT 
              IR1=REACT(I,1) 
              IR2=REACT(I,2) 
              IR3=REACT(I,3) 
              IPROD1=REACT(I,4) 
              IPROD2=REACT(I,5) 
              IPROD3=REACT(I,6) 
              IPROD4=REACT(I,7) 
                                                                        
              DO JSPACE=1,NX 
                      IF (IR3.NE.NS+1) THEN 
                      RATE=K(I,JSPACE)*Y(IR1,JSPACE)*Y(IR2,JSPACE)      &
     &                   *Y(IR3,JSPACE)*NHTOT(JSPACE)*NHTOT(JSPACE)     

                     IF (IPROD1.NE.NS+1) UP(IPROD1,JSPACE)=UP(IPROD1,JSPACE)+RATE 
                     IF (IPROD2.NE.NS+1) UP(IPROD2,JSPACE)=UP(IPROD2,JSPACE)+RATE 
                     IF (IPROD3.NE.NS+1) UP(IPROD3,JSPACE)=UP(IPROD3,JSPACE)+RATE 
                     IF (IPROD4.NE.NS+1) UP(IPROD4,JSPACE)=UP(IPROD4,JSPACE)+RATE 
                                                                        
                      DOWN(IR1,JSPACE)=DOWN(IR1,JSPACE)+RATE 
                      DOWN(IR2,JSPACE)=DOWN(IR2,JSPACE)+RATE 
                      DOWN(IR3,JSPACE)=DOWN(IR3,JSPACE)+RATE 

                     ENDIF 
                                                                        
                      IF ((IR3.EQ.NS+1).AND.(IR2.NE.NS+1)) THEN 
                      RATE=K(I,JSPACE)*Y(IR1,JSPACE)*Y(IR2,JSPACE)* &
     &                NHTOT(JSPACE)                                

                     IF (IPROD1.NE.NS+1) UP(IPROD1,JSPACE)=UP(IPROD1,JSPACE)+RATE 
                     IF (IPROD2.NE.NS+1) UP(IPROD2,JSPACE)=UP(IPROD2,JSPACE)+RATE 
                     IF (IPROD3.NE.NS+1) UP(IPROD3,JSPACE)=UP(IPROD3,JSPACE)+RATE 
                     IF (IPROD4.NE.NS+1) UP(IPROD4,JSPACE)=UP(IPROD4,JSPACE)+RATE 
                                                                         
                      DOWN(IR1,JSPACE)=DOWN(IR1,JSPACE)+RATE 
                      DOWN(IR2,JSPACE)=DOWN(IR2,JSPACE)+RATE 
 
 
                     ENDIF 
                                                                        
                      IF (IR2.EQ.NS+1) THEN 
                            RATE=K(I,JSPACE)*Y(IR1,JSPACE) 

                     IF (IPROD1.NE.NS+1) UP(IPROD1,JSPACE)=UP(IPROD1,JSPACE)+RATE 
                     IF (IPROD2.NE.NS+1) UP(IPROD2,JSPACE)=UP(IPROD2,JSPACE)+RATE 
                     IF (IPROD3.NE.NS+1) UP(IPROD3,JSPACE)=UP(IPROD3,JSPACE)+RATE 
                     IF (IPROD4.NE.NS+1) UP(IPROD4,JSPACE)=UP(IPROD4,JSPACE)+RATE 
                                                                         
                      DOWN(IR1,JSPACE)=DOWN(IR1,JSPACE)+RATE 
 
                     ENDIF 
                                                                                                                                                
              ENDDO 
                                                                        
       ENDDO 
                                                                        
                           
	do I=1,ns
		DO JSPACE=1,NX
			YDOT(I,JSPACE)=(UP(I,JSPACE)-DOWN(I,JSPACE))
		ENDDO	
	enddo
 
	                                                                                                  
       RETURN 
      END                                           
                                                                        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
!       this subroutine computes the rates of formation and destruction 
!       for each species and put them in the increasing order           
                                                                        
       SUBROUTINE COMP_RATES (Y,SPEC,NUMBER) 
                                                                        
       IMPLICIT NONE 
                                                                        
       INCLUDE 'header.f90' 
                                                                        
       INTEGER J,IR1,IR2,IR3,IPROD1,IPROD2,IPROD3,                      &
     &               IPROD4,IY,M,I,JSPACE,NUMBER                               
       DOUBLE PRECISION Y,FORMTOT,DESTOT,RATE,X,FORM,DEST 
       DIMENSION Y(NS,NX),FORMTOT(NS),DESTOT(NS),FORM(NS,NRTOT),        &
     &              DEST(NS,NRTOT),IY(NRTOT),X(NRTOT),NUMBER(NRTOT)                   
       CHARACTER*10 SPEC(NS+1) 
                                                                        
       INTEGER REACT(NRTOT,7) 
       DOUBLE PRECISION K(NRTOT,NX),NHTOT(NX),RADIUS(2,NX)
       COMMON/EQUA/K,NHTOT,RADIUS,REACT
                                                                        
       DO JSPACE=1,NX 
                                                                        
       DO J=1,NS 
              FORMTOT(J)=0.D0 
              DESTOT(J)=0.D0 
       ENDDO 
                                                                        
       DO I=1,NRTOT 
                                                                        
              IR1=REACT(I,1) 
              IR2=REACT(I,2) 
              IR3=REACT(I,3) 
              IPROD1=REACT(I,4) 
              IPROD2=REACT(I,5) 
              IPROD3=REACT(I,6) 
              IPROD4=REACT(I,7) 
	                                                                
               IF (IR3.NE.NS+1) RATE=K(I,JSPACE)*Y(IR1,JSPACE)*         &
     &              Y(IR2,JSPACE)*Y(IR3,JSPACE)                              
                                                                        
               IF ((IR3.EQ.NS+1).AND.(IR2.NE.NS+1))        &
     &              RATE=K(I,JSPACE)*Y(IR1,JSPACE)*Y(IR2,JSPACE)                                     
                                                                        
               IF (IR2.EQ.NS+1) RATE=K(I,JSPACE)*Y(IR1,JSPACE) 
                                                                        
                                                                        
              IF (IPROD1.NE.NS+1) FORM(IPROD1,I)=RATE 
              IF (IPROD2.NE.NS+1) FORM(IPROD2,I)=RATE 
              IF (IPROD3.NE.NS+1) FORM(IPROD3,I)=RATE 
              IF (IPROD4.NE.NS+1) FORM(IPROD4,I)=RATE 
              IF (IR1.NE.NS+1) DEST(IR1,I)=RATE 
              IF (IR2.NE.NS+1) DEST(IR2,I)=RATE 
              IF (IR3.NE.NS+1) DEST(IR3,I)=RATE 
                                                                        
              IF (IPROD1.NE.NS+1) FORMTOT(IPROD1)=FORMTOT(IPROD1)+RATE 
              IF (IPROD2.NE.NS+1) FORMTOT(IPROD2)=FORMTOT(IPROD2)+RATE 
              IF (IPROD3.NE.NS+1) FORMTOT(IPROD3)=FORMTOT(IPROD3)+RATE 
              IF (IPROD4.NE.NS+1) FORMTOT(IPROD4)=FORMTOT(IPROD4)+RATE 
              IF (IR1.NE.NS+1) DESTOT(IR1)=DESTOT(IR1)+RATE 
              IF (IR2.NE.NS+1) DESTOT(IR2)=DESTOT(IR2)+RATE 
              IF (IR3.NE.NS+1) DESTOT(IR3)=DESTOT(IR3)+RATE 
                                                                        
       ENDDO 
                                                                        
!       put FORM and DEST in the increasing order                       
       DO M=1,NS 
       DO I=1,NRTOT 
              X(I)=FORM(M,I) 
              IY(I)=NUMBER(I) 
       ENDDO 
                                                                        
       CALL SSORT (X, IY, NRTOT) 
                                                                        
       WRITE (11,*) SPEC(M) 
       WRITE (11,*) 'REACTIONS OF PRODUTION    Total flux (cm-3s-1) =',FORMTOT(M) 
       DO I=1,5 
              IF (X(I).NE.0D0) WRITE(11,*)  IY(I), X(I) 
       ENDDO 
                                                                        
       DO I=1,NRTOT 
              X(I)=DEST(M,I) 
              IY(I)=NUMBER(I) 
       ENDDO 
                                                                        
       CALL SSORT (X, IY, NRTOT) 
                                                                        
       WRITE (11,*) 'REACTIONS OF DESTRUCTION    Total flux (cm-3s-1) =',DESTOT(M) 
       DO I=1,5 
              IF (X(I).NE.0D0) WRITE(11,*)  IY(I), -X(I) 
       ENDDO 
                                                                        
       ENDDO 
                                                                        
       ENDDO 
                                                                        
       RETURN 
      END                                           
                                                                        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
!       this subroutine checks that the densities of charged species    
!       and elements are conservative                                   
                                                                        
       SUBROUTINE CHECKING (ELEMENT,SN,IGRAIN0,IGRAINN) 
                                                                        
       IMPLICIT NONE 
                                                                        
       include 'header.f90' 
                                                                        
       INTEGER ELEMENT,I,J,IGRAIN0,IGRAINN,jspace 
       DOUBLE PRECISION SN,PLUS,MINUS,DPM,CHECK_ELEM 
       DIMENSION SN(NS,NX),ELEMENT(NELEM,NS+1),                           &
     &              CHECK_ELEM(NELEM-1)                                 
                                                                        
       PLUS=0.D0 
       MINUS=0.D0 
                                                                        
       DO JSPACE=1,NX 
       WRITE(*,*) 'jspace ',JSPACE 
                                                                        
       DO I=1,NELEM-1 
              CHECK_ELEM(I)=0.D0 
       ENDDO 
                                                                        
       DO I=1,NS 
              IF (ELEMENT(1,I).EQ.1) PLUS=PLUS+SN(I,JSPACE) 
              IF (ELEMENT(1,I).EQ.-1) MINUS=MINUS+SN(I,JSPACE) 
       ENDDO 
                                                                        
       DPM=(PLUS-MINUS)/PLUS 
                                                                        
       WRITE (*,*) '(positive) - (negative)/(positive)',DPM 
       WRITE (*,*) 'grains',SN(IGRAIN0,JSPACE)+SN(IGRAINN,JSPACE) 
       DO I=1,NELEM-1 
              DO J=1,NS 
                     CHECK_ELEM(I)=CHECK_ELEM(I)+                       &
     &                     FLOAT(ELEMENT(I+1,J))*SN(J,JSPACE)           
              ENDDO 
       ENDDO 
                                                                        
       WRITE(11,*) '( H  He C  N  O  SI S  FE NA MG CL P F)' 
       WRITE(11,*) CHECK_ELEM 
       WRITE(*,*) '( H  He C  N  O  SI S  FE NA MG CL P F)' 
       WRITE(*,*) CHECK_ELEM 
       enddo 
                                                                        
       RETURN 
      END                                           
                                                                        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
!       dummy subroutine                                                
!       this subroutine is needed because we are not supplying          
!       the jacobian matrix                                             
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
                                                                        
       SUBROUTINE DUMMY 
       IMPLICIT NONE 
       INTEGER N,J 
       DOUBLE PRECISION T,Y,IAN,JAN,PDJ 
                                                                        
       ENTRY JAC (N, T, Y, J, IAN, JAN, PDJ) 
                                                                        
       RETURN 
      END                                           
                                                                        
!       SUBROUTINE DUMMY                                                
                                                                        
!       ENTRY JAC (NEQ, T, Y, ML, MU, PD, NROWPD)                       
                                                                        
!       RETURN                                                          
!       END                                                             
                                                                        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
!       subroutine to sort the numbers                                  
                                                                        
      SUBROUTINE SSORT (X, IY, N) 
      IMPLICIT NONE 
                                                                        
!                                                                       
!    Example of an Insertion Sort                                       
!                                                                       
!***BEGIN PROLOGUE  SSORT                                               
!***PURPOSE  Sort an array and make the same interchanges in            
!            an auxiliary array.  The array is sorted in                
!            decreasing order.                                          
!***TYPE      SINGLE PRECISION                                          
!***KEYWORDS  SORT, SORTING                                             
!                                                                       
!   Description of Parameters                                           
!      X - array of values to be sorted   (usually abscissas)           
!      IY - array to be carried with X (all swaps of X elements are     
!          matched in IY .  After the sort IY(J) contains the original  
!          postition of the value X(J) in the unsorted X array.         
!      N - number of values in array X to be sorted                     
!      KFLAG - Not used in this implementation                          
!                                                                       
!***REVISION HISTORY  (YYMMDD)                                          
!   950310  DATE WRITTEN                                                
!   John Mahaffy                                                        
!***END PROLOGUE  SSORT                                                 
!     .. Scalar Arguments ..                                            
      INTEGER N 
!     .. Array Arguments ..                                             
      DOUBLE PRECISION X(*) 
      INTEGER IY(*) 
!     .. Local Scalars ..                                               
      DOUBLE PRECISION TEMP 
      INTEGER I, J, L, ITEMP 
!     .. External Subroutines ..                                        
!     None                                                              
!     .. Intrinsic Functions ..                                         
!     None                                                              
!                                                                       
!***FIRST EXECUTABLE STATEMENT  SSORT                                   
!                                                                       
      DO 100 I=2,N 
         IF ( X(I).GT.X(I-1) ) THEN 
            DO 50 J=I-2,1,-1 
              IF(X(I).LT.X(J)) go to 70 
   50         CONTINUE 
            J=0 
   70       TEMP=X(I) 
            ITEMP=IY(I) 
            DO 90 L=I,J+2,-1 
              IY(L)=IY(L-1) 
   90         X(L)=X(L-1) 
            X(J+1)=TEMP 
            IY(J+1)=ITEMP 
         ENDIF 
  100 END DO 
      RETURN 
      END                                           

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       SUBROUTINE JACFH(NEQ,T,Y, ZJ,IAN,JAN, PDJ)

       IMPLICIT NONE

       include 'header.f90'

       INTEGER NEQ,jspace

       INTEGER IR1,IR2,IR3,DES1,DES2,IPROD1,IPROD2,IPROD3,IPROD4,I,J, ZJ,ipdj,jpdj

       DOUBLE PRECISION Y, YDOT, UP, DOWN, RATE,T,DD,nu2
       DIMENSION Y(NS,NX), YDOT(NS,NX), UP(NS,NX), DOWN(NS,NX),DD(NX)
       real(kind=8),dimension(ns,nx) :: PDJ
       real(kind=8), dimension(nx) :: DY, radius2
       real(kind=8), dimension(2,nx) :: radius
       REAL(KIND=8), dimension(NEQ) :: IAN, JAN


       real(kind=8), dimension(1:nx) :: ajac, bjac, cjac


       INTEGER REACT(NRTOT,7)
       DOUBLE PRECISION K(NRTOT,NX),NHTOT(NX)
       COMMON/EQUA/K,NHTOT,radius, REACT


! Find the relevant i and j corresponding to zj

       ipdj=mod(zj,ns)
       jpdj=int(zj/ns)+1
       if (ipdj.eq.0) then
       ipdj=ns
       jpdj=jpdj-1
       endif

! Initialize PDJ

       PDJ=0.d0

! Flee if ipdj=NS
!       if (ipdj.eq.NS) return


       UP=0.d0
       DOWN=0.d0


!       the differential equations are calcultaed in a loop here
       DO I=1,NRTOT
              IR1=REACT(I,1)
              IR2=REACT(I,2)
              IR3=REACT(I,3)
              IPROD1=REACT(I,4)
              IPROD2=REACT(I,5)
              IPROD3=REACT(I,6)
              IPROD4=REACT(I,7)

                if (IR3.ne.NS+1) then
!                RATE=K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)

       if (IR1.eq.IPDJ) then
       IF (IPROD1.NE.NS+1) PDJ(IPROD1,jpdj)=PDJ(IPROD1,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD2.NE.NS+1) PDJ(IPROD2,jpdj)=PDJ(IPROD2,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD3.NE.NS+1) PDJ(IPROD3,jpdj)=PDJ(IPROD3,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD4.NE.NS+1) PDJ(IPROD4,jpdj)=PDJ(IPROD4,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR1,jpdj)=PDJ(IR1,jpdj)-K(I,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR2,jpdj)=PDJ(IR2,jpdj)-K(I,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR3,jpdj)=PDJ(IR3,jpdj)-K(I,JPDJ)*Y(IR2,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       endif

       if (IR2.eq.IPDJ) then
       IF (IPROD1.NE.NS+1) PDJ(IPROD1,jpdj)=PDJ(IPROD1,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD2.NE.NS+1) PDJ(IPROD2,jpdj)=PDJ(IPROD2,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD3.NE.NS+1) PDJ(IPROD3,jpdj)=PDJ(IPROD3,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD4.NE.NS+1) PDJ(IPROD4,jpdj)=PDJ(IPROD4,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR1,jpdj)=PDJ(IR1,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR2,jpdj)=PDJ(IR2,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR3,jpdj)=PDJ(IR3,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR3,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       endif

       if (IR3.eq.IPDJ) then
       IF (IPROD1.NE.NS+1) PDJ(IPROD1,jpdj)=PDJ(IPROD1,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD2.NE.NS+1) PDJ(IPROD2,jpdj)=PDJ(IPROD2,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD3.NE.NS+1) PDJ(IPROD3,jpdj)=PDJ(IPROD3,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       IF (IPROD4.NE.NS+1) PDJ(IPROD4,jpdj)=PDJ(IPROD4,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR1,jpdj)=PDJ(IR1,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR2,jpdj)=PDJ(IR2,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       PDJ(IR3,jpdj)=PDJ(IR3,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)*NHTOT(JPDJ)
       endif

                endif
                if ((IR3.eq.NS+1).and.(IR2.ne.NS+1)) then
!                RATE=K(I,JPDJ)*Y(IR1,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)

       if (IR1.eq.IPDJ) then
       IF (IPROD1.NE.NS+1) PDJ(IPROD1,jpdj)=PDJ(IPROD1,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)
       IF (IPROD2.NE.NS+1) PDJ(IPROD2,jpdj)=PDJ(IPROD2,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)
       IF (IPROD3.NE.NS+1) PDJ(IPROD3,jpdj)=PDJ(IPROD3,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)
       IF (IPROD4.NE.NS+1) PDJ(IPROD4,jpdj)=PDJ(IPROD4,jpdj)+K(I,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)
       PDJ(IR1,jpdj)=PDJ(IR1,jpdj)-K(I,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)
       PDJ(IR2,jpdj)=PDJ(IR2,jpdj)-K(I,JPDJ)*Y(IR2,JPDJ)*NHTOT(JPDJ)
       endif

       if (IR2.eq.IPDJ) then
       IF (IPROD1.NE.NS+1) PDJ(IPROD1,jpdj)=PDJ(IPROD1,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*NHTOT(JPDJ)
       IF (IPROD2.NE.NS+1) PDJ(IPROD2,jpdj)=PDJ(IPROD2,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*NHTOT(JPDJ)
       IF (IPROD3.NE.NS+1) PDJ(IPROD3,jpdj)=PDJ(IPROD3,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*NHTOT(JPDJ)
       IF (IPROD4.NE.NS+1) PDJ(IPROD4,jpdj)=PDJ(IPROD4,jpdj)+K(I,JPDJ)*Y(IR1,JPDJ)*NHTOT(JPDJ)
       PDJ(IR1,jpdj)=PDJ(IR1,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*NHTOT(JPDJ)
       PDJ(IR2,jpdj)=PDJ(IR2,jpdj)-K(I,JPDJ)*Y(IR1,JPDJ)*NHTOT(JPDJ)
       endif

       endif

       if (IR2.eq.NS+1) then

       if (IR1.eq.IPDJ) then
       IF (IPROD1.NE.NS+1) PDJ(IPROD1,jpdj)=PDJ(IPROD1,jpdj)+K(I,JPDJ)
       IF (IPROD2.NE.NS+1) PDJ(IPROD2,jpdj)=PDJ(IPROD2,jpdj)+K(I,JPDJ)
       IF (IPROD3.NE.NS+1) PDJ(IPROD3,jpdj)=PDJ(IPROD3,jpdj)+K(I,JPDJ)
       IF (IPROD4.NE.NS+1) PDJ(IPROD4,jpdj)=PDJ(IPROD4,jpdj)+K(I,JPDJ)
       PDJ(IR1,jpdj)=PDJ(IR1,jpdj)-K(I,JPDJ)
       endif

                endif
       ENDDO

       RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine computeIAJA(IA,JA,liw,SPEC,Y)
      implicit none
      include 'header.f90'
      integer :: liw
      integer, dimension(NS*NX+1) :: IA
      integer, dimension(liw) :: JA
      real(kind=8), dimension(NS,NX) :: DUMMYY, Y
      real(kind=8), dimension(NS*NX) :: PDJ
      CHARACTER(len=8), dimension(NS) :: SPEC

      integer :: i,j,kjac
      real(kind=8) :: T, hop
      REAL(KIND=8), dimension(NS*NX) :: IAN, JAN

      INTEGER REACT(NRTOT,7)
      DOUBLE PRECISION K(NRTOT,NX),NHTOT(NX)
       real(kind=8), dimension(2,nx) :: radius
      COMMON/EQUA/K,NHTOT,radius, react

      DUMMYY(:,1)=1.d-5
      do i=1,NX-1
      DUMMYY(:,I+1)=DUMMYY(:,I)*1.1d0
      enddo

      kjac=1

      do j=1,NS*NX
      call JACFH(NS*NX,T,Y,j,IAN,JAN,PDJ)

!      print *, j, maxval(PDJ)

      IA(j)=kjac

      hop=1.d-4*maxval(abs(PDJ))
!      hop=0.D0
      do i=1,NS*NX
!      if (PDJ(i).ne.0.d0) then
!      if (abs(PDJ(i)).gt.1.d-9) then
      if (abs(PDJ(i)).gt.hop) then
      JA(kjac)=i
      kjac=kjac+1
      endif
      enddo

      enddo

!      stop

      IA(NS*NX+1)=kjac

      return
      end


