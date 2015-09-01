!       NS = NUMBER OF SPECIES
!       NRTOT = NUMBER OF REACTIONS IN THE NETWORK (THE REAL ONE!!!)
!       NX = NUMBER OF SPATIAL POINTS
!       NTIME = NUMBER OF OUTPUT TIMES
!       NELEM = NUMBER OF ELEMENTS + 1 (FOR CHARGES)
!       XCO_0 and XH2_0 are the abundances of CO and H2 at the edge of 
!       the cloud to compute shelf sheilding of CO and H2. ONE NEEDS TO PAY
!       ATTENTION TO THESE PARAMETERS SINCE THEY CAN STRONGLY AFFECT THE 
!       CHEMISTRY
!       
       integer NS, NRTOT, NX, NTIME, NELEM 
       double precision XCO_0,XH2_0
       parameter (NS=489, NRTOT=7499, NTIME=124,NX=1,NELEM=14)	
       parameter (XCO_0=1.D-4, XH2_0=5.D-1)
