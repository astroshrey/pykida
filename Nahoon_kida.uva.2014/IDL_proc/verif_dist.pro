;	pro verif_dist
;
; This IDL procedure will help you checking the 
; distribution of the rate coefficients and gas 
; temperature and density

aa=' '
q1=' '	&	q2=q1	&	q3=q1
r=0.d0
bin1=0.d0
n=0
iK=0
nK=0
nrun=0
rep=' '

read,'Do you want to check the distribution of k(i) (Y or N)',q1

read,'Do you want to check the distribution of the temperature (Y or N)',q2

read,'Do you want to check the distribution of H density (Y or N)',q3

read,'Number of runs ?',nrun

read,'Repertory where the files are ?',rep

if q1 eq 'Y' then n=n+1
if q2 eq 'Y' then n=n+1
if q3 eq 'Y' then n=n+1

!p.multi=[0,1,n]

if q1 eq 'Y' then begin
	read,'Number of reactions ?',nK
	read,'Which reaction do you want to see ?',iK
	
	k=dblarr(nK,nrun)

	for i=0,nrun-1 do begin
		char3=strcompress(rep+'/Kout'+string(i)+'.dat', /remove_all)
		print,'Opening ',char3
	
		openr,1,char3
		readf,1,aa
		readf,1,aa
		for j=0,nK-1 do begin
			readf,1,n2,n1
			K(j,i)=alog10(n1)
		endfor
		close,1
	endfor



bin1=(max(K(iK-1,*))-min(K(iK-1,*)))/9.
histogram_ez,K(iK-1,*),binsize=bin1,xtitle='log(k)',charsize=2.5,$
xrange=[min(K(iK-1,*))-0.1,max(K(iK-1,*))+0.1],xstyle=1

endif

if q2 eq 'Y' or q3 eq 'Y' then begin
	rep2=rep
;	read,'Repertory where the plotxxxx.dat files are ?',rep2
	
	T=dblarr(nrun)
	nHtot=dblarr(nrun)
	
	for i=0,nrun-1 do begin
		char3=strcompress(rep+'/plot'+string(i)+'.dat', /remove_all)
		print,'Opening ',char3
	
		openr,1,char3
		readf,1,format='(23x,e14.8)',r
		T(i)=r
		readf,1,format='(23x,e14.8)',r
		nHtot(i)=r
		close,1
	endfor
print,T
if q2 eq 'Y' then begin
bin1=(max(T)-min(T))/9.
histogram_ez,T,binsize=bin1,xtitle='T (K)',charsize=2.5
endif

if q3 eq 'Y' then begin
bin1=(max(nHtot)-min(nHtot))/9.
histogram_ez,nHtot,binsize=bin1,xtitle='nH cm-3',charsize=2.5
endif

endif

end
 
