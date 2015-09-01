aa=' '
n=0.d0
aa2=' '

nrun=1999
ns=474
ntime=123
spec=strarr(ns)	
ab=dblarr(ntime,ns)
tyr=dblarr(ntime)
n2=dblarr(nrun)
mean_ab=dblarr(ntime)
stand_dev=dblarr(ntime)
istart=0
close,1
close,2

spec=strarr(ns)
openr,1,'cond_initial.dat'
for i=0,0 do readf,1,aa
for i=0,ns-1 do begin
	readf,1,aa,format='(5x,a10)'
	spec(i)=aa
endfor
close,1


for i=istart,ns-1 do begin
;for i=istart,istart do begin
	mean_ab(*)=0.d0
	stand_dev(*)=0.d0

	char1=string(i+1,format='(i4)')
	char2=strcompress('ab_spec/'+spec(i)+'.dat', /remove_all)

	print,'Opening ',char2
;	print,'nombre de lignes:',file_lines(char2)
;	aaa=''
;	read,'?',aaa
	openr,1,char2
	n2(*)=0.d0
	readf,1,format='(e10.3,2499(e10.3))',n,n2
;	readf,1,n,n2
	for k=0,122 do begin
		n=0.d0
		n2(*)=0.d0
		readf,1,format='(e10.3,2499(e10.3))',n,n2
;		readf,1,n,n2
;		print,'toto',k
;		print,n
;		print,n2
		n2=alog10(n2)
		tyr(k)=alog10(n)
		mean_ab(k)=mean(n2)
		for j=0,nrun-1 do stand_dev(k)=stand_dev(k)+(mean_ab(k)-n2(j))*(mean_ab(k)-n2(j))
		stand_dev(k)=sqrt(stand_dev(k)/(nrun-1))
	endfor	
	close,1
	
	
	char3=strcompress('stand_dev/'+spec(i)+'.dat', /remove_all)
	openw,2,char3
	print,'Writing ',char3
	for k=0,ntime-1 do printf,2,tyr(k),mean_ab(k),stand_dev(k)
	close,2

endfor			

	
end	
