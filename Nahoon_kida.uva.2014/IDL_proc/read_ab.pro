;	pro read_ab

; parameters to modify:
nrun=100		; number of runs
nt=124			; number of times
rep='ab_spec'		; directory where the species.dat files are

nt=nt-1
aa=' '
n=0.d0
tyr=fltarr(nt)
ab=fltarr(nrun,nt)
nn=fltarr(nrun+1)

spec_selec=' '
read,'Which species ?',spec_selec

char1=strcompress(rep+'/'+spec_selec+'.dat', /remove_all)
	
openr,1,char1
for j=0,nt-1 do begin
	readf,1,nn
	tyr(j)=nn(0)
	ab(*,j)=nn(1:nrun)
endfor
close,1
range=[min(alog10(ab(*,1:nt-1)))-0.1,max(alog10(ab(*,1:nt-1)))+0.1]

plot,alog10(tyr(1:nt-1)),alog10(ab(0,1:nt-1)),xtitle='log(time (yr))',$
ytitle='log(Abundance)',yrange=range,ystyle=1,charsize=2.5
for i=1,nrun-1 do oplot,alog10(tyr(1:nt-1)),alog10(ab(i,1:nt-1))

end
