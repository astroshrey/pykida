;	pro read_plots

aa=' '
n=0.d0
nt=124
tyr=dblarr(nt)
ns=454
ab=dblarr(ns,nt)
nn=dblarr(nt)
spec=strarr(ns)


char1='../plot.dat'
	
openr,1,char1
print,'Opening ',char1
readf,1,aa
readf,1,aa
readf,1,format='(23x,e10.2)',n
readf,1,tyr
for j=0,ns-1 do begin
	readf,1,format='(a10,124(d16.8))',aa,nn
	spec(j)=strcompress(aa, /remove_all)
	ab(j,*)=nn(*)	
endfor
	close,1

spec_selec=' '
read,'Which species ?',spec_selec

ispec=where(spec eq spec_selec)

range=[min(alog10(ab(ispec,1:nt-1)))-0.1,max(alog10(ab(ispec,1:nt-1)))+0.1]

plot,alog10(tyr(1:nt-1)),alog10(ab(ispec,1:nt-1)),xtitle='log(time (yr))',$
ytitle='log(Abundance)',yrange=range,ystyle=1,charsize=2.5,thick=2,$
charthick=2


end
