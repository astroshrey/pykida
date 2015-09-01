;	pro read_plots



read,'ps (0) ou ecran (1)',aff
if aff eq 0 then begin
set_plot,'ps'
device,filename='tot.ps',scale_factor=1,/landscape
endif 


spec_select1='HCS+'

aa=' '
n=0.d0
nt=121
tyr=fltarr(nt)
ns=454
nx=1
ab=fltarr(nx,ns,nt)
nn=fltarr(nt)
spec=strarr(ns)
tau=fltarr(nx)
tau2=fltarr(nx)
n2=fltarr(5)
ligne=fltarr(4)
table=fltarr(4,nt)

char1=strcompress('Error/'+spec_select1+'.dat', /remove_all)	
	
openr,1,char1
print,'Opening ',char1
for j=0,nt-1 do begin
	readf,1,ligne
	table(*,j)=ligne(*)
endfor
close,1


ispec1=where(spec eq spec_select1)

plot,table(0,*),table(1,*),xtitle='log(time (yr))',$
ytitle='log(Abundance)',ystyle=1,charsize=2.5,thick=2,$
charthick=2,xrange=[2,7],yrange=[-14,-10.5]
oplot,table(0,*),smooth(table(2,*),3),linestyle=1
oplot,table(0,*),smooth(table(3,*),3),linestyle=1


if aff eq 0 then begin
device,/close
set_plot,'x'
endif


end
