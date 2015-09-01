; pro rewrite_kida_species

kida_file='../gas_species_kida.uva.2014.in'
output_file='../cond_initial.dat'
ns=489
aa=' '
ii=intarr(14)
d=0.

openr,1,kida_file
openw,2,output_file
printf,2,'JSPACE = 0'
for i=0,ns-1 do begin
	readf,1,format='(a11,14(i3))',aa,ii
	printf,2,format='(i4,1x,a10,14(i3),e16.8)',i+1,aa,ii,d
endfor
close,1
close,2

end

