echo "#input parameters
===============================================================================
vpfile=vel //velocity file
n1=200 //nz of input FD model
n2=300 //nx of input FD model
out=vini //wavefield output file
r1=10
r2=10
repeat=1
" >inputpar.txt

../bin/smooth $(cat inputpar.txt)
python main.py $(cat inputpar.txt)
