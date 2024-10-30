echo "#input parameters
===============================================================================
nter=1
verb=y
vpfile=vini //velocity file
shots=shots
wtfile=wt //wavelet file
nt=2000 //number of time steps
dt=0.001 //temporal sampling
lft=30 //left boundary
rht=30 //right boundary
top=30 //top boundary
bot=30 //bottom boundary
n1=200 //nz of input FD model
n2=300 //nx of input FD model
d1=10 //dz spacing of the input FD model
d2=10 //dx spacing of the input FD model
jsx=15 
jrx=1
jsz=0
jrz=0
sx=0
rx=0
sz=5
rz=0
ns=20
nr=300
mode=0 //type of EAL 0,EAL 1,MEAL
out=vinv //wavefield output file

" >inputpar.txt

../bin/fwi $(cat inputpar.txt)
