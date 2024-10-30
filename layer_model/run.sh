#!/bin/bash
vpfile=vel #velocity file
wtfile=wt #wavelet file
nt=2000 #number of time steps
dt=0.001 #temporal sampling
lft=30 #left boundary
rht=30 #right boundary
top=30 #top boundary
bot=30 #bottom boundary
n1=200 #nz of input FD model
n2=300 #nx of input FD model
d1=10 #dz spacing of the input FD model
d2=10 #dx spacing of the input FD model
jsx=30 
nzx=`expr $n1 \* $n2`
jrx=1
jsz=0
jrz=0
sx=0
rx=0
sz=5
rz=0
ns=10
nr=300
mode=0 #type of EAL 0,EAL 1,MEAL
out=wfd #wavefield output file
binfile='../bin/'

#mpirun -n 20 $binfile/fdmodeling2 vpfile=$vpfile wtfile=$wtfile nt=$nt dt=$dt lft=$lft rht=$rht top=$top bot=$bot n1=$n1 n2=$n2 d1=$d1 d2=$d2\
                                jsx=$jsx jrx=$jrx jsz=$jsz jrz=$jrz sx=$sx rx=$rx sz=$sz rz=$rz ns=$ns nr=$nr node=$mode out=$out 
#python main.py $(cat inputpar.txt)
