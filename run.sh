#!/bin/bash
nz=200
nx=300
nzx=`expr $nz \* $nx`
nter=10
vel='vini'
for((i=1;i<=$nter;i++))
do
    grad='grad'$i
    rcd='rcd'$i
    dir='dir'$i
    velnew='vel'$i
    alpha='alpha'$i
    #forward modeling
    #objective function
    #updating direction
    #update
    echo update n=$nzx fdir=$dir out=$velnew in=$vel falpha=$alpha
    vel=$velnew
done