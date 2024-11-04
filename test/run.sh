#!/bin/bash
# 真实速度模型、初始速度模型、地震子波、真实地震记录、目标函数在当前文件夹
# 
export PATH=../bin:$PATH # 临时环境变量

prefix="./dat" #数据文件夹
suffix=".dat"  #文件后缀
read -p "please check ${prefix} is empty[y/n]:" ok
if [[ $ok != "y" ]];then
    echo "stop"
    exit
elif [ -e ${prefix} ];then
    echo "exist ${prefix}"
    #rm -rf ${prefix}
else
    mkdir ${prefix}
fi
niter=10 #反演迭代次数
nodes=4 #mpi节点数
# 网格参数
nz=92
nx=334
dz=10
dx=10
top=30
bot=30
lft=30
rht=30
# 炮点参数
ns=20
sz=5
sx=0
jsx=15
jsz=0
# 检波点参数
nr=334
rz=10
rx=0
jrx=1
jrz=0
# 子波参数
fm=10
amp=10
nt=2000
dt=0.001
t0=0.1
# 平滑参数
rect1=5
rect2=5
repeat=3
# 文件定义

wt=ricker${suffix} # 子波
vp=vel # 速度模型
dobs=dobs${suffix} # 观测地震记录
#地震子波
# ricker n1=$nt d1=$dt fm=$fm amp=$amp t0=$t0 out=$wt
# 制作观测地震记录
# mpirun -np $nodes fdmodeling2 n1=$nz n2=$nx nt=$nt dt=$dt top=$top bot=$bot lft=$lft rht=$rht ns=$ns nr=$nr sz=$sz \
#         sx=$sx jsx=$jsx jsz=$jsz rx=$rx rz=$rz jrx=$jrx jrz=$jrz d1=$dz d2=$dx vpfile=$vp \
#         wtfile=$wt out=$dobs mode=0
# 制作初始速度模型
vini="vini${suffix}"
# smooth repeat=$repeat cut=10 n1=$nz n2=$nx r1=$rect1 r2=$rect2 vpfile=$vp out=$vini

# 开始反演
vinv="vinv${suffix}"
mpirun -np $nodes fwi verb=y niter=$niter fobs=$dobs n1=$nz n2=$nx nt=$nt dt=$dt top=$top bot=$bot lft=$lft rht=$rht ns=$ns nr=$nr sz=$sz \
        sx=$sx jsx=$jsx jsz=$jsz rx=$rx rz=$rz jrx=$jrx jrz=$jrz d1=$dz d2=$dx vpfile=$vini \
        wtfile=$wt out=$vinv mode=0

