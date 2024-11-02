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
nter=10 #反演迭代次数
nodes=20 #mpi节点数
# 网格参数
nz=200
nx=300
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
nr=300
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

# 文件定义

wt=ricker${suffix} # 子波
vp=vel.dat #速度模型
dobs=dobs${suffix} # 观测地震记录
#地震子波
ricker n1=$nt d1=$dt fm=$fm amp=$amp t0=$t0 out=$wt
# 制作观测地震记录
mpirun -np $nodes fdmodel n1=$nz n2=$nx nt=$nt dt=$dt top=$top bot=$bot lft=$lft rht=$rht ns=$ns nr=$nr sz=$sz \
        sx=$sx jsx=$jsx jsz=$jsz rx=$rx rz=$rz jrx=$jrx jrz=$jrz d1=$dz d2=$dx vpfile=$vp \
        wtfile=$wt out=$dobs
exit
# 目标函数文件
fobj="obj.txt"
if [ -e ${fobj} ]&&[ -f ${fobj} ];then
    echo "remove ${fobj}"
    rm ${fobj}
elif [ -d ${fobj} ];then
    echo "${fobj} is directory"
    exit
else
    echo "create ${fobj}"
    touch ${fobj}
fi
# 开始反演
vini="vini${suffix}"
echo $vini
for((i=1;i<${nter};i++))
do
    echo "starting iteration$i"
    dcal="${prefix}/dcal$i${suffix}"
    dcal="${prefix}/dcaltest$i${suffix}"
    grad="${prefix}/grad$i${suffix}"
    veltest="${prefix}/veltest$i${suffix}"
# 正演计算dcal
# 计算目标函数
    fwi_obj ns=$ns nt=$nt ng=$nr out=$fobj fcal=$dcal fobs=$dobs
# 计算负梯度
# 计算测试模型
# 计算测试模型dcal
# 更新模型
    vini="${prefix}/velnew${i}${suffix}"
done
