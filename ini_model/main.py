#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys
import struct
import numpy as np
import matplotlib.pyplot as plt  
import matplotlib as mpl  



if __name__ == "__main__":
    print('脚本名为：', sys.argv[0])
    dic = {}
    for i in range(1, len(sys.argv)):
        if "=" in sys.argv[i]:
            k=sys.argv[i].split('=')
            dic[k[0]]=k[1]

    nx=int(dic["n2"])
    nz=int(dic["n1"])  
    vel=np.fromfile(dic["vpfile"],np.float32,nx*nz)
    vel=vel.reshape((nz,nx),order='F')
    vini=np.fromfile(dic["out"],np.float32,nx*nz)
    vini=vini.reshape((nz,nx),order='F')
    plt.figure(1)
    #plt.title("速度模型")
    plt.imshow(vel,cmap="jet",aspect="auto",vmin=1500,vmax=3500)
    plt.figure(2)
    #plt.title("波场记录")
    plt.imshow(vini,cmap="jet",aspect="auto",vmin=1500,vmax=3500)

    plt.show()
