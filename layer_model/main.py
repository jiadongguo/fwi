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
    nr=int(dic["nr"])
    nt=int(dic["nt"])
    ns=int(dic["ns"])
    lft=int(dic["lft"])
    rht=int(dic["rht"])
    top=int(dic["top"])
    bot=int(dic["bot"])
    vel=np.fromfile(dic["vpfile"],np.float32,nx*nz)
    vel=vel.reshape((nz,nx),order='F')
    rcd=np.fromfile(dic["out"],np.float32,nr*nt,offset=10*4*nr*nt)
    rcd=rcd.reshape((nt,nr),order='F')
    plt.figure(1)
    #plt.title("速度模型")
    plt.imshow(vel,cmap="jet",aspect="auto")
    plt.figure(2)
    #plt.title("波场记录")
    plt.imshow(rcd,cmap="bwr",aspect="auto",vmin=-1,vmax=1)

    plt.show()
