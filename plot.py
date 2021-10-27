#!/usr/bin/env python3

"""
Script to plot data generated with Worm.x.
The output file can be given as an argument, otherwise
it must be explicitly entered. Note that the displayed
error bars are the one computed on the fly by Worm.x, and
are just for rough reference; proper statistical analysis
(reblocking) must be performed separately!
"""

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from typing import Tuple


def LoadData(quantity:str,filename:str) -> npt.ArrayLike:
    """
    Reads data regarding the phyisical observable 'quantity' from a
    Worm.x output file, and stores it in a numpy array with the block
    value, the overall average and the raw errorbar'.
    """
    os.system('grep ' + quantity + ' ' + filename +' > filetemp.tmp')
    data = np.loadtxt('filetemp.tmp',usecols=(1,2,3),unpack=False)
    os.system('rm filetemp.tmp')
    return data


def PlotTimeSeries(quantity:str,filename:str,x_label:str,
    y_label:str,y_lim:Tuple=None,h_line:bool=False) -> None:
    """"
    Loads data using LoadData and plot a physical quantity as
    a function of the MC block; it shows the block value and
    the overall average, with the raw error bar. A horizontal
    line at y=1 can be added, as a reference e.g. for the SF
    fraction.
    """
    data = LoadData(quantity,filename)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if y_lim:
        plt.ylim(y_lim)
    plt.plot(range(len(data)),data[:,0])
    plt.errorbar(range(len(data)),data[:,1],data[:,2])
    if h_line:
        plt.axhline(y=1.0,color='black',linestyle='-')
    plt.show()


def Plot_XY(filename:str,title:str) -> None:
    data = np.loadtxt(filename,usecols=(0,1,2))
    """
    Reads and plots data in the xy plane, showing a colormap
    (useful e.g. for density).
    """
    x = []
    for aa in data[:,0]:
        if aa not in x:
            x.append(aa)
    x = np.array(x)
    y = []
    for aa in data[:,1]:
        if aa not in y:
            y.append(aa)
    y = np.array(y)
    D = max(data[:,2])
    dat = data[:,2].reshape(len(x),len(y)).transpose()/D
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.gca().set_aspect('equal','box')
    plt.pcolormesh(x,y,dat,cmap='jet',shading='gouraud')
    plt.show()


def Plot_Fr(filename:str,title:str=None,label:str='r',x_lim:bool=None) -> None:
    """
    Reads and plots data as a function of a single spatial dimention,
    with errorbars (useful e.g. for g(r)).
    """
    data = np.loadtxt(filename,usecols=(0,2,3))
    plt.xlabel(label)
    plt.title(title)
    if x_lim:
        plt.xlim(x_lim)
    plt.errorbar(data[:,0],data[:,1],data[:,2])
    plt.show()


def Main()->None:
    """
    Reads the input file name as an argument or requests it
    from the user, calls function to plot SF fraction, total
    energy, rho(x,y), g(x,y), g(r). Edit accordingly to add/remove
    quantities to be shown.
    """
    if len(sys.argv) == 1:
        fname = input("Please insert file name: ")
    else:
        fname = sys.argv[1]
    PlotTimeSeries('Superfluid_fraction',fname,'Step','Superfluid fraction',(0.,1.25),True)
    PlotTimeSeries('Total_energy',fname,'Step','Total energy')
    Plot_XY('rho_xy.dat','Density')
    Plot_XY('g_xy.dat','Pair correlation function')
    Plot_Fr('gofr.dat','g(r)')


if __name__=="__main__":
    Main()
