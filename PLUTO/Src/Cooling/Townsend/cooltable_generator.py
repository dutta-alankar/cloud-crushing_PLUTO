'''
Last modified  : 2021/06/18, 17:18:50
@author        : alankar dutta, ritali ghosh
Description    : generating cooltable_townsend.dat. It calculates the Townsend Y function
(contains stacked value of T, LAMBDA(T), Y, Flipped value of Y and Flipped value of T)
'''
 
import numpy as np
from scipy.interpolate import interp1d
#from scipy.integrate import quad
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from numba import jit

LAMBDA = np.loadtxt('cooltable.dat')
Tstart, Tstop = np.log10(np.min(LAMBDA[:,0])+0.5), np.log10(np.max(LAMBDA[:,0])-0.5)
LAMBDA = interp1d(LAMBDA[:,0], LAMBDA[:,1])
T_ref = 1.e9 #Reference Temperature

#@jit()
def townsendY(temperature):
    dlogT = 0.0001 
    temp = np.logspace(np.log10(temperature), np.log10(T_ref), 
        int(np.abs(np.log10(temperature)-np.log10(T_ref))/dlogT))
    if (temp.shape[0]==0): return 0.
    Y = cumtrapz(LAMBDA(T_ref)/LAMBDA(temp), temp, initial=0)[-1]
    #Y = LAMBDA(T_ref)*quad(lambda T :1/LAMBDA(T), temperature, T_ref, limit=100)[0] //Quad causes error in calculation of Y
    return Y/T_ref
    
temperature = np.logspace(Tstart,Tstop,1000)
Y = np.array([townsendY(T) for T in temperature])
data = np.vstack((temperature,LAMBDA(temperature),Y, np.flip(Y), np.flip(temperature))).T
np.savetxt('cooltable_townsend.dat',data)
