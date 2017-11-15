'''
Created on Feb 5, 2015

@author: wange
'''
import pylab
import numpy as np
from pylab import show
from scipy.constants import *
import matplotlib.pyplot as plt
from scipy.stats import maxwell
import math
from pprint import pprint
from numpy import genfromtxt, hstack
from scipy.interpolate import interp1d
import inspect
import os


def main():
    a0=3
    xshift=5
    def fit_func(x,a,t1,t2):
        print(a0,xshift)
        return a*np.exp(-(x-xshift)/t1)+(a0-a)*np.exp(-(x-xshift)/t2)
    print(fit_func(1,1,1,1))
if __name__ == '__main__':
    main()