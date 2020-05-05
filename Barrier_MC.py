#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 11:14:24 2020

@author: jeremymarck
"""

import numpy as np
import numpy.random as npr
import math
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import xlwings as xw
from random import gauss

class barrier_options_pricer_MC:
    
    # Constructeur :
    def __init__(self, stock, barrier, strike, time, r, sigma, Type, Option, nb_simul, M):
        self.stock = stock  #Spot price of the underlying.
        self.barrier = barrier  #Barrier level.
        self.strike = strike  #Strike level.
        self.time = time  #Time to maturity.
        self.r = r  #Free risk interest rate.
        self.sigma = sigma  #Volatility.
        self.Type = Type  #Option type (Call or put)
        self.Option = Option  #Knock out or Knck in option
        self.nb_simul = nb_simul  #Number of simulations used in the Monte Carlo process.
        self.M = M  #Time step. 
    
    def generate_paths(self):
        dt = self.time/self.M #Time step setup.
        # Black-Scholes dynamics for the underlying.
        # dS_t/S_t = r*dt + sigma*dW_t (risk neutral world)
        # Discretization using basic Euler Scheme: 
        # S_{t+dt} = S_t[1+ r*dt + sigma*dW_t]
        paths = []
        for i in range(self.nb_simul):
            path = [self.stock]
            for j in range(1,self.M):
                z = gauss(0,1)
                #S = path[-1]*(1 + self.r*dt + self.sigma*math.sqrt(dt)*z)
                S = path[-1] + path[-1]*self.r*dt + path[-1]*self.sigma*dt*z
                path.append(S)
            paths.append(path)
        return paths
    
    def options_value(self):
        
        # Generating paths using the previous function?
        paths = self.generate_paths()
        
        ###################################
        # Options à barrière desactivante #
        ###################################
        
        # Up and Out Call.
        ##################
        # Mechanism: if any point in the trajectory reaches the barrier,
        # the value of the path is set to zero.
        if self.Type == 'Call' and self.Option == 'Out':
            for path in paths:
                maximum = max(path)
                if maximum > self.barrier:
                    n = len(path)
                    for i in range(n):
                        path[i] = 0
        
        # Down and Out Call.
        ####################
        # Mechanism: if any point in the trajectory reaches the barrier,
        # the value of the path is set to zero.
        if self.Type == 'Put' and self.Option == 'Out':
            for path in paths:
                minimum = min(path)
                if min < self.barrier:
                    n = len(path)
                    for i in range(n):
                        path[i] = 0
                        
        ################################
        # Options à barrière activante #
        ################################
        
        # Up and In Call.
        #################
        # Mechanism: if any point in the trajectory doesn't reach the barrier,
        # the value of the path is set to zero.
        if self.Type == 'Call' and self.Option == 'In':
            for path in paths:
                maximum = max(path)
                if max < self.barrier:
                    n = len(path)
                    for i in range(n):
                        path[i] = 0
        
        # Down and In Call.
        ###################
        # Mechanism: if any point in the trajectory doesn't reach the barrier,
        # the value of the path is set to zero.
        if self.Type == 'Put' and self.Option == 'In':
            for path in paths:
                minimum = min(path)
                if minimum > self.barrier:
                    n = len(path)
                    for i in range(n):
                        path[i] = 0
        
        ##############
        # Evaluation #
        ##############
        if self.Type == 'Call':
            payoff = 0
            for path in paths:
                payoff += max(path[-1] - self.strike, 0)
            payoff = payoff/(self.nb_simul)
            option_value = np.exp(-self.r * self.time)*payoff
            return option_value
        
        if self.Type == 'Put':
            payoff = 0
            for path in paths:
                payoff += max(self.strike - path[-1], 0)
            payoff = payoff/(self.nb_simul)
            option_value = np.exp(-self.r * self.time)*payoff
            return option_value
                
    
a = barrier_options_pricer_MC(100,130,100,1,0.02,1,'Call','Out',100000,252)
b = a.generate_paths()
c = a.options_value()
print('Price',c)
