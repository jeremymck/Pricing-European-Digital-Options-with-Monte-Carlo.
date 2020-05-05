#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 19:51:58 2020

@author: jeremymarck
"""

from math import *
import numpy as np
from scipy.stats import norm
import pandas as pd
import matplotlib.pyplot as plt
from random import gauss, seed

class BlackScholes(object):
    
    # Constructeur.
    def __init__(self, stock, strike, time, rfr, sigma):
        self.stock = stock
        self.strike = strike
        self.time = time
        self.rfr = rfr
        self.sigma = sigma
        self.d1num = (log(self.stock / self.strike) + (self.rfr + 0.5 * self.sigma *self.sigma) * self.time)
        self.d1 = self.d1num / (self.sigma * sqrt(self.time))
        self.d2 = self.d1 - self.sigma * sqrt(self.time)

    # Définitions des fonctions.
    ############################
    
    # Valeur du Call.
    def value_call(self):
        value = - self.stock * norm.cdf(-self.d1) + self.strike * exp(-self.rfr*self.time) * norm.cdf(-self.d2)
        res = value + self.stock - self.strike*np.exp(-self.rfr*self.time)
        return res
    
    # Valeur du Put.
    def value_put(self):
        value = - self.stock * norm.cdf(-self.d1) + self.strike * exp(-self.rfr*self.time) * norm.cdf(-self.d2)
        return value
    
    # Delta du Call.
    def delta_call(self):
        delta = norm.cdf(self.d1)
        return delta

    # Delta du Put.
    def delta_put(self):
        delta = -(norm.cdf(-self.d1))
        return delta
    
    # Gamma Call/Put.
    def gamma(self):
        gamma = norm.pdf(self.d1) / (self.stock * self.sigma * sqrt(self.time))
        return gamma

    # Theta du Put.
    def theta_call(self):
        theta = -(self.stock * norm.pdf(self.d1) * self.sigma / (2 * sqrt(self.time))) - \
        (self.rfr * self.strike * exp(-self.rfr * self.time) * norm.cdf(self.d2))
        return theta
    
    # Theta du Put.
    def theta_put(self):
        theta = -(self.stock * norm.pdf(self.d1) * self.sigma / (2 * sqrt(self.time))) + \
        (self.rfr * self.strike * exp(-self.rfr * self.time) * norm.cdf(-self.d2))
        return theta

    # Rho du Call.
    def rho_call(self):
        rho = (self.time * self.strike * exp(-self.rfr*self.time) * norm.cdf(self.d2))
        return rho
    
    # Rho du Put.
    def rho_put(self):
        rho = (-self.time * self.strike * exp(-self.rfr*self.time) * norm.cdf(-self.d2))
        return rho

    # Vega du Put.
    def vega(self):
        vega = (self.stock * norm.pdf(self.d1) * sqrt(self.time))
        return vega
    
a = BlackScholes(100,100,0.001,0.02,1)
a.value_call()
a.delta_call()
a.gamma()
a.theta_call()
a.vega()
