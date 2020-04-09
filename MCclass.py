#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 12:47:15 2020

@author: jeremymarck
"""

import numpy as np
import pandas as pd
from math import *
from scipy.stats import norm
import matplotlib.pyplot as plt
from random import gauss, seed


class monte_carlo_pricer:
    
    # Constructeur: 
    def __init__(self, stock, strike, time, r, sigma, nb_simul, M):
        self.stock = stock
        self.strike = strike
        self.time = time
        self.r = r
        self.sigma = sigma
        self.nb_simul = nb_simul
        self.M = M
    
    # Simulating random paths
    def random_paths_euler(self,is_antith):
        dt = self.time/self.M
        res = []
        if is_antith == True:
            for i in range(self.nb_simul):
                path = [self.stock]
                for j in range(1,self.M):
                    z = gauss(0,1)
                    path.append(path[-1]*np.exp((self.r - 0.5*self.sigma**2)*dt +\
                                self.sigma*np.sqrt(dt)*z))
                res.append(path)
        else:
            for i in range(self.nb_simul):
                path = [self.stock]
                for j in range(1,self.M):
                    z = gauss(0,1)
                    path.append(path[-1]*np.exp((self.r - 0.5*self.sigma**2)*dt -\
                                self.sigma*np.sqrt(dt)*z))
                res.append(path)
        return res
    
    # Price of European calls and Puts.
    def european_price(self, is_call, is_antith):
        random_paths = self.random_paths_euler(True)
        res = {}
        l_price = []
        vol = []
        # Choosing the type of the option.
        if is_call == True:
            c = 1
        else:
            c = -1
        # Getting the associated price/delta.
        # Without antithetic.
        if is_antith == False:
            if is_call == True:
                for simul in random_paths:
                    price = np.exp(-self.r*self.time) * max(c*(simul[-1] - self.strike),0)
                    l_price.append(price)
            else:
                for simul in random_paths:
                    price = np.exp(-self.r*self.time) * max(c*(simul[-1] - self.strike),0)
                    l_price.append(price)
            price = round(np.mean(l_price),2)
            res['Price'] = price
            # Getting the confidence interval.
            for simul in random_paths:
                vol.append((simul[-1] - price)**2)
            var = round(float(np.mean(vol)/len(random_paths)),2)
            IC = 1.96 * np.sqrt(var)
            conf_int = [round(price - IC,2), round(price + IC,2)]
            res['Variance'] = var
            res['95% conf.interv.'] = conf_int
            tableau = pd.DataFrame(list(res.items()), columns = ['Parameters', 'Values'])
            print(tableau)
        # With antithetic.
        else:
            if is_call == True:
                random_paths_antith = self.random_paths_euler(False)
                for simul in random_paths:
                    price = np.exp(-self.r*self.time) * max(c*(simul[-1] - self.strike),0)
                    l_price.append(price)
                for simul in random_paths_antith:
                    price = np.exp(-self.r*self.time) * max(c*(simul[-1] - self.strike),0)
                    l_price.append(price)
            else:
                random_paths_antith = self.random_paths_euler(False)
                for simul in random_paths:
                    price = np.exp(-self.r*self.time) * max(c*(simul[-1] - self.strike),0)
                    l_price.append(price)
                for simul in random_paths_antith:
                    price = np.exp(-self.r*self.time) * max(c*(simul[-1] - self.strike),0)
                    l_price.append(price)
            # Computing final results.
            price = round(np.mean(l_price),2)
            res['Price'] = price
            for simul in random_paths:
                vol.append((simul[-1] - price)**2)
            for simul in random_paths_antith:
                vol.append((simul[-1] - price)**2)
            var = round(float(np.mean(vol)/len(2*random_paths)),2)
            IC = 1.96 * np.sqrt(var)
            conf_int = [round(price - IC,2), round(price + IC,2)]
            res['Variance'] = var
            res['95% conf.interv.'] = conf_int
            tableau = pd.DataFrame(list(res.items()), columns = ['Parameters', 'Values'])
            print(tableau)
    
    def digital_option(self, coupon):
        # Generating random paths.
        random_paths = self.random_paths_euler(False)
        # Getting the price.
        l_payoff = []
        vol = []
        res = {}
        for simul in random_paths:
            if simul[-1] >= self.strike:
                price = coupon
                l_payoff.append(price)
            else:
                price = 0
                l_payoff.append(0)
        price = round(np.exp(-self.r*self.time)*np.mean(l_payoff),2)
        res['Price'] = price
        for simul in random_paths:
            vol.append((simul[-1] - price)**2)
        var = round(float(np.mean(vol)/len(random_paths)),2)
        IC = 1.96 * np.sqrt(var)
        conf_int = [round(price - IC,2), round(price + IC,2)]
        res['Variance'] = var
        res['95% conf.interv.'] = conf_int
        tableau = pd.DataFrame(list(res.items()), columns = ['Parameters', 'Values'])
        print(tableau)
       
            
        
            
                       
a = monte_carlo_pricer(100,100,1,0.02,0.25,50000,1000)
#print('Antithétique')
#a.european_price(is_call = True, is_antith = True)
#print('Normal')
#a.european_price(is_call = True, is_antith = False)
#print('Digitale')
a.digital_option(10)

                
                