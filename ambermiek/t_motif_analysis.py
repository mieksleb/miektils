#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 15:18:47 2022

analysis and plot of t_motif data

@author: michaelselby
"""
import numpy as np
import matplotlib.pyplot as plt


polarities = [3]
bulges = [3,4,5]

list_of_dics = []

for polarity in polarities:
    for bulge in bulges:
        direc = "/Users/michaelselby/OneDrive - Nexus365/DPhil/Behnam/T_motifs_michael/"
        alphas, betas, phis = [], [], []
        direc += str(polarity) + "E/" + str(bulge) + "/"
        with open( direc + "angles.txt", "r") as file:
            for line in file:
                alpha, beta, phi = [float(x) for x in line.split()[:3]]
                alphas.append(alpha)
                betas.append(beta)
                phis.append(phi)
                
        dic = {}
        dic["polarity"] = polarity
        dic["bulge"] = bulge
        dic["av_alpha"], dic["std_alpha"] = np.mean(alphas), np.std(alphas)
        dic["av_beta"], dic["std_beta"] = np.mean(betas), np.std(betas)
        dic["av_phi"], dic["std_phi"] = np.mean(phis), np.std(phis)
        list_of_dics.append(dic)

        
                
cap = 3
fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
for ax in [ax1,ax2]:
    ax.set_xticks([3,4,5,6,7,8,9])
    ax.set_yticks([-90,0,90,180])
    ax.set_xlim([2.5,9.5])
    ax.set_ylim([-90,180])
    ax.set_xticklabels([3,4,5,6,7,8,9])
    ax.set_yticklabels([-90,0,90,180])
    ax.set_xlabel("L/nt")

for dic in list_of_dics:
    if dic["polarity"] == 3:
        ax = ax2
        # ax.set_title(r'\bold{T}_{5\'}^L\left[L,\bar{s}\right]')
    else:
        ax = ax1
        ax.set_ylabel("Angle")

    ax.errorbar(dic["bulge"],dic["av_alpha"],np.linspace(-dic["std_alpha"],dic["std_alpha"],1), color="green", marker="o", capsize=cap)
    ax.errorbar(dic["bulge"],dic["av_beta"],np.linspace(-dic["std_beta"],dic["std_beta"],1), color="crimson", marker="s", capsize=cap)
    ax.errorbar(dic["bulge"],dic["av_phi"],np.linspace(-dic["std_phi"],dic["std_phi"],1), color="blue", marker="v",capsize=cap)
    
plt.show()
    
                
        