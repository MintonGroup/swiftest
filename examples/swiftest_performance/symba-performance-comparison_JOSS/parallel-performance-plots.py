#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
swomp = {}
swest = {} 
swift = {}
npl =  [1,  2,  4,  8,  16, 32]
ymax = [24, 24, 24, 24, 24, 24]

for n in npl:
    swomp[n] = pd.read_csv(f"swifter-omp/pl{n:02d}k/swifter-{n}k-timehist.log")
    swest[n] = pd.read_csv(f"swiftest/pl{n:02d}k/swiftest-{n}k-timehist.log")
    swift[n] = pd.read_csv(f"swift/pl{n:02d}k/swift-{n}k-timehist.log")


CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
y1style =(0, (5, 1)) # densly dashed line

x1line = np.arange(1,24)
y1line = np.full_like(x1line, 1.0)

for i,n in enumerate(npl):
    nsteps = 1000
    axes_fontsize = 32
    legend_fontsize = 24
    fig = plt.figure(1, figsize=(14,10), facecolor="white")
    ax = fig.add_subplot(111)
    plt.setp(ax.get_xticklabels(), fontsize=axes_fontsize)
    plt.setp(ax.get_yticklabels(), fontsize=axes_fontsize)
    ax.set_ylim([0,ymax[i]])
    ax.set_xlabel("Number of CPUs", fontsize=axes_fontsize)
    ax.set_ylabel("Speedup (relative to Swift-SyMBA)", fontsize=axes_fontsize)
    ax.set_facecolor("white")
    plt.plot(swest[n]['N cores'], swift[n]['wall time(s)'][0] / swest[n]['wall time(s)'], alpha=0.5, linewidth=6, label="Swiftest")  
    plt.plot(swomp[n]['N cores'], swift[n]['wall time(s)'][0] / swomp[n]['wall time(s)'], c="darkgreen", alpha=0.5, linewidth=6, label="Swifter-OMP")
    plt.plot(x1line, y1line, c='k', linestyle=y1style, linewidth = 3, label ="Swift")
    idealx = swest[n]['N cores']
    idealy = swest[n]['N cores'] * swift[n]['wall time(s)'][0] / swest[n]['wall time(s)'][0] 
    plt.plot(idealx, idealy, c='k', linestyle=':', linewidth = 3, label ="ideal")
    if n < 0:
        plt.legend(loc='upper right', fontsize=legend_fontsize, markerscale=20)
    else:
        plt.legend(loc='upper left', fontsize=legend_fontsize, markerscale=20)
    plt.title(f"{n}k fully interacting bodies with SyMBA", fontsize=axes_fontsize)
    fig.savefig(f"swiftest_vs_swifter_{n}k.png")
    plt.close(fig)
