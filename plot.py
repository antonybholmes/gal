# -*- coding: utf-8 -*-
"""
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.
This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with 
this program. If not, see <https://www.gnu.org/licenses/>. 

Copyright (C) 2022 Antony Holmes.
"""

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot
import matplotlib.font_manager
import numpy as np
import scipy.stats.kde
import scipy.interpolate
#from scipy.interpolate import spline
from scipy.interpolate import CubicSpline

import sys
import re

PLOT_W = 12
PLOT_H = 8
FONT_PATH = '/ifs/scratch/cancer/Lab_RDF/ngs/Arial.ttf'


def plot_tss(file, color, xlabel, ylabel, prom_ext_5p=2000, prom_ext_3p=1000):
  prop = matplotlib.font_manager.FontProperties(fname=FONT_PATH)
  
  matplotlib.rcParams['font.family'] = prop.get_name()
  matplotlib.rcParams['font.size'] = 20
  #matplotlib.rcParams['font.sans-serif'] = ['Arial']
  
  #import matplotlib.font_manager
  #print matplotlib.font_manager.findSystemFonts(fontpaths=None)
  
  #file = sys.argv[1]
  #id = sys.argv[2]
  #color = sys.argv[3] # '#2c5aa0'
  
  out = re.sub("\.tsv", ".pdf", file)
  out = re.sub("closest", "TSS", out)
  out_smooth = re.sub("\.pdf", "_Sm.pdf", out)
  
  data = np.loadtxt(file, skiprows=1, usecols=range(0,1))
  
  data = data[np.where(np.logical_and(data >= -prom_ext_5p, data <= prom_ext_3p))]
  data = np.divide(data, 1000)
  
  
  #my_pdf = scipy.stats.kde.gaussian_kde(data)
  
  # -5/+4 kb = 90 bins
  hist, bin_edges = np.histogram(data, bins=90)
  
  #bar
  
  fig = matplotlib.pyplot.figure(figsize=[PLOT_W, PLOT_H])
  
  ax = fig.add_subplot(111)
  ax.bar(bin_edges[0:(len(bin_edges) - 1)], hist, width=0.1, facecolor=color, edgecolor='none')
  
  matplotlib.pyplot.xlabel(xlabel)
  matplotlib.pyplot.ylabel(ylabel)
  
  ymax = 100 * (np.ceil(np.max(hist) / 100.0))
  
  ax.set_xlim([-5, 4])
  ax.set_xticks(range(-5, 5))
  ax.set_ylim([0, ymax])
  
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')
  
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(out)
  
  #smooth
  
  
  x = bin_edges[0:(len(bin_edges) - 1)]
  y = hist

  cs = CubicSpline(x, y)

  x_smooth = np.arange(-5, 4, 0.01) #np.linspace(x_sm.min(), x_sm.max(), 200)
  y_smooth = cs(x_smooth) #spline(x, y, x_smooth)
  
  fig = matplotlib.pyplot.figure(figsize=[PLOT_W, PLOT_H])
  
  ax = fig.add_subplot(111)
  
  # curves using kernel density
  #my_pdf = scipy.stats.kde.gaussian_kde(data)
  #xr = np.arange(-5, 4, 0.01)
  #yr = my_pdf(xr)
  #yr = np.divide(yr, np.max(yr))
  #yr = np.multiply(yr, np.max(hist))
  
  ax.plot(x_smooth, y_smooth, linewidth=1, color=color)
  ax.fill_between(x_smooth, y_smooth, edgecolor="none", facecolor=color, alpha=0.8)
  
  matplotlib.pyplot.xlabel(xlabel)
  matplotlib.pyplot.ylabel(ylabel)
  
  ax.set_xlim([-5, 4])
  ax.set_xticks(range(-5, 5))
  ax.set_ylim([0, ymax])
  
  # remove border
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')
  
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(out_smooth)


def plot_log10_tss(file, color, xlabel, ylabel):
  prop = matplotlib.font_manager.FontProperties(fname=FONT_PATH)
  
  matplotlib.rcParams['font.family'] = prop.get_name()
  matplotlib.rcParams['font.size'] = 20
  #matplotlib.rcParams['font.sans-serif'] = ['Arial']
  
  #import matplotlib.font_manager
  #print matplotlib.font_manager.findSystemFonts(fontpaths=None)
  
  #file = sys.argv[1]
  #id = sys.argv[2]
  #color = sys.argv[3] #'#2c5aa0'
  
  out = re.sub("\.tsv", "_Log10_Abs_Dist.pdf", file)
  out_smooth = re.sub("\.tsv", "_Log10_Abs_Dist_Sm.pdf", file)
  
  data = np.loadtxt(file, skiprows=1, usecols=range(0, 1)) 
  data = data[np.where(data != 0)]
  data = np.abs(data)
  data[np.where(data < 1)] = 1
  data = np.log10(data)
  #data[np.where(data < 0)] = 0
  data = data[np.where(data <= 7)]
  
  labels = []
  
  for i in range(0, 8):
    labels.append("1E" + str(i)) #"{:.0E}".format(10 ** i))
  

  # log bins up to 8
  
  hist, bin_edges = np.histogram(data, bins=70)
  
  # calculate percentage of peaks less than 1000 and more than
  # 1000 from tss
  
  tc = 0.0
  
  for h in hist:
    tc += h
    
  tl = 0.0

  for i in range(0, len(bin_edges) - 1):
    if bin_edges[i] <= 3:
      tl += hist[i]
  
  tr = 0.0
  
  for i in range(0, len(bin_edges) - 1):
    if bin_edges[i] > 3:
      tr += hist[i]
  
  #sys.stderr.write("tc\t" + str(tc) + "\n")
  sys.stderr.write("<= 1000\t" + str(tl / tc) + "\n")
  sys.stderr.write("> 1000\t" + str(tr / tc) + "\n")
  
  # bar
  
  fig = matplotlib.pyplot.figure(figsize=[PLOT_W, PLOT_H])
  
  ax = fig.add_subplot(111)
  
  ax.bar(bin_edges[0:(len(bin_edges) - 1)], hist, width=0.1, facecolor=color, edgecolor='none')
  
  matplotlib.pyplot.xlabel(xlabel)
  matplotlib.pyplot.ylabel(ylabel)
  
  ax.set_xlim([0, 7])
  
  ymax = 100 * (np.ceil(np.max(hist) / 100.0))
  
  ax.set_xticks(range(0, 8))
  ax.set_xticklabels(labels)
  ax.set_ylim([0, ymax])
  
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')
  
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(out)
  
  #
  # Smoothed using splines
  #
  
  x = bin_edges[0:(len(bin_edges) - 1)]
  y = hist
  
  sys.stderr.write("has " + str(len(x)) + " " + str(len(y)) + "\n")

  cs = CubicSpline(x, y)

  x_smooth = np.arange(0, 7, 0.01) #np.linspace(x_sm.min(), x_sm.max(), 200)
  y_smooth = cs(x_smooth) #spline(x, y, x_smooth)
  
  
  #my_pdf = scipy.stats.kde.gaussian_kde(data)
  #xr = np.arange(0, 8, 0.01)
  #yr = my_pdf(xr)
  #yr = np.divide(yr, np.max(yr))
  #yr = np.multiply(yr, np.max(hist))
  
  fig = matplotlib.pyplot.figure(figsize=[PLOT_W, PLOT_H])
  ax = fig.add_subplot(111)
  ax.plot(x_smooth, y_smooth, linewidth=1, color=color)
  ax.fill_between(x_smooth, y_smooth, edgecolor="none", facecolor=color, alpha=0.8)
  
  matplotlib.pyplot.xlabel(xlabel)
  matplotlib.pyplot.ylabel(ylabel)
  
  ax.set_xlim([0, 7])
  
  ax.set_xticks(range(0, 8))
  ax.set_xticklabels(labels)
  ax.set_ylim([0, ymax])
  
  # remove border
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')
  
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(out_smooth)
  
