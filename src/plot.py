"""Plotting module for Fesslix.

"""

# Fesslix - Stochastic Analysis
# Copyright (C) 2010-2025 Wolfgang Betz
#
# Fesslix is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fesslix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Fesslix.  If not, see <http://www.gnu.org/licenses/>. 


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

##################################################
# Constants                                      #
##################################################

color_seq = [ 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan' ]

color_era_1 = (0.35,0.58,0.36)
color_era_2 = (0.36,0.45,0.33)
color_era_3 = (0.35,0.31,0.27)
color_era_grey_1 = '#939598'
color_era_grey_2 = '#636466'
color_era_grey_3 = '#3B3B3C'
color_tumblue = "#0065BD"
color_tumgreen = "#A2AD00"
color_tumorange = "#E37222"
color_tumlightblue = "#98C6EA"
color_tumdarkblue = "#005293"
color_seq_era = [ color_era_1, color_era_2, color_era_3, color_era_grey_1, color_era_grey_2, color_era_grey_3 ]
color_seq_tum = [ color_tumblue, color_tumgreen, color_tumorange, color_tumlightblue, color_tumdarkblue ]

def gen_cmap_gr():
    c = [ "darkgreen", "lawngreen", "yellow", "orange", "red", "darkred"]
    v = [0,.2,.5,0.8,.9,1.]
    l = list(zip(v,c))
    cmap=mpl.colors.LinearSegmentedColormap.from_list('rg',l, N=256)
    return cmap

cmap_gr = gen_cmap_gr()





