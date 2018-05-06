#!/usr/bin/python26

from __future__ import print_function
from getdist import plots, MCSamples, loadMCSamples
import getdist, IPython

A611=loadMCSamples('chains/examp/A611')
new_chains = raw_input('Path to new chains:').strip()
A6112=loadMCSamples(new_chains)

g=plots.getSubplotPlotter(width_inch=8)
g.settings.axes_fontsize=8
g.settings.alpha_filled_add=0.4
g.triangle_plot([A611, A6112], ['1', '2', '3', '4', '5', '6'],
    filled_compare=True, 
    legend_labels=['Example', 'New chains'], 
    legend_loc='upper right', 
    line_args=[{'ls':'--', 'color':'green'},
               {'lw':2, 'color':'darkblue'}], 
    contour_colors=['green','darkblue'])
g.export('A611_tri.png')

