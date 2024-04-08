# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 22:03:47 2022

@author: Igor
"""
import matplotlib.pyplot as plt
from obspy import read_inventory

input_dir = 'H:\Central_Kamchatka\Meta_data'

inv1_path = (input_dir + '\inventory_perm_stations.xml')
inv1 = read_inventory(inv1_path)
inv2_path = (input_dir + '\inventory_temp_stations.xml')
inv2 = read_inventory(inv2_path)


fig = inv1.plot(projection='local', resolution='h', marker='^', color='r',
                method='basemap', size=80, label=True, )
inv2.plot(projection='local', resolution='h', marker='^', color='b',
          method='basemap', size=80, label=True, fig=fig)

# legend
plt.scatter([], [], s=80, c='r', marker='^', label='permanent')
plt.scatter([], [], s=80, c='b', marker='^', label='temporary')
plt.legend(loc='upper right', frameon=True, framealpha=1, scatterpoints=1,
           labelspacing=0.5, title='Station type:', title_fontsize=13,
           fontsize=13, )
