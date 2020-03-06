#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cmocean as cmo

red_colormap = color_map=sns.blend_palette(["0.9", sns.xkcd_rgb["bright red"]], as_cmap=True)
tab20a = sns.mpl_palette("tab20_r", 20)
tab20b = sns.mpl_palette("tab20b_r", 20)
tab20c = sns.mpl_palette("tab20c_r", 20)
main_palette = (
    tab20c[17:] + [tab20a[1]] + tab20c[6:8] + tab20b[17:] +
    tab20b[5:7] + [tab20a[8]] + [tab20b[2]] + [tab20b[1]] + tab20b[9:11] +
    tab20c[9:12] + [tab20c[14]]
)
nonneuronal_palette = (
    tab20c[8:12] + tab20b[14:16] + tab20b[4:7] + tab20c[18:] + tab20b[9:11] + tab20b[1:3] + 
    tab20b[16:19:2] + [tab20a[8]]
)
neuronal_palette = (
    tab20c[13:16:2] + tab20c[9:12] + tab20b[1:3] + tab20b[4:7] + tab20c[16:] + [tab20a[1]] +
    [tab20b[9]] + [tab20a[5]] + tab20b[13:15] + [tab20a[9]]
)
neuronal_nonneuronal_palette = sns.xkcd_palette(["cerulean", "kelly green"])

heatmap_cmap = cmo.tools.cmap(cmo.cm.balance(np.linspace(0,1,256), 0.9))
