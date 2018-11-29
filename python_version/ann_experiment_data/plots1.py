#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt


if __name__ == "__main__":
    ds = pd.read_csv("more_ann_results.txt", header=0, sep=';', dtype = { 'Epsilon' : object })
    # ds = pd.read_csv("ann_results_single.txt", header=0, sep=';', dtype = { 'Epsilon' : object })
    index_cols = ['Dimension', 'Epsilon',  'point_gen_method', 'query_method', 'Npoints']

    ds.drop(columns=['MinFraction', 'MaxFraction'], inplace=True)

    print(ds.describe())
    print("********************************************************************************")
    print(ds.head(20))
    print("********************************************************************************")
    ds.set_index(index_cols, inplace=True)

    ds_to_plot = ds.loc[(6, '0.01', 'get_uniform_points', 'query_2')]
    ds_to_plot.reset_index(inplace=True)
    x = ds_to_plot.Npoints
    y6_01 = ds_to_plot.AvgNumberOfDistances

    ds_to_plot = ds.loc[(6, '0.005', 'get_uniform_points', 'query_2')]
    ds_to_plot.reset_index(inplace=True)
    y6_005 = ds_to_plot.AvgNumberOfDistances

    ds_to_plot = ds.loc[(6, '0.001', 'get_uniform_points', 'query_2')]
    ds_to_plot.reset_index(inplace=True)
    y6_001 = ds_to_plot.AvgNumberOfDistances

    ds_to_plot = ds.loc[(6, '0.0005', 'get_uniform_points', 'query_2')]
    ds_to_plot.reset_index(inplace=True)
    y6_0005 = ds_to_plot.AvgNumberOfDistances

    ds_to_plot = ds.loc[(6, '0.0001', 'get_uniform_points', 'query_2')]
    ds_to_plot.reset_index(inplace=True)
    y6_0001 = ds_to_plot.AvgNumberOfDistances

    # fit = np.polyfit(x, y, 1)
    # fit_fn = np.poly1d(fit)
    # plt.plot(x, y, 'x', x, fit_fn(x), '--r')

    plt.plot(x, y6_01, x, y6_005, x, y6_001, x, y6_0005, x, y6_0001)
    plt.xlabel('#points')
    plt.ylabel('#computed distances')

    # ds_to_plot = ds_to_plot.set_index(["Npoints"])
    # ds_to_plot = ds_to_plot.sort_index(ascending=True)
    # # print(ds_to_plot.index)
    # ds_to_plot.plot()

    plt.show()
