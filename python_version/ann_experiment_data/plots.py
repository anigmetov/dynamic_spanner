#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

if __name__ == "__main__":
    if True:
        ds = pd.read_csv("averaged_ann_data.csv", header=0, sep=';', dtype={'epsilon': object})
        # ds.drop(columns=["epsilon", "point_gen_method", "query_method"], inplace=True)
        index_cols = ['dim']
        ds.set_index(index_cols, inplace=True)

        ds_to_plot2 = ds.loc[(2)]
        ds_to_plot5 = ds.loc[(5)]
        ds_to_plot10 = ds.loc[(10)]
        # ds_to_plot20 = ds.loc[(20)]

        ds_to_plot2.set_index(["Points"], inplace=True)
        ds_to_plot5.set_index(["Points"], inplace=True)
        ds_to_plot10.set_index(["Points"], inplace=True)


        ds_to_plot2.rename(index=str, columns={"Computed distances": "Computed distances (dim 2)"}, inplace=True)
        ds_to_plot5.rename(index=str, columns={"Computed distances": "Computed distances (dim 5)"}, inplace=True)
        ds_to_plot10.rename(index=str, columns={"Computed distances": "Computed distances (dim 10)"}, inplace=True)
        # ds_to_plot20.set_index(["Points"], inplace=True)

        ds_all_dims = ds_to_plot2.join(ds_to_plot5)
        ds_all_dims = ds_all_dims.join(ds_to_plot10)
        # ds_all_dims = ds_all_dims.join(ds_to_plot20, rsuffix=" (dim 20)")

        print(ds_all_dims.index.values)
        ds_all_dims.plot(use_index = True)

        plt.show()
        sys.exit(0)

    if False:
        ds = pd.read_csv("averaged_ann_data.csv", header=0, sep=';', dtype={'epsilon': object})

        # index_cols = ['dim', 'epsilon',  'point_gen_method', 'query_method']
        # ds.set_index(index_cols, inplace=True)
        # ds_to_plot2 = ds.loc[(2, '0.01', 'get_uniform_points', 'query_1')]
        # ds_to_plot10 = ds.loc[(10, '0.01', 'get_uniform_points', 'query_1')]
        # ds_to_plot20 = ds.loc[(20, '0.01', 'get_uniform_points', 'query_1')]

        index_cols = ['dim']
        ds.drop(["epsilon", "point_gen_method", "query_method"], inplace=True)
        ds.set_index(index_cols, inplace=True)

        ds_to_plot2 = ds.loc[(2)]
        ds_to_plot10 = ds.loc[(10)]
        ds_to_plot20 = ds.loc[(20)]

        ds_to_plot2.set_index(["npoints"], inplace=True)

        x = ds_to_plot2.n_points
        y = ds_to_plot2.dist_num

        y = y / np.log(x)
        # x = np.log(x)
        print(x)
        print(y)
        fit = np.polyfit(x, y, 1)
        fit_fn = np.poly1d(fit)
        print(fit_fn(x))
        plt.xlabel("#points")
        plt.ylabel("#computed distances / log(#points)")
        plt.plot(x, y)
        # plt.plot(x, y, 'ro', x, fit_fn(x), '--k')
        plt.show()
        sys.exit(0)

    ds = pd.read_csv("more_ann_results.txt", header=0, sep=';', dtype={'Epsilon': object})
    # ds = pd.read_csv("ann_results_single.txt", header=0, sep=';', dtype = { 'Epsilon' : object })
    index_cols = ['Dimension', 'Epsilon', 'point_gen_method', 'query_method', 'Npoints']

    ds.drop(columns=['MinFraction', 'MaxFraction'], inplace=True)

    print(ds.describe())
    print("********************************************************************************")
    print(ds.head(20))
    print("********************************************************************************")
    ds.set_index(index_cols, inplace=True)

    ds_to_plot = ds.loc[(2, '0.01', 'get_uniform_points', 'query_1')]
    ds_to_plot.reset_index(inplace=True)
    x = ds_to_plot.Npoints
    y2_01 = ds_to_plot.MaxNumberOfDistances

    ds_to_plot = ds.loc[(4, '0.01', 'get_uniform_points', 'query_1')]
    ds_to_plot.reset_index(inplace=True)
    x = ds_to_plot.Npoints
    y4_01 = ds_to_plot.MaxNumberOfDistances

    ds_to_plot = ds.loc[(6, '0.01', 'get_uniform_points', 'query_1')]
    ds_to_plot.reset_index(inplace=True)
    x = ds_to_plot.Npoints
    y6_01 = ds_to_plot.MaxNumberOfDistances

    # fit = np.polyfit(x, y, 1)
    # fit_fn = np.poly1d(fit)
    # plt.plot(x, y, 'x', x, fit_fn(x), '--r')

    plt.plot(x, y2_01, x, y4_01, x, y6_01)
    plt.xlabel('#points')
    plt.ylabel('#computed distances log(#points)')

    # ds_to_plot = ds_to_plot.set_index(["Npoints"])
    # ds_to_plot = ds_to_plot.sort_index(ascending=True)
    # # print(ds_to_plot.index)
    # ds_to_plot.plot()

    plt.show()
