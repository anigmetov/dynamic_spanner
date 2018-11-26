import pandas as pd
import sys
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # ds = pd.read_csv("ann_results.txt", header=0, sep=';', dtype = { 'Epsilon' : object })
    ds = pd.read_csv("ann_results.txt", header=0, sep=';', dtype = { 'Epsilon' : object })
    index_cols = ['Dimension', 'Epsilon',  'point_gen_method', 'query_method', 'Npoints']

    ds.drop(columns=['MinFraction', 'MaxFraction', 'MinNumberOfDistances', 'MaxNumberOfDistances'], inplace=True)

    print(ds.describe())
    print("********************************************************************************")
    print(ds.head(20))
    print("********************************************************************************")
    ds.set_index(index_cols, inplace=True)
    # sys.exit(0)


    # ds_non_blind_greedy.drop(columns=['spanner_method'], inplace=True)
    # ds_blind_greedy.drop(columns=['spanner_method'], inplace=True)
    # ds_blind_random.drop(columns=['spanner_method'], inplace=True)
    # ds_quasi_greedy.drop(columns=['spanner_method'], inplace=True)
    # ds_quasi_shaker.drop(columns=['spanner_method'], inplace=True)
    # ds_blind_random_bad_ratio.drop(columns=['spanner_method'], inplace=True)
    # ds_blind_random_bad_ratio_lbf.drop(columns=['spanner_method'], inplace=True)
    # ds_blind_random_bad_ratio_cf.drop(columns=['spanner_method'], inplace=True)
    # ds_blind_random_bad_ratio_cf_lbf.drop(columns=['spanner_method'], inplace=True)
    #
    # ds_all_dims = ds_non_blind_greedy.join(ds_blind_greedy, lsuffix=' (non-blind greedy)', rsuffix=' (blind greedy)', how="outer")
    # ds_all_dims = ds_all_dims.join(ds_blind_random, rsuffix=' (blind random)')
    # ds_all_dims = ds_all_dims.join(ds_quasi_greedy, rsuffix=' (quasi-sorted greedy').join(ds_quasi_shaker, rsuffix=' (quasi-sorted shaker)')
    # ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio, rsuffix=' (blind random bad ratio)')
    # ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio_lbf, rsuffix=' (blind random bad ratio, force lower bound)')
    # ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio_cf, rsuffix=' (blind random bad ratio, force connectedness)')
    # ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio_cf_lbf, rsuffix=' (blind random bad ratio, force connectedness and lower bound)')

    # ycols = ["sparseness_nbg", "sparseness_bg", "sparseness_rbr_lbf"]
    # ycols = ["n_edges_nbg", "n_edges_bg", "n_edges_rbr_lbf", "n_edges_qs", "n_edges_qg"]
    # ycols = ["edges_to_points_ratio_nbg", "edges_to_points_ratio_bg", "edges_to_points_ratio_rbr_lbf",
    #          "edges_to_points_ratio_qs", "edges_to_points_ratio_qg"]

    # ycols = ["Edges (non-blind greedy)", "Edges (blind greedy)", "Edges (blind random bad ratio, force lower bound)", "Edges (blind random bad ratio)"]
    # ycols = ["Edges to points ratio (non-blind greedy)", "Edges to points ratio (blind greedy)", "Edges to points ratio (blind random bad ratio, force lower bound)", "Edges to points ratio (blind random bad ratio)"]
    # # # ycols = ["sparseness_random", "sparseness_greedy", "sparseness_blind_greedy"]
    ds_to_plot = ds.loc[(2, '0.01', 'normal', 'query_1')]
    ds_to_plot.reset_index(inplace=True)
    ds_to_plot = ds_to_plot.set_index(["Npoints"])
    ds_to_plot = ds_to_plot.sort_index(ascending=True)
    # print(ds_to_plot.index)
    ds_to_plot.plot()
    plt.show()
