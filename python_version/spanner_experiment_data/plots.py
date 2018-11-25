import pandas as pd
import sys
import matplotlib.pyplot as plt

# ds = pd.read_csv("all_data.txt", header=0, sep=';')
#
# index_cols = ['dim', 'n_points', 'point_method', 'epsilon']
# ds_greedy = ds.loc[(ds['spanner_method'] == 'greedy')].set_index(index_cols, drop=False)
# ds_random = ds.loc[(ds['spanner_method'] == 'random')].set_index(index_cols, drop=False)
# ds_quasi = ds.loc[(ds['spanner_method'] == 'quasi')].set_index(index_cols, drop=False)
# ds_blind_greedy = ds.loc[(ds['spanner_method'] == 'blind_greedy')].set_index(index_cols, drop=False)
#
# ds_greedy.drop(columns=['spanner_method'], inplace=True)
# ds_random.drop(columns=['spanner_method'], inplace=True)
# ds_quasi.drop(columns=['spanner_method'], inplace=True)
# ds_blind_greedy.drop(columns=['spanner_method'], inplace=True)
#
# ds_all_dims = ds_greedy.join(ds_random, lsuffix='_greedy', rsuffix='_random').join(ds_quasi, rsuffix='_quasi').join(
#     ds_blind_greedy, rsuffix='_blind_greedy')
#
# print(ds_all_dims.head())
# print(ds_all_dims.describe(include='all').to_string())
#
# ycols = ["sparseness_random", "sparseness_greedy", "sparseness", "sparseness_blind_greedy"]
# # ycols = ["sparseness_random", "sparseness_greedy", "sparseness_blind_greedy"]
# mask_1 = (ds_all_dims['dim'] == 2) & (ds_all_dims['point_method'] == 'get_exponential_points') & (ds_all_dims['epsilon'] == 0.01)
# ds_all_dims.loc[mask_1].plot(x="n_points", y=ycols)
# # mask_1 = (ds_all_dims['n_points'] == 200) & (ds_all_dims['point_method'] == 'get_exponential_points') & (ds_all_dims['epsilon'] == 0.01)
# # ds_all_dims.loc[mask_1].plot(x="dim", y=ycols)
# plt.show()


if __name__ == "__main__":
    ds = pd.read_csv("spanner_cpp_logs.txt", header=0, sep=';')
    index_cols = ['dim', 'point_gen_method', 'epsilon', 'input_file', 'n_points']

    ds_non_blind_greedy = ds.loc[(ds['spanner_method'] == 'greedy-non-blind')].set_index(index_cols, drop=True)
    ds_blind_greedy = ds.loc[(ds['spanner_method'] == 'blind-greedy')].set_index(index_cols, drop=True)
    ds_blind_random = ds.loc[(ds['spanner_method'] == 'blind-random')].set_index(index_cols, drop=True)
    ds_quasi_greedy = ds.loc[(ds['spanner_method'] == 'blind-quasi-sorted-greedy')].set_index(index_cols, drop=True)
    ds_quasi_shaker = ds.loc[(ds['spanner_method'] == 'blind-quasi-sorted-shaker')].set_index(index_cols, drop=True)
    ds_blind_random_bad_ratio = ds.loc[(ds['spanner_method'] == 'blind-random-bad-ratio')].set_index(index_cols,
                                                                                                     drop=True)
    ds_blind_random_bad_ratio_lbf = ds.loc[
        (ds['spanner_method'] == 'blind-random-bad-ratio-lower-bound-first')].set_index(index_cols, drop=True)
    ds_blind_random_bad_ratio_cf = ds.loc[(ds['spanner_method'] == 'blind-random-bad-ratio-connect-first')].set_index(
        index_cols, drop=True)
    ds_blind_random_bad_ratio_cf_lbf = ds.loc[
        (ds['spanner_method'] == 'blind-random-bad-ratio-connect-first-lower-bound-first')].set_index(index_cols,
                                                                                                      drop=True)

    print(ds_non_blind_greedy.size)
    # print(ds_quasi_shaker.head(20))
    print("********************************************************************************")
    # sys.exit(0)

    ds_non_blind_greedy.drop(columns=['spanner_method'], inplace=True)
    ds_blind_greedy.drop(columns=['spanner_method'], inplace=True)
    ds_blind_random.drop(columns=['spanner_method'], inplace=True)
    ds_quasi_greedy.drop(columns=['spanner_method'], inplace=True)
    ds_quasi_shaker.drop(columns=['spanner_method'], inplace=True)
    ds_blind_random_bad_ratio.drop(columns=['spanner_method'], inplace=True)
    ds_blind_random_bad_ratio_lbf.drop(columns=['spanner_method'], inplace=True)
    ds_blind_random_bad_ratio_cf.drop(columns=['spanner_method'], inplace=True)
    ds_blind_random_bad_ratio_cf_lbf.drop(columns=['spanner_method'], inplace=True)

    ds_all_dims = ds_non_blind_greedy.join(ds_blind_greedy, lsuffix='_nbg', rsuffix='_bg', how="outer")
    ds_all_dims = ds_all_dims.join(ds_blind_random, rsuffix='_br')
    ds_all_dims = ds_all_dims.join(ds_quasi_greedy, rsuffix='_qg').join(ds_quasi_shaker, rsuffix='_qs')
    ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio, rsuffix='_rbr')
    ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio_lbf, rsuffix='_rbr_lbf')
    ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio_cf, rsuffix='_rbr_cf')
    ds_all_dims = ds_all_dims.join(ds_blind_random_bad_ratio_cf_lbf, rsuffix='_rbr_cf_lbf')

    # ycols = ["sparseness_nbg", "sparseness_bg", "sparseness_rbr_lbf"]
    # ycols = ["n_edges_nbg", "n_edges_bg", "n_edges_rbr_lbf", "n_edges_qs", "n_edges_qg"]
    ycols = ["edges_to_points_ratio_nbg", "edges_to_points_ratio_bg", "edges_to_points_ratio_rbr_lbf",
             "edges_to_points_ratio_qs", "edges_to_points_ratio_qg"]

    ycols = ["edges_to_points_ratio_nbg", "edges_to_points_ratio_bg", "edges_to_points_ratio_rbr_lbf" ]
     # ycols = ["n_edges_nbg", "n_edges_bg", "n_edges_rbr_lbf", "n_edges_rbr"]
    # # ycols = ["sparseness_random", "sparseness_greedy", "sparseness_blind_greedy"]
    ds_to_plot = ds_all_dims.loc[(2, "normal", 0.1)]
    ds_to_plot.reset_index(inplace=True)
    ds_to_plot.drop(columns="input_file", inplace=True)
    ds_to_plot = ds_to_plot.set_index(["n_points"])
    ds_to_plot = ds_to_plot.sort_index(ascending=True)
    print(ds_to_plot.index)
    print(ds_to_plot.head(n=20))
    ds_to_plot.plot(y=ycols)
    plt.show()
