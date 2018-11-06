import pandas as pd
import matplotlib.pyplot as plt

ds = pd.read_csv("all_data.txt", header=0, sep=';')

index_cols = ['dim', 'n_points', 'point_method', 'epsilon']
ds_greedy = ds.loc[(ds['spanner_method'] == 'greedy')].set_index(index_cols, drop=False)
ds_random = ds.loc[(ds['spanner_method'] == 'random')].set_index(index_cols, drop=False)
ds_quasi = ds.loc[(ds['spanner_method'] == 'quasi')].set_index(index_cols, drop=False)
ds_blind_greedy = ds.loc[(ds['spanner_method'] == 'blind_greedy')].set_index(index_cols, drop=False)

ds_greedy.drop(columns=['spanner_method'], inplace=True)
ds_random.drop(columns=['spanner_method'], inplace=True)
ds_quasi.drop(columns=['spanner_method'], inplace=True)
ds_blind_greedy.drop(columns=['spanner_method'], inplace=True)

ds_all_dims = ds_greedy.join(ds_random, lsuffix='_greedy', rsuffix='_random').join(ds_quasi, rsuffix='_quasi').join(
    ds_blind_greedy, rsuffix='_blind_greedy')

print(ds_all_dims.head())
print(ds_all_dims.describe(include='all').to_string())

ycols = ["sparseness_random", "sparseness_greedy", "sparseness", "sparseness_blind_greedy"]
# ycols = ["sparseness_random", "sparseness_greedy", "sparseness_blind_greedy"]
mask_1 = (ds_all_dims['dim'] == 2) & (ds_all_dims['point_method'] == 'get_exponential_points') & (ds_all_dims['epsilon'] == 0.01)
ds_all_dims.loc[mask_1].plot(x="n_points", y=ycols)
# mask_1 = (ds_all_dims['n_points'] == 200) & (ds_all_dims['point_method'] == 'get_exponential_points') & (ds_all_dims['epsilon'] == 0.01)
# ds_all_dims.loc[mask_1].plot(x="dim", y=ycols)
plt.show()
