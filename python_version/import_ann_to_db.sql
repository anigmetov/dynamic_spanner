delete from ann_results_euclidean
where n_points >= 20000;

copy ann_results_euclidean(dim, epsilon, n_points, points_generation_method, query_generation_method, min_fraction, max_fraction, avg_fraction, min_distances_number, max_distances_number, average_distances_number) 
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/ann_results.txt' delimiter ';' csv header 1;