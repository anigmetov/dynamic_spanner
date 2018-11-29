--delete from ann_results_euclidean;
--where n_points >= 20000;

copy ann_results_euclidean_single(dim, epsilon, n_points, points_generation_method, query_generation_method, min_fraction, max_fraction, avg_fraction, min_distances_number, max_distances_number, average_distances_number) 
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/ann_experiment_data/ann_results_single_tu_server.txt' delimiter ';' csv header


select count(*)
from ann_results_euclidean


--delete from ann_results_euclidean;
--where n_points >= 20000;

copy ann_results_euclidean_density(dim, epsilon, n_points, points_generation_method, query_generation_method, min_fraction, max_fraction, avg_fraction, min_distances_number, max_distances_number, average_distances_number) 
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/ann_experiment_data/ann_density_dependency.txt' delimiter ';' csv 


select count(*)
from ann_results_euclidean