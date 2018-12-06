--delete from spanner_results_euclidean;


copy spanner_results_euclidean(spanner_method, dim, n_points, epsilon, points_generation_method, n_edges, sparseness)
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/spanner_experiment_data/blind_greedy_spanner_from_server.txt' delimiter ';' ;



delete from spanner_results_cpp;


copy spanner_results_cpp(dim, points_generation_method, input_file, spanner_method,  epsilon, n_points, n_edges, sparseness, n_edges_to_n_points_ratio)
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/spanner_experiment_data/spanner_cpp_logs.txt' delimiter ';' csv header;


copy spanner_results_cpp(dim, points_generation_method, input_file, spanner_method,  epsilon, n_points, n_edges, sparseness, n_edges_to_n_points_ratio)
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/spanner_experiment_data/diff_densities_log.txt' delimiter ';' csv header;


copy spanner_results_euclidean_wspd(dim, points_generation_method, input_file, spanner_method,  epsilon, n_points, n_edges, n_distances_computed, sparseness, n_edges_to_n_points_ratio)
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/spanner_experiment_data/all_wspd_logs_server.txt' delimiter ';' csv header;

