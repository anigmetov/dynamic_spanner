--delete from spanner_results_euclidean;


copy spanner_results_euclidean(spanner_method, dim, n_points, epsilon, points_generation_method, n_edges, sparseness)
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/spanner_experiment_data/blind_greedy_spanner_from_server.txt' delimiter ';' ;


