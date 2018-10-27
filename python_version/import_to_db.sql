delete from nn_results_euclidean
where n_points >= 20000;

	copy nn_results_euclidean(dim, n_points, points_generation_method, query_generation_method, min_fraction, max_fraction, avg_fraction) 
	from '/home/narn/hdd_data/code/dynamic_spanner/python_version/nn_more_results.txt' delimiter ';' ;