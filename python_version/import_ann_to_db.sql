--delete from ann_results_euclidean;
--where n_points >= 20000;

copy ann_results_euclidean_single(dim, epsilon, n_points, points_generation_method, query_generation_method, min_fraction, max_fraction, avg_fraction, min_distances_number, max_distances_number, average_distances_number) 
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/ann_experiment_data/ann_results_single_tu_server.txt' delimiter ';' csv header


select count(*)
from ann_results_euclidean


--delete from ann_results_euclidean;
--where n_points >= 20000;

copy ann_results_euclidean_density(dim, epsilon, n_points, points_generation_method, query_generation_method, min_fraction, max_fraction, avg_fraction, min_distances_number, max_distances_number, average_distances_number) 
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/ann_experiment_data/nn_mcgill_original_results.csv' delimiter ';' csv 


copy ann_results_mcgill(fname, epsilon, n_points, min_fraction, max_fraction, avg_fraction, min_distances_number, max_distances_number, average_distances_number) 
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/ann_experiment_data/nn_mcgill_original_results.csv' delimiter ';' csv 


update ann_results_mcgill t
set wasserstein_degree = case 
                           when fname like '%q_1%' then 1 
                           when fname like '%q_2%' then 2 
                           when fname like '%q_3%' then 3 
                         end,
   dim = case                when fname like '%dim_0%' then 0
                           when fname like '%dim_1%' then 1 
                           when fname like '%dim_2%' then 2 
                         end
where t.fname like '%/%'

select *
from ann_results_mcgill t
where t.fname like '%/%'
set 


select count(*)
from ann_results_euclidean