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


truncate table public.spanner_results_wasserstein_mcgill
copy  public.spanner_results_wasserstein_mcgill(
            input_file, spanner_method, epsilon, n_points,  n_edges, sparseness, n_edges_to_n_points_ratio)
from '/home/narn/hdd_data/code/dynamic_spanner/data_files/blind_spanner_mcgill_original_results_450_max.txt' delimiter ';' csv;         


copy  public.spanner_results_wasserstein_mcgill(
            input_file, spanner_method, epsilon, n_points,  n_edges, sparseness, n_edges_to_n_points_ratio)
from '/home/narn/hdd_data/code/dynamic_spanner/data_files/mcgill_original_correct_results_450_all.txt' delimiter ';' csv;         


update  public.spanner_results_wasserstein_mcgill
set q = case
when input_file like '%q_1%' then 1
when input_file like '%q_2%' then 2
when input_file like '%q_3%' then 3
when input_file like '%q_4%' then 4
end,
dim = case
when input_file like '%dim_1%' then 1
when input_file like '%dim_2%' then 2
when input_file like '%dim_0%' then 0
end

