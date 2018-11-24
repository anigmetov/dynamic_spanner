SELECT spanner_method, dim, n_points, epsilon, points_generation_method, 
       input_file, n_edges, sparseness, n_edges_to_n_points_ratio
  FROM public.spanner_results_cpp;


select distinct spanner_method
from public.spanner_results_cpp;

select distinct points_generation_method
from public.spanner_results_cpp;

"exp"
"clustered"
"unif"
"normal"
"unif_100"


SELECT spanner_method, dim, n_points, epsilon, points_generation_method, 
       n_edges, sparseness, n_edges_to_n_points_ratio
FROM public.spanner_results_cpp t
where spanner_method = 'greedy-non-blind'
and abs(epsilon - 0.1) < 0.001
and dim = 2
and t.points_generation_method = 'unif'


SELECT spanner_method, dim, epsilon, points_generation_method, 
       corr(n_edges, n_points),
       regr_slope(n_edges, n_points)
FROM public.spanner_results_cpp t
where spanner_
group by spanner_method, dim, epsilon, points_generation_method
order by epsilon, points_generation_method, spanner_method