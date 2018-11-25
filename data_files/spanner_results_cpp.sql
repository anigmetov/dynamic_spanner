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


SELECT spanner_method, dim, n_points, epsilon, points_generatiSELECT spanner_method, dim, n_points, epsilon, points_generation_method, 
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
       corr(n_edges, n_points) as corr_points,
       regr_slope(n_edges, n_points) as coeff_points,
       corr(n_edges, n_points*(n_points) / 2.0) as corr_edges,
       regr_slope(n_edges, n_points*(n_points) / 2.0) as coeff_edges,
       count(*) as cnt
FROM public.spanner_results_cpp t
--where spanner_method = 'greedy-non-blind'
group by spanner_method, dim, epsilon, points_generation_method
order by epsilon, points_generation_method, spanner_method

select *
from (
SELECT spanner_method, dim, epsilon, points_generation_method, 
       corr(n_edges, n_points) as corr_points,
       regr_slope(n_edges, n_points) as coeff_points,
       corr(n_edges, n_points*(n_points) / 2.0) as corr_edges,
       regr_slope(n_edges, n_points*(n_points) / 2.0) as coeff_edges,
       count(*) as cnt,
       min(sparseness) as min_sparseness,
       max(sparseness) as max_sparseness
FROM public.spanner_results_cpp t
--where spanner_method = 'greedy-non-blind'
group by spanner_method, dim, epsilon, points_generation_method) tt
where tt.corr_points < corr_edges 
order by epsilon, points_generation_method, spanner_methodon_method, 
       n_edges, sparseness, n_edges_to_n_points_ratio
FROM public.spanner_results_cpp t
where spanner_method = 'greedy-non-blind'
and abs(epsilon - 0.1) < 0.001
and dim = 2
and t.points_generation_method = 'unif'


SELECT spanner_method, dim, epsilon, points_generation_method, 
       corr(n_edges, n_points) as corr_points,
       regr_slope(n_edges, n_points) as coeff_points,
       corr(n_edges, n_points*(n_points) / 2.0) as corr_edges,
       regr_slope(n_edges, n_points*(n_points) / 2.0) as coeff_edges,
       count(*) as cnt
FROM public.spanner_results_cpp t
--where spanner_method = 'greedy-non-blind'
group by spanner_method, dim, epsilon, points_generation_method
order by epsilon, points_generation_method, spanner_method

select *
from (
SELECT spanner_method, dim, epsilon, points_generation_method, 
       corr(n_edges, n_points) as corr_points,
       regr_slope(n_edges, n_points) as coeff_points,
       corr(n_edges, n_points*(n_points) / 2.0) as corr_edges,
       regr_slope(n_edges, n_points*(n_points) / 2.0) as coeff_edges,
       count(*) as cnt,
       min(sparseness) as min_sparseness,
       max(sparseness) as max_sparseness,
       max(n_points) as max_points,
       min(n_points) as min_points
       
FROM public.spanner_results_cpp t
--where spanner_method = 'greedy-non-blind'
group by spanner_method, dim, epsilon, points_generation_method
having count(*) > 5) tt
where tt.corr_points < corr_edges 
order by epsilon, points_generation_method, spanner_method

select regr_slope(log(x*x), log(x))
from (
select generate_series as x
from generate_series(1, 1600)
) t

SELECT spanner_method, dim, epsilon, points_generation_method, 
       regr_slope(log(n_edges), log(n_points)) as coeff_points,
       count(*) as cnt,
       min(sparseness) as min_sparseness,
       max(sparseness) as max_sparseness,
       max(n_edges_to_n_points_ratio) as max_n_edges_to_n_points_ratio,
       min(n_edges_to_n_points_ratio) as min_n_edges_to_n_points_ratio
FROM public.spanner_results_cpp t
where spanner_method != 'greedy-non-blind'
group by spanner_method, dim, epsilon, points_generation_method
having count(*) > 5
order by points_generation_method, epsilon, 5


