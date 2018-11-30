SELECT spanner_method, dim, n_points, epsilon, points_generation_method, 
       input_file, n_edges, sparseness, n_edges_to_n_points_ratio
  FROM public.spanner_results_cpp;


select distinct spanner_method
from public.spanner_results_cpp;

select distinct points_generation_method
from public.spanner_results_cpp;


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
       corr(n_edges, n_points*(n_points) / 2.0) as cdelete from spanner_results_cpp;

copy spanner_results_cpp(dim, points_generation_method, input_file, spanner_method,  epsilon, n_points, n_edges, sparseness, n_edges_to_n_points_ratio)
from '/home/narn/hdd_data/code/dynamic_spanner/python_version/spanner_experiment_data/spanner_cpp_logs.txt' delimiter ';' csv header;


orr_edges,
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
where 1 = 1
      --and spanner_method != 'greedy-non-blind'
      and abs(epsilon - 0.1) < 0.001
group by spanner_method, dim, epsilon, points_generation_method
having count(*) > 5
order by points_generation_method, epsilon, 5


select distinct *
from public.spanner_results_cpp t
where 1 = 1
---and abs(epsilon - 0.5) < 0.001
and t.spanner_method = 'blind-greedy'
and t.points_generation_method = 'unif'
and t.n_points = 700
order by dim
--and t.dim = 6




select distinct t2.n_points as "# points", 
t2.n_edges as "# edges (dim =2)", 
t3.n_edges as "# edges (dim =3)", 
t4.n_edges as "# edges (dim =4)"
from
(select n_points, n_edges
from public.spanner_results_cpp t
where 1 = 1
and abs(epsilon - 0.1) < 0.001
and t.spanner_method = 'blind-greedy'
and t.points_generation_method = 'unif'
and dim = 2) as t2,
(select n_points, n_edges
from public.spanner_results_cpp t
where 1 = 1
and abs(epsilon - 0.1) < 0.001
and t.spanner_method = 'blind-greedy'
and t.points_generation_method = 'unif'
and dim = 3) as t3,
(select n_points, n_edges
from public.spanner_results_cpp t
where 1 = 1
and abs(epsilon - 0.1) < 0.001
and t.spanner_method = 'blind-greedy'
and t.points_generation_method = 'unif'
and dim = 4) as t4
where 1 = 1
and t2.n_points = t3.n_points
and t3.n_points = t4.n_points
--and t4.n_points = t5.n_points
order by 1





SELECT spanner_method, dim, epsilon, substr(points_generation_method, 1, 4),
       regr_slope(log(n_edges), log(n_points)) as coeff_points,
       count(*) as cnt
FROM public.spanner_results_cpp t
where 1 = 1
--and dim = 2
      and spanner_method in ( 'greedy-non-blind', 'blind-greedy')
      and points_generation_method like 'unif%'
      and abs(epsilon - 0.1) < 0.001
      and n_points > 150
group by spanner_method, dim, epsilon, substr(points_generation_method, 1, 4)
--having count(*) > 5
order by 1,2, 3

