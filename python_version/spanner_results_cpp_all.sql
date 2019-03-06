SELECT distinct spanner_method
       
  FROM public.spanner_results_cpp
 ;

select coalesce(t_gnb.dim, t_wspd.dim, t_bg.dim) as dim, 
coalesce(t_gnb.n_points, t_wspd.n_points, t_bg.n_points) as n_points,
coalesce(t_gnb.epsilon, t_wspd.epsilon, t_bg.epsilon) as epsilon, 
coalesce(t_gnb.pgm, t_wspd.pgm, t_bg.pgm) as pgm, 
    t_gnb.n_edges as edges_greedy, t_bg.n_edges as edges_blind_greedy, 
    t_wspd.n_edges as edges_wspd,
    t_qsg.n_edges as edges_quasi_sorted_greedy,
    t_qss.n_edges as edges_quasi_sorted_shaker,
    
    t_gnb.sparseness as sparseness_greedy,
    t_wspd.sparseness as sparseness_wspd,
    t_bg.sparseness as sparseness_blind_greedy
    
from 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges,
       avg(sparseness) as sparseness,
       avg(n_edges_to_n_points_ratio) as n_edges_to_n_points_ratio
FROM public.spanner_results_cpp t
where t.spanner_method = 'greedy-non-blind'
group by dim, n_points, epsilon, points_generation_method ) t_gnb
full outer join
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges,
       avg(sparseness) as sparseness,
       avg(n_edges_to_n_points_ratio) as n_edges_to_n_points_ratio
FROM public.spanner_results_euclidean_wspd t
--where t.spanner_method = 'greedy-non-blind'
group by dim, n_points, epsilon, points_generation_method ) t_wspd on (t_gnb.dim = t_wspd.dim and t_gnb.n_points = t_wspd.n_points and t_gnb.pgm = t_wspd.pgm and t_gnb.epsilon = t_wspd.epsilon)
full outer join (SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges,
       avg(sparseness) as sparseness,
       avg(n_edges_to_n_points_ratio) as n_edges_to_n_points_ratio
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-greedy'
group by dim, n_points, epsilon, points_generation_method ) t_bg on (t_gnb.dim = t_bg.dim and t_gnb.n_points = t_bg.n_points and t_gnb.pgm = t_bg.pgm and t_gnb.epsilon = t_bg.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-quasi-sorted-greedy'
group by dim, n_points, epsilon, points_generation_method ) t_qsg on (t_gnb.dim = t_qsg.dim and t_gnb.n_points = t_qsg.n_points and t_gnb.pgm = t_qsg.pgm and t_gnb.epsilon = t_qsg.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-quasi-sorted-shaker'
group by dim, n_points, epsilon, points_generation_method ) t_qss on (t_gnb.dim = t_qss.dim and t_gnb.n_points = t_qss.n_points and t_gnb.pgm = t_qss.pgm and t_gnb.epsilon = t_qss.epsilon)
where t_wspd.dim = 4
order by 4, 1, 3, 2

 
select coalesce(t_gnb.dim, t_wspd.dim, t_bg.dim) as dim, 
coalesce(t_gnb.n_points, t_wspd.n_points, t_bg.n_points) as n_points,
coalesce(t_gnb.epsilon, t_wspd.epsilon, t_bg.epsilon) as epsilon, 
coalesce(t_gnb.pgm, t_wspd.pgm, t_bg.pgm) as pgm, 
    --t_gnb.n_edges  as edges_greedy, 
    --t_bg.n_edges as edges_blind_greedy, 
    --t_wspd.n_edges as edges_wspd,
    --t_qsg.n_edges as edges_quasi_sorted_greedy,
    --t_qss.n_edges as edges_quasi_sorted_shaker,
    --t_rbr.n_edges as edges_blind_random,
    --t_rbr_cf.n_edges as edges_blind_random_connect_first,
    --t_rbr_lb.n_edges as edges_blind_random_lower_bound_first,
    --t_rbr_cf_lb.n_edges as edges_blind_random_connect_first_lower_bound_first,

    t_gnb.n_edges / t_gnb.n_points as ratio_greedy, 
    t_bg.n_edges / t_gnb.n_points as ratio_blind_greedy, 
    t_wspd.n_edges / t_gnb.n_points  as ratio_wspd,
    t_qsg.n_edges / t_gnb.n_points as ratio_quasi_sorted_greedy,
    t_qss.n_edges / t_gnb.n_points as ratio_quasi_sorted_shaker,
    t_rbr.n_edges / t_gnb.n_points as ratio_blind_random,
    t_rbr_cf.n_edges / t_gnb.n_points as ratio_blind_random_connect_first,
    t_rbr_lb.n_edges / t_gnb.n_points as ratio_blind_random_lower_bound_first,
    t_rbr_cf_lb.n_edges / t_gnb.n_points as ratio_blind_random_connect_first_lower_bound_first
    
from 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'greedy-non-blind'
group by dim, n_points, epsilon, points_generation_method ) t_gnb
join
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_euclidean_wspd t
group by dim, n_points, epsilon, points_generation_method ) t_wspd on (t_gnb.dim = t_wspd.dim and t_gnb.n_points = t_wspd.n_points and t_gnb.pgm = t_wspd.pgm and t_gnb.epsilon = t_wspd.epsilon)
full outer join (SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-greedy'
group by dim, n_points, epsilon, points_generation_method ) t_bg on (t_gnb.dim = t_bg.dim and t_gnb.n_points = t_bg.n_points and t_gnb.pgm = t_bg.pgm and t_gnb.epsilon = t_bg.epsilon)
join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-quasi-sorted-greedy'
group by dim, n_points, epsilon, points_generation_method ) t_qsg on (t_gnb.dim = t_qsg.dim and t_gnb.n_points = t_qsg.n_points and t_gnb.pgm = t_qsg.pgm and t_gnb.epsilon = t_qsg.epsilon)
join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-quasi-sorted-shaker'
group by dim, n_points, epsilon, points_generation_method ) t_qss on (t_gnb.dim = t_qss.dim and t_gnb.n_points = t_qss.n_points and t_gnb.pgm = t_qss.pgm and t_gnb.epsilon = t_qss.epsilon)
join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-random-bad-ratio-connect-first-lower-bound-first'
group by dim, n_points, epsilon, points_generation_method ) t_rbr_cf_lb on (t_gnb.dim = t_rbr_cf_lb.dim and t_gnb.n_points = t_rbr_cf_lb.n_points and t_gnb.pgm = t_rbr_cf_lb.pgm and t_gnb.epsilon = t_rbr_cf_lb.epsilon)
join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-random-bad-ratio-connect-first'
group by dim, n_points, epsilon, points_generation_method ) t_rbr_cf on (t_gnb.dim = t_rbr_cf.dim and t_gnb.n_points = t_rbr_cf.n_points and t_gnb.pgm = t_rbr_cf.pgm and t_gnb.epsilon = t_rbr_cf.epsilon)
join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-random-bad-ratio-lower-bound-first'
group by dim, n_points, epsilon, points_generation_method ) t_rbr_lb on (t_gnb.dim = t_rbr_lb.dim and t_gnb.n_points = t_rbr_lb.n_points and t_gnb.pgm = t_rbr_lb.pgm and t_gnb.epsilon = t_rbr_lb.epsilon)
join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-random-bad-ratio'
group by dim, n_points, epsilon, points_generation_method ) t_rbr on (t_gnb.dim = t_rbr.dim and t_gnb.n_points = t_rbr.n_points and t_gnb.pgm = t_rbr.pgm and t_gnb.epsilon = t_rbr.epsilon)
where t_rbr.dim = 4 
--and t_rbr.pgm = 'normal'
--and abs(t_rbr.epsilon - 0.1) < 0.0001
order by 4, 1, 3, 2



select coalesce(t_gnb.dim, t_wspd.dim, t_bg.dim) as dim, 
coalesce(t_gnb.n_points, t_wspd.n_points, t_bg.n_points) as n_points,
coalesce(t_gnb.epsilon, t_wspd.epsilon, t_bg.epsilon) as epsilon, 
coalesce(t_gnb.pgm, t_wspd.pgm, t_bg.pgm) as pgm, 
    t_gnb.n_edges  as edges_greedy, 
    t_bg.n_edges as edges_blind_greedy, 
    t_wspd.n_edges as edges_wspd
from 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'greedy-non-blind'
group by dim, n_points, epsilon, points_generation_method ) t_gnb
full outer join
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_euclidean_wspd t
group by dim, n_points, epsilon, points_generation_method ) t_wspd on (t_gnb.dim = t_wspd.dim and t_gnb.n_points = t_wspd.n_points and t_gnb.pgm = t_wspd.pgm and t_gnb.epsilon = t_wspd.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp tselect coalesce(t_gnb.dim, t_wspd.dim, t_bg.dim) as dim, 
coalesce(t_gnb.n_points, t_wspd.n_points, t_bg.n_points) as n_points,
coalesce(t_gnb.epsilon, t_wspd.epsilon, t_bg.epsilon) as epsilon, 
coalesce(t_gnb.pgm, t_wspd.pgm, t_bg.pgm) as pgm, 
    t_gnb.n_edges as edges_greedy, t_bg.n_edges as edges_blind_greedy, t_wspd.n_edges as edges_wspd,
    t_qsg.n_edges as edges_quasi_sorted_greedy,
    t_qss.n_edges as edges_quasi_sorted_shaker
    
from 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'greedy-non-blind'
group by dim, n_points, epsilon, points_generation_method ) t_gnb
full outer join
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_euclidean_wspd t
--where t.spanner_method = 'greedy-non-blind'
group by dim, n_points, epsilon, points_generation_method ) t_wspd on (t_gnb.dim = t_wspd.dim and t_gnb.n_points = t_wspd.n_points and t_gnb.pgm = t_wspd.pgm and t_gnb.epsilon = t_wspd.epsilon)
full outer join (SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-greedy'
group by dim, n_points, epsilon, points_generation_method ) t_bg on (t_gnb.dim = t_bg.dim and t_gnb.n_points = t_bg.n_points and t_gnb.pgm = t_bg.pgm and t_gnb.epsilon = t_bg.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-quasi-sorted-greedy'
group by dim, n_points, epsilon, points_generation_method ) t_qsg on (t_gnb.dim = t_qsg.dim and t_gnb.n_points = t_qsg.n_points and t_gnb.pgm = t_qsg.pgm and t_gnb.epsilon = t_qsg.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_cpp t
where t.spanner_method = 'blind-quasi-sorted-shaker'
group by dim, n_points, epsilon, points_generation_method ) t_qss on (t_gnb.dim = t_qss.dim and t_gnb.n_points = t_qss.n_points and t_gnb.pgm = t_qss.pgm and t_gnb.epsilon = t_qss.epsilon)
where t_wspd.dim = 4
order by 4, 1, 3, 2


where t.spanner_method = 'blind-greedy'
group by dim, n_points, epsilon, points_generation_method ) t_bg on (t_gnb.dim = t_bg.dim and t_gnb.n_points = t_bg.n_points and t_gnb.pgm = t_bg.pgm and t_gnb.epsilon = t_bg.epsilon)
where t_bg.dim = 2 
--and abs(t_bg.epsilon - 0.1) < 0.0001
and t_bg.n_points > 900
order by 4, 1, 3, 2




select t_unif.dim, t_unif.n_points, t_unif.epsilon, t_normal.n_edges as edges_normal, t_unif.n_edges as edges_unif, t_normal.n_points * (t_normal.n_points - 1) / 2
from (
SELECT dim, n_points, epsilon, 
--points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_euclidean_wspd t
where t.points_generation_method = 'normal'
group by t.dim, t.n_points, t.epsilon, t.points_generation_method ) t_normal join
(
SELECT dim, n_points, epsilon, 
--points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_euclidean_wspd t
where t.points_generation_method = 'unif'
group by t.dim, t.n_points, t.epsilon, t.points_generation_method ) t_unif on (t_normal.dim = t_unif.dim and t_normal.n_points = t_unif.n_points and t_normal.epsilon = t_unif.epsilon)
order by 1, 3, 2



SELECT *
FROM public.spanner_results_cpp t
order by n_edges - n_points


select regr_slope(log(t_normal.n_edges), log(t_normal.n_points))
from (
SELECT dim, n_points, epsilon, 
--points_generation_method as pgm, 
       avg(n_edges) as n_edges
FROM public.spanner_results_euclidean_wspd t
where t.points_generation_method = 'normal'
group by t.dim, t.n_points, t.epsilon, t.points_generation_method ) t_normal
where n_points > 700


update spanner_results_euclidean_wspd
set sparseness = n_edges / (n_points * (n_points - 1) / 2),
    n_edges_to_n_points_ratio = n_edges / n_points



 

select coalesce(t_gnb.dim, t_bg.dim) as dim, 
coalesce(t_gnb.n_points,  t_bg.n_points) as n_points,
coalesce(t_gnb.epsilon, t_bg.epsilon) as epsilon, 
coalesce(t_gnb.q, t_bg.q) as q, 
    --t_gnb.n_edges  as edges_greedy, 
    --t_bg.n_edges as edges_blind_greedy, 
    --t_wspd.n_edges as edges_wspd,
    --t_qsg.n_edges as edges_quasi_sorted_greedy,
    --t_qss.n_edges as edges_quasi_sorted_shaker,
    --t_rbr.n_edges as edges_blind_random,
    --t_rbr_cf.n_edges as edges_blind_random_connect_first,
    --t_rbr_lb.n_edges as edges_blind_random_lower_bound_first,
    --t_rbr_cf_lb.n_edges as edges_blind_random_connect_first_lower_bound_first,

    t_gnb.n_edges / t_gnb.n_points as ratio_greedy, 
    t_bg.n_edges / t_gnb.n_points as ratio_blind_greedy, 
    --t_qsg.n_edges / t_gnb.n_points as ratio_quasi_sorted_greedy,
    --t_qss.n_edges / t_gnb.n_points as ratio_quasi_sorted_shaker,
    --t_rbr.n_edges / t_gnb.n_points as ratio_blind_random,
    --t_rbr_cf.n_edges / t_gnb.n_points as ratio_blind_random_connect_first,
    t_rbr_lb.n_edges / t_gnb.n_points as ratio_blind_random_lower_bound_first,
    --t_rbr_cf_lb.n_edges / t_gnb.n_points as ratio_blind_random_connect_first_lower_bound_first,

    t_gnb.n_edges / t_gnb.max_edges as sparseness_greedy, 
    t_bg.n_edges / t_gnb.max_edges as sparseness_blind_greedy, 
    --t_qsg.n_edges / t_gnb.max_edges as sparseness_quasi_sorted_greedy,
    --t_qss.n_edges / t_gnb.max_edges  as sparseness_quasi_sorted_shaker,
    --t_rbr.n_edges / t_gnb.max_edges  as sparseness_blind_random,
    --t_rbr_cf.n_edges / t_gnb.max_edges  as sparseness_blind_random_connect_first,
    t_rbr_lb.n_edges / t_gnb.max_edges  as sparseness_blind_random_lower_bound_first
    --, t_rbr_cf_lb.n_edges / t_gnb.max_edges  as sparseness_blind_random_connect_first_lower_bound_first

    
from 
(SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges,
       n_points * (n_points -1 ) / 2 as max_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'greedy-non-blind'
group by dim, n_points, epsilon, q ) t_gnb
full outer join (SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'blind-greedy'
group by dim, n_points, epsilon, q ) t_bg on (t_gnb.dim = t_bg.dim and t_gnb.n_points = t_bg.n_points and t_gnb.q = t_bg.q and t_gnb.epsilon = t_bg.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'blind-quasi-sorted-greedy'
group by dim, n_points, epsilon, q ) t_qsg on (t_gnb.dim = t_qsg.dim and t_gnb.n_points = t_qsg.n_points and t_gnb.q = t_qsg.q and t_gnb.epsilon = t_qsg.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'blind-quasi-sorted-shaker'
group by dim, n_points, epsilon, q ) t_qss on (t_gnb.dim = t_qss.dim and t_gnb.n_points = t_qss.n_points and t_gnb.q = t_qss.q and t_gnb.epsilon = t_qss.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'blind-random-bad-ratio-connect-first-lower-bound-first'
group by dim, n_points, epsilon, q ) t_rbr_cf_lb on (t_gnb.dim = t_rbr_cf_lb.dim and t_gnb.n_points = t_rbr_cf_lb.n_points and t_gnb.q = t_rbr_cf_lb.q and t_gnb.epsilon = t_rbr_cf_lb.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'blind-random-bad-ratio-connect-first'
group by dim, n_points, epsilon, q ) t_rbr_cf on (t_gnb.dim = t_rbr_cf.dim and t_gnb.n_points = t_rbr_cf.n_points and t_gnb.q = t_rbr_cf.q and t_gnb.epsilon = t_rbr_cf.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'blind-random-bad-ratio-lower-bound-first'
group by dim, n_points, epsilon, q ) t_rbr_lb on (t_gnb.dim = t_rbr_lb.dim and t_gnb.n_points = t_rbr_lb.n_points and t_gnb.q = t_rbr_lb.q and t_gnb.epsilon = t_rbr_lb.epsilon)
full outer join 
(SELECT dim, n_points, epsilon, q as q, 
       avg(n_edges) as n_edges
FROM public.spanner_results_wasserstein_mcgill t
where t.spanner_method = 'blind-random-bad-ratio'
group by dim, n_points, epsilon, q ) t_rbr on (t_gnb.dim = t_rbr.dim and t_gnb.n_points = t_rbr.n_points and t_gnb.q = t_rbr.q and t_gnb.epsilon = t_rbr.epsilon)
where 1 = 1
  and t_gnb.dim =  0
and t_gnb.q = 1
--and t_gnb.n_points > 200
and abs(t_rbr.epsilon - 0.2) < 0.0001
order by 3, 1, 4, 2