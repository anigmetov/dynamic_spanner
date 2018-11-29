SELECT dim, epsilon, points_generation_method, query_generation_method, 
	regr_slope(min_distances_number, log(n_points)) as min_dist_regr_degree,
	regr_slope(max_distances_number, log(n_points)) as max_dist_regr_degree,
	regr_slope(average_distances_number, log(n_points)) as avg_dist_regr_degree,
	count(*)
FROM public.ann_results_euclidean
where points_generation_method like '%unif%' and query_generation_method = 'query_1'
group by dim, epsilon, points_generation_method, query_generation_method
order by dim, epsilon, points_generation_method, query_generation_method


SELECT dim, epsilon, points_generation_method, query_generation_method, 
	--regr_slope(min_distances_number, log(n_points)) as min_dist_regr_degree,
	--regr_slope(max_distances_number, log(n_points)) as max_dist_regr_degree,
	regr_slope(average_distances_number, log(n_points)) as avg_dist_regr_degree1,
	corr(average_distances_number, log(n_points)) as avg_dist_corr,
	regr_slope(log(average_distances_number), log(log(n_points))) as avg_dist_regr_degree2,
	count(*)
FROM public.ann_results_euclidean
where abs(epsilon - 0.01) < 0.00001
--points_generation_method = 'unif' and query_generation_method = 'query_1'
group by dim, epsilon, points_generation_method, query_generation_method
order by 7, dim, epsilon, points_generation_method, query_generation_method



SELECT *
FROM public.ann_results_euclidean
where dim = 20 
      and epsilon < 0.49 
      and points_generation_method = 'unif'
      

SELECT dim, epsilon, points_generation_method, query_generation_method, 
	--regr_slope(min_distances_number, log(n_points)) as min_dist_regr_degree,
	--regr_slope(max_distances_number, log(n_points)) as max_dist_regr_degree,
	regr_slope(average_distances_number, log(n_points)) as avg_dist_regr_degree1,
	corr(average_distances_number, log(n_points)) as avg_dist_corr,
	regr_slope(log(average_distances_number), log(log(n_points))) as avg_dist_regr_degree2,
	count(*)
FROM public.ann_results_euclidean
--where
-- points_generation_method = 'unif' and query_generation_method = 'query_1'
group by dim, epsilon, points_generation_method, query_generation_method
order by epsilon, dim, points_generation_method, query_generation_method

SELECT dim, epsilon, points_generation_method, query_generation_method, 
count(*),
min(average_distances_number / n_points),
max(average_distances_number / n_points)
FROM public.ann_results_euclidean
group by dim, epsilon, points_generation_method, query_generation_method
order by dim, epsilon, points_generation_method, query_generation_method


SELECT dim, epsilon, points_generation_method as points_method, query_generation_method as query_method, n_points, avg(average_distances_number) as dist_num
FROM public.ann_results_euclidean_single
where 1 = 1 
and dim in (2, 5, 10)
and n_points < 21000
group by dim, epsilon, points_generation_method, query_generation_method, n_points
having count(*) > 5
order by dim, n_points


SELECT dim, n_points, avg(average_distances_number) as dist_num
FROM public.ann_results_euclidean_single
where 1 = 1 
and dim in (2, 5, 10)
and n_points < 21000
group by dim, epsilon, points_generation_method, query_generation_method, n_points
having count(*) > 5
order by dim, n_points


SELECT epsilon, n_points, points_generation_method, avg(average_distances_number) as dist_num
FROM public.ann_results_euclidean_density
where n_points = 1000
group by epsilon, n_points, points_generation_method
order by 1, 2, 3, 4




SELECT dim, epsilon, points_generation_method, query_generation_method, 
	regr_slope(log(average_distances_number), log(log(n_points))) as avg_dist_regr_degree2,
	count(*)
FROM public.ann_results_euclidean
where points_generation_method like '%unif%' and query_generation_method = 'query_1'
group by dim, epsilon, points_generation_method, query_generation_method
order by epsilon, dim, points_generation_method, query_generation_method
