SELECT dim, epsilon, points_generation_method, query_generation_method, 
	regr_slope(min_distances_number, log(n_points)) as min_dist_regr_degree,
	regr_slope(max_distances_number, log(n_points)) as max_dist_regr_degree,
	regr_slope(average_distances_number, log(n_points)) as avg_dist_regr_degree,
	count(*)
FROM public.ann_results_euclidean
where points_generation_method = 'unif' and query_generation_method = 'query_1'
group by dim, epsilon, points_generation_method, query_generation_method
order by dim, epsilon, points_generation_method, query_generation_method


SELECT dim, epsilon, points_generation_method, query_generation_method, 
	regr_slope(min_distances_number, log(n_points)) as min_dist_regr_degree,
	regr_slope(max_distances_number, log(n_points)) as max_dist_regr_degree,
	regr_slope(average_distances_number, log(n_points)) as avg_dist_regr_degree,
	count(*)
FROM public.ann_results_euclidean
where points_generation_method = 'unif' and query_generation_method = 'query_1'
group by dim, epsilon, points_generation_method, query_generation_method
order by dim, epsilon, points_generation_method, query_generation_method



SELECT *
FROM public.ann_results_euclidean
where dim = 20 
      and epsilon < 0.49 
      and points_generation_method = 'unif'
      

