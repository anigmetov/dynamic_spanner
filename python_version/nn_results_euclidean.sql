SELECT dim, n_points, points_generation_method, query_generation_method, 
       min_fraction, max_fraction, avg_fraction,
       corr(log(n_points), avg_fraction)
  FROM public.nn_results_euclidean
  where dim = 20
  and query_generation_method = 'query_1'
  --and points_generation_method = 'get_points'
  and points_generation_method = 'get_uniform_points'
  --and points_generation_method = 'get_exponential_points'
  



SELECT dim, points_generation_method, query_generation_method, 
       corr(log(n_points), avg_fraction), corr(log(n_points), max_fraction)
  FROM public.nn_results_euclidean
  group by dim, points_generation_method, query_generation_method
  


SELECT dim, points_generation_method, query_generation_method, 
       corr(log(n_points), n_points * avg_fraction), corr(log(n_points), n_points * max_fraction),
       regr_slope(n_points * avg_fraction, log(n_points))
  FROM public.nn_results_euclideanSELECT dim, points_generation_method, query_generation_method, 
       corr(log(n_points), n_points * avg_fraction), corr(log(n_points), n_points * max_fraction),
       regr_slope(n_points * avg_fraction, log(n_points))
  FROM public.nn_results_euclidean
  group by dim, points_generation_method, query_generation_method
  order by 3, 2, 1


  group by dim, points_generation_method, query_generation_method
  order by 3, 2, 1


select distinct points_generation_method
  FROM public.nn_results_euclidean


2;"get_exponential_points";"query_1";0.920514898259124;0.926059226453816;4.08750514370556
3;"get_exponential_points";"query_1";0.943680436311807;0.92073517737512;10.9891526368338
10;"get_exponential_points";"query_1";0.991860745211884;0.975062552092248;8.14279102550976
20;"get_exponential_points";"query_1";0.987920944940354;0.960812608175291;19.5838544591107
2;"get_points";"query_1";0.977276319330191;0.877371259137558;2.82563206941339
3;"get_points";"query_1";0.966237490379959;0.851722583358071;4.62125992774077
10;"get_points";"query_1";0.992248229116422;0.98380891627036;73.1010084066743
20;"get_points";"query_1";0.915240090628727;0.888394576355053;1054.87925843387
2;"get_uniform_points";"query_1";0.931185099567011;0.743761627999598;4.29327036092317
3;"get_uniform_points";"query_1";0.887219918928688;0.554577089051824;5.81368762484344
10;"get_uniform_points";"query_1";0.979059595492524;0.974992267659215;98.8873577158166
20;"get_uniform_points";"query_1";0.938325029801143;0.924018434326418;1056.29840206586
2;"get_exponential_points";"query_2";0.981269571297237;0.973023709206457;3.92905692080768
3;"get_exponential_points";"query_2";0.987190025556414;0.982863208104619;6.46233401295457
10;"get_exponential_points";"query_2";0.991822312833497;0.98662616875532;8.11738122346986
20;"get_exponential_points";"query_2";0.991677168868015;0.971317216883968;18.5500768487854
2;"get_points";"query_2";0.930263514977635;0.921447399301354;2.44655362256516
3;"get_points";"query_2";0.95738367161963;0.967489134003678;3.82322437203775
10;"get_points";"query_2";0.995009608497848;0.938284351112862;68.6844901495992
20;"get_points";"query_2";0.940954857646094;0.943681368056374;972.356746429621
2;"get_uniform_points";"query_2";0.89216341003453;0.715902415280257;3.70360000411543
3;"get_uniform_points";"query_2";0.895780737365563;0.725625295033562;10.546289459134
10;"get_uniform_points";"query_2";0.955330480716474;0.897197044233469;138.46071712456
20;"get_uniform_points";"query_2";0.954626789039536;0.945112934511715;930.428542999842
  
