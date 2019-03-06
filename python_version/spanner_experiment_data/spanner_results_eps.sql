select tg.dim, tg.epsilon as eps, tg.n_edges as greedy_non_blind_edges, tb.n_edges as blind_greedy_edges,
tb.n_edges  / tg.n_edges as ratio

from (
SELECT dim, epsilon, avg(n_edges) as n_edges
  FROM public.spanner_results_eps
  where spanner_method = 'greedy-non-blind'
  group by dim, epsilon ) tg join
(
SELECT dim, epsilon, avg(n_edges) as n_edges
  FROM public.spanner_results_eps
  where spanner_method = 'blind-greedy'
  group by dim, epsilon ) tb on (tg.dim = tb.dim and tb.epsilon = tg.epsilon)
order by tg.dim, tg.epsilon  