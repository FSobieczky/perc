# perc

perc.R

# Contains:

make_perc(N, d, p);
run_DRW(g, n_steps, initial_vertex, gauge_color = 2);
run_SRW(g, n_steps, initial_vertex, gauge_color = 2);


# Run like this (e.g.):                                                                               
#                                                                                                     
g <- make_perc(10, 2, 0.95)#                                                                        
h <- run_DRW(g, 10, 1)                                                                              
tkplot(h)                                                                                           
# This produces graphs like in 'Heat.png'.

#                                                                                                     
# If Graph is large (N large), then                                                                   
# larger value of gauge_color might                                                                   
# be picked larger.                                                                                   

# Note: run_DRW and run_DRW can be used                                                               
#  for arbitrary undirected graphs.         
#
# 'SRW_vs_DRW.png' contains comparison of
run_SRW(g,10,1)
# versus:
run_DRW(g,10,1)
# on a 5 by 5 bond-percolation subgraph.
#
# One can see well the parity issues due to
# the lattice being a bipartite graph.

