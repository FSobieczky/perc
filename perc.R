##################################################
# perc.R
#
# Florian Sobieczky - Nov 17, 2015
#
##################################################

library(igraph)   # graph library
library(expm)     # powers of matrices


   # This makes a bond-percolative subgraph of a lattice on {1, ..., N}^d
   #  where the edge-retention probability is p:
make_perc <- function(N, d, p, colored=FALSE){   
    g<-make_lattice(length = N, dim = d);

    a<-rep("keep", length(E(g)));
    for(e in E(g)) a[e]<-ifelse(runif(1)<p,"black","white");
    E(g)$color <- a;

    if(colored==TRUE){
        return(g);
    }
    else{
        h<-delete_edges(g, which(a=="white"));
        return(h);
    }
}

   # This lets a delayed random walk (DRW) progress n-steps from
   # initial-vertex. The graph is vertex-colored according to
   # the heat-distribution (=transition probability of the DRW):
run_DRW<-function(g,n_steps,initial_vertex,gauge_color=2){
    A<-as.matrix(as_adj(g));
    Delta<-max(degree(g));
    id<-diag(1,length(V(g)));
    D<-diag(0,length(V(g)));
    for(v in V(g)){D[v,v]<-degree(g,v)};
    P<-as.matrix((1/Delta)*(A+(Delta*id-D)));  # DRW

    # Note: P^n_xy = Prob[ X_n = y | X_0 = y ]

    u<-rep(0.0,length(V(g)));
    u[initial_vertex]<-1.0;     # Initial Distribution:
                                # Singleton at initial vertex
    
    w<-rep(0, length(V(g)));
    w<- t(u) %*% (P%^%n_steps);  # This will be the transition
                                 # probab after n_steps steps.


    color_map<-c("yellow","orange","red","purple","blue")

    for(v in V(g)){              # Coloring ~ 2^(-n)
        V(g)$color[v]<-"yellow";
        for(i in 1:5){           # Use gauge_color to adjust coloring
            if(w[v]<0.5^(i+gauge_color)){ 
                V(g)$color[v]<-color_map[i];
            }
        }
    }

    return(g);
}

run_SRW<-function(g,n_steps,initial_vertex,gauge_color=2){
    A<-as.matrix(as_adj(g));
    Delta<-max(degree(g));
    id<-diag(1,length(V(g)));
    DD<-diag(0,length(V(g)));
    for(v in V(g)){
        if(degree(g,v)==0){
            DD[v,v]<-0;
        }
        else{
            DD[v,v]<-1/degree(g,v);
        }
    }

    
    P<-as.matrix(DD %*% A );      # SRW


    # Note: P^n_xy = Prob[ X_n = y | X_0 = y ]

    u<-rep(0.0,length(V(g)));
    u[initial_vertex]<-1.0;     # Initial Distribution:
                                # Singleton at initial vertex
    
    w<-rep(0, length(V(g)));
    w<- t(u) %*% (P%^%n_steps);  # This will be the transition
                                 # probab after n_steps steps.

    color_map<-c("yellow","orange","red","purple","blue")

    for(v in V(g)){              # Coloring ~ 2^(-n)
        V(g)$color[v]<-"yellow";
        for(i in 1:5){           # Use gauge_color to adjust coloring
            if(w[v]<0.5^(i+gauge_color)){
                V(g)$color[v]<-color_map[i];
            }
        }
    }

    return(g);
}


# Run like this (e.g.):
#
# g <- make_perc(10, 2, 0.95)#
# h <- run_DRW(g, 10, 1)
# tkplot(h)
#
#
# If Graph is large (N large), then
# larger value of gauge_color might
# be picked larger.

# Note: run_DRW and run_DRW can be used
#  for arbitrary undirected graphs.

