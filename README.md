# treeCentrality

treeCentrality is a package for computing network science statistics on (rooted or unrooted) phylogenetic trees in linear 
time and space. The statistics include diameter, mean path length, degree assortativity and three node-based centrality values: 
betweenness centrality, closeness centrality, and eigenvector centrality. These statistics, useful for distinguishing 
evolutionary scenarios, are computed in time and space linear in the tree size, and can take branch length into account.  

In addition, this package can compute the spectra of the adjacency, Laplacian, and distance matrices, although this computation 
may be cubic in the tree size. The package can also compute a number of classical topology statistics, such as the Colless and 
Sackin indices, numbers of small configurations such as cherries, the heights of tree nodes, and “staircaseness” statistics. 
The methods are described in “Network science inspires novel tree shape statistics”.
The paper is currently available at https://www.biorxiv.org/content/10.1101/608646v1.

You can install the package from github via:

devtools::install_github('Leonardini/treeCentrality')
