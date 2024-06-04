# Single-Cell Protein Normalization

minimizes the variance of library sizes across cells + the absolute values of covariances between every protein and the normalization factor

the second constraint basically assumes the normalization factor to be independent of the protein counts.

So far this only works for small problems with the help of the cplex library, and it is a vry simple python program.
=> Update: we will implement the problem in c++ and solve it with e.g. CGAL
