# Reduced Ordered Binary decision diagrams: distribution with respect to size 

## Dependencies

- The GNU MPFR Library: for arbitrary precision integers
  https://www.mpfr.org/

- FLINT : Fast Library for Number Theory FLINT
  Used for efficient polynomials arithmetics
  https://flintlib.org/

- The OpenMP API specification for parallel programming
  For parallelization 
  https://www.openmp.org/

## Compilation
With GCC:
```
gcc bdd_distribution.c -o bdd_distribution -lflint -lmpfr -lgmp -lm -Wall -fopenmp
```

# Usage
Under linux, in a terminal, in order to get distribution of ROBDDs for size for 5 variables, type 
```
./bdd_distribution 5
```
Computation is distributed on available cores thanks to library OPENMP.

On can also specify the number of cores allocated to the problem (when available) by adding a second argument. To specify 8 cores for instance, type

```
./bdd_distribution 5 8
```

## Example
```
$ ./bdd_distribution 5
max size=17, max profile: [1 2 4 8 2 ]
level 1 (T=15, B=2, nb_nodes=2) ................ in 0.00047 seconds (approx)
level 2 (T=7, B=10, nb_nodes=8) ........ in 9.1e-05 seconds (approx)
level 3 (T=3, B=14, nb_nodes=4) .... in 3.3e-05 seconds (approx)
level 4 (T=1, B=16, nb_nodes=2) .. in 1.5e-05 seconds (approx)
level 5 (T=0, B=17, nb_nodes=1) . in 1e-05 seconds (approx)
0	2
1	10
2	80
3	580
4	3920
5	24940
6	148832
7	819274
8	4077440
9	18038498
10	69381840
11	223877520
12	572592240
13	1074728520
14	1281360960
15	806420160
16	223534080
17	19958400
```
