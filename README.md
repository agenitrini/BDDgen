# BDDgen
BDD generator
(C) J. Clément and A. Genitrini under the GNU Public license v.3 (cf. LICENSE.txt)


BDDgen is an unranking algorithm for Reduced Ordered Binary Decision Diagrams (ROBDDs) of a given size and a given number of variables. The tool is based on the classical recursive method of Nijenhuis and Wilf, but the ROBDDs are not decomposable in their sense. Thus an extra set of structures called pool-BDDs must first be counted before the generation step.
The algorithm is based on the paper Binary Decision Diagrams: from Tree Compaction to Sampling that appeared in LATIN'20 conference (from J. Clément and A. Genitrini).

The first step of the algorithm is a pre-computation (done only once). Then the generation, using either the unranking approach, or a uniform random approach takes place.

----
![Tree example](https://github.com/agenitrini/BDDgen/blob/master/tests_dot/example.png)
----

On a standard PC the precomputation for the ROBDDs with at most 7 variables needs 6s and then each ROBDD generation uses around 1 ms. To deal with ROBDDs with 8 variables, the pre-computation step takes 12min (and 2.5 GiB of RAM) and then the generation of any structure necessitates around 8ms.

This python programs should be used (on a standard PC) for the generation of ROBDDs with at most 8 variables. With our prototype C++ implementation (under packaging) we reach the generation of ROBDDs with 9 variables (within several hours).  

----
The main file is Count_Gen.py. The others contain tool-functions.
