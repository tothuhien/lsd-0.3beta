# lsd-0.3beta
LSD: LEAST-SQUARES METHODS TO ESTIMATE RATES AND DATES FROM SERIAL PHYLOGENIES - v0.3beta by Thu-Hien TO

If you use this software, please cite: “ Fast dating using least-squares criteria and algorithms”, T-H. To, M. Jung, S. Lycett, O. Gascuel, Syst Biol. 2016 Jan;65(1):82-97.


How To Compile LSD:
     Use C++ compiler and library support for the ISO C++ 2011 to compile the program from the code source. From the folder src, type 'make':
     
How to Run LSD:
	After compiling the program, you have the lsd executable file in the src folder, then run it from a terminal. You can also use the binary files in the bin folder but sometimes it might not work on some computers that have different configuration as the machine compiled them. From the directory that contains the executable file:
	
		if you want to use use the interface, type ./lsd
		
		if you want to use the command line, type ./lsd -i <"input_tree_file"> -d <"input_date_file"> (to estimate absolute dates)
		
		                                       or ./lsd -i <"input_tree_file"> -a root_date -z leaves_date (to estimate relative dates)
		                                       
			further options can be specified (-f -o -c -v -n -r -g ...). Type "./lsd -h" for help. Option -c is recommended to take into account the temporal constraints (date of a node >= date of its ancestors).

Some examples:

1) Examples of command line:

for rooted tree, constrained mode, and using variances
./lsd -i rootedtree_file -d date_file -c -v 1

for rooted tree, constrained mode, re-estimate the root position around the given root, and non variances
./lsd -i rootedtree_file -d date_file -c -r l

with the previous example, if you want to calculate confidence intervals from 100 simulated trees
./lsd -i rootedtree_file -d date_file -c -r l -f 100

for unrooted tree without outgroups, without constraints, estimate the root position, and using variances
./lsd -i unrootedtree_file -d date_file -c -v 1 -r a

for unrooted tree with outgroups, constrained mode, using variances from the estimated branch lengths (run LSD twice), remove outgroups to obtain the root
./lsd -i unrootedtree_file -d date_file -g outgroup_file -c -v 2

similary but keep outgroups and estimate the root on the branch defined by them
./lsd -i unrootedtree_file -d date_file -g outgroup_file -k -c -v 2

for rooted tree, constrained mode, using variances, and using given rates to estimate dates
./lsd -i rootedtree_file -d date_file -w given_rate_file -c -v 1

for rooted tree, estimating relative dates with date root=0 and date of all leaves=1, using variances, with constraint
./lsd -i tree_file -c -v 1 -a 0 -z 1


2) Examples of input files:

Example of Input_tree_file format (newick, can be either binary or polytomy)

((a:0.12,d:0.12):0.3,(b:0.3,c:0.5):0.4);

((a:0.12,b:0.3):0.7,(c:0.5,d:0.8):0.1);

Example of Input_date_file format (it's not necessary to give the temporal constraints for all tips):

  4                     # number of temporal constraints
  a 1999                # the date of a is 1999
  b 2000                # the date of b is 2000
  c l(1780)             # the date of c is at least 1780
  mrca(a,b,c) u(2000)   # the date of the most recent ancestor of a,b, and c is at most 2000
  d b(1998,2000)        # the date of d is between 1998 and 2000

You can also define the labels for internal nodes and use them to define their temporal constraints.
For example you have an input tree: ((a:0.12,d:0.12)n1:0.3,(b:0.3,c:0.5)n2:0.4)root;
You can have an input date file as follows:

  5
  a 2000
  n1 l(2001)
  c b(2001,2004)
  n2 u(2003)
  root b(1998,1999)

Example of Outgroup_file format:
	2
	outgroup1
	outgroup2

Example of given_rate_file format:
	0.0068
	0.0052

3) Output files: 

.result : contain the estimated rate, root date and the value of the objective function.

.newick : trees in newick format with the new branch length (re-estimated by the program).

.date.newick : trees in newick format where branch lengths are measured rescaled to time unit by multiplying with the estimated rate. 

.nexus : trees in nexus format which contain information about the dates of internal nodes (named date), branch lengths, and the confidence intervals (named CI) if option -f was used.
