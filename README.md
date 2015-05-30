#pcluster#

Parallel hierarchical single-linkage clustering for DNA sequences. Input reading and verification, K-mer comparison, sequence alignment 
with SIMD, and some small utility methods were using the code from [swarm](https://github.com/torognes/swarm) with some modifications 
within these classes: `db`, `cpu_info`, `util`, `qgram`, `search8`, `search16`, and `scan`.

##Compiling##
To compile executable "main" for general usage (debug log disabled) use:

`make clean main`

To compile executable "test" for testing (debug log enabled and results are sort alphabetically by header) use:

`make clean test`

To compile both version use:

`make clean all`


##Unit Test##
To run short test, execute `test.sh`
To run full test (with valgrind), execute `run.sh`

The unit test compares the result from `pcluster` with the result from `swarm` that are both sorted alphabetically by header. They
are expected to be exactly the same.

##Run:##
`./main -t[threadcount] -z filename.fas`

Add this option `-y=[depth]` to enable economic search (skipping comparisons based on triangle inequality). 

Depth is the level of connections that will be using the skipping scheme. To skip only first level of connections use `-y=2`, while to skip up to sixth level of connections use `-y=7`. 
Based on our experiment input with short sequences (less than or equal to 125 nucleotides in average) do not benefit from multiple level connections when running clustering for `d>1`. Instead, the first level setting `-y=2` is advised to use for this type of short sequence input.

The program will run as brute force (comparing all possible combinations) when option `-y` is not specified, or if it is set as `-y=1`.

Other options for this program:

* `-d, --differences INTEGER          ` resolution (default 1)
* `-h, --help                         ` display this help and exit
* `-o, --output-file FILENAME         ` output result filename (default is `result.log`)
* `-t, --threads INTEGER              ` number of threads to use (default 1)
* `-m, --match-reward INTEGER         ` reward for nucleotide match (default 5)
* `-p, --mismatch-penalty INTEGER     ` penalty for nucleotide mismatch (default 4)
* `-g, --gap-opening-penalty INTEGER  ` gap open penalty (default 12)
* `-e, --gap-extension-penalty INTEGER` gap extension penalty (default 4)
* `-y, --depth INTEGER                ` maximum level for economic search, as explained before

##Some notes about preparing fasta file:##
Same as in swarm, fasta file needs to be sorted decreasing by abundance and marked with abundance information, e.g.

`>sequence-header_12345`


