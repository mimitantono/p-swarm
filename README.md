#pcluster#

Parallel hierarchical single-linkage clustering for DNA sequences, adapted from **swarm**.

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


##Run:##
`./main -t[threadcount] -z filename.fas`

`-z` to enable speed improvement with alternative algorithm (skipping comparisons based on triangle inequality)
without -z program will run as brute force (comparing all possible combinations)

##Some notes about preparing fasta file:##
Same as in swarm, fasta file needs to be presorted and marked with abundance information.

Add these lines to `~/.bash_profile` and then restart shell to use : `prepareFasta filename.fas`

```
function prepareFasta() {
   awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' $1 > temp
   grep -v "^>" temp | \
   grep -v [^ACGTacgt] | sort -d | uniq -c | \
   while read abundance sequence ; do
     hash=$(printf "${sequence}" | shasum)
     hash=${hash:0:40}
     printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
   done | sort -t "_" -k2,2nr -k1.2,1d | \
   sed -e 's/\_/\n/2' > $2
}
```
