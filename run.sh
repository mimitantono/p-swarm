#!/bin/bash

make clean
make all
./main -t 1 test.fas
diff debug.log testdata/test.log

