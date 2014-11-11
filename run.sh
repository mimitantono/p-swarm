#!/bin/bash

make clean
make all
./main -t 1 -o run.log test.fas
