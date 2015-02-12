#!/bin/bash

function compareResult() {
  echo ""
#  echo   "Comparing $1 with $2:"
  if diff $1 $2 > $3; then
    echo "OK    : **** $4 ****";
    rm $3;
  else
    echo "FAIL  : ==== $4 ====";
    echo "      check $3 file for diff details";
  fi  
}

rm *.log

make all

./test -t 1 test.fas
head -n 1476 result.log > result1.log
./test -t 2 test.fas
head -n 1476 result.log > result2.log
./test -t 1 -z test.fas
head -n 1476 result.log > result3.log
./test -t 2 -z test.fas
head -n 1476 result.log > result4.log
./test -t 8 -z test.fas
head -n 1476 result.log > result5.log
./test -t 8 -y 4 -z test.fas
head -n 1476 result.log > result6.log

compareResult result1.log testdata/data2 temp1 "Test single thread disable flag"
compareResult result2.log testdata/data2 temp2 "Test 2 threads disable flag"
compareResult result3.log testdata/data2 temp3 "Test single thread enable flag depth default (2)"
compareResult result4.log testdata/data2 temp4 "Test 2 threads enable flag depth default (2)"
compareResult result5.log testdata/data2 temp5 "Test 8 threads enable flag depth default (2)"
compareResult result6.log testdata/data2 temp6 "Test 8 threads enable flag depth 4"
