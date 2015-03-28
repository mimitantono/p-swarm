#!/bin/bash

function compareResult() {
  echo ""
#  echo   "Comparing $1 with $2:"
  if diff $1 $2 > $3; then
    echo "OK    : $4";
    rm $3;
  else
    echo "FAIL  : $4";
    echo "      check $3 file for diff details";
  fi  
}

rm *.log

make all

./test -t 1 test.fas
head -n 1476 result.log > result1.log
./test -t 2 test.fas
head -n 1476 result.log > result2.log
./test -t 1  test.fas
head -n 1476 result.log > result3.log
./test -t 2  test.fas
head -n 1476 result.log > result4.log
./test -t 8  test.fas
head -n 1476 result.log > result5.log
./test -t 8 -y 7  test.fas
head -n 1476 result.log > result6.log
./test -t 1 -d 2 test.fas
head -n 1441 result.log > result7.log
./test -t 8 -d 2 -y 7  test.fas
head -n 1441 result.log > result8.log
./test -t 8 -d 8 -y 4  test.fas
head -n 1414 result.log > result9.log

compareResult result1.log testdata/data2 temp1 "./test -t 1 test.fas"
compareResult result2.log testdata/data2 temp2 "./test -t 2 test.fas"
compareResult result3.log testdata/data2 temp3 "./test -t 1  test.fas"
compareResult result4.log testdata/data2 temp4 "./test -t 2  test.fas"
compareResult result5.log testdata/data2 temp5 "./test -t 8  test.fas"
compareResult result6.log testdata/data2 temp6 "./test -t 8 -y 7  test.fas"
compareResult result7.log testdata/data3 temp7 "./test -t 1 -d 2 test.fas"
compareResult result8.log testdata/data3 temp8 "./test -t 8 -d 2 -y 7  test.fas"
compareResult result9.log testdata/data4 temp9 "./test -t 8 -d 8 -y 4  test.fas"
