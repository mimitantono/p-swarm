#!/bin/bash

function compareResult() {
  echo ""
  echo   "Comparing $1 with $2:"
  if diff $1 $2 > $3; then
    echo "**** OK! ****";
    rm $3;
  else
    echo "===== test failed ====";
    echo "      check $3 file for diff details";
  fi  
}

make clean
rm *.log

make all

valgrind --dsymutil=yes --leak-check=yes --log-file=leak --show-possibly-lost=no ./main -t 1 test.fas
mv debug.log debug1.log
mv result.log result1.log
valgrind --dsymutil=yes --leak-check=yes --log-file=leak2 --show-possibly-lost=no ./main -t 2 test.fas
mv result.log result2.log
./main -t 3 test.fas
mv result.log result3.log
./main -t 4 test.fas
mv result.log result4.log
./main -t 5 test.fas
mv result.log result5.log
./main -t 6 test.fas
mv result.log result6.log
./main -t 7 test.fas
mv result.log result7.log
./main -t 8 test.fas
mv result.log result8.log


compareResult debug1.log testdata/data1 temp0
compareResult result1.log testdata/data2 temp1
compareResult result2.log testdata/data2 temp2
compareResult result3.log testdata/data2 temp3
compareResult result4.log testdata/data2 temp4
compareResult result5.log testdata/data2 temp5
compareResult result6.log testdata/data2 temp6
compareResult result7.log testdata/data2 temp7
compareResult result8.log testdata/data2 temp8


