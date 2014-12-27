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
#mv debug.log debug1.log
#./main -t1 test.fas
mv result.log result1.log
#valgrind --dsymutil=yes --leak-check=yes --log-file=leak2 --show-possibly-lost=no ./main -t 2 test.fas
#mv result.log result2.log

#compareResult debug1.log testdata/data1 temp0
compareResult result1.log testdata/data2 temp1
#compareResult result2.log testdata/data2 temp2

