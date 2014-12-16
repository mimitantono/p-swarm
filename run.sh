#!/bin/bash

function compareResult() {
  echo ""
  echo "Comparing $1 with $2:"
  if diff $1 $2 > $3; then
    echo "**** OK! ****";
    rm $3;
  else
    echo "===== test failed ====";
    echo "check $3 file for diff details";
  fi  
}

make clean
rm *.log

make all
#./main -t 1 test.fas
valgrind --dsymutil=yes --leak-check=yes --log-file=leak --show-possibly-lost=no ./main -t 1 test.fas
mv debug.log debug1.log
mv result.log result1.log
#./main -t 2 test.fas
#mv debug.log debug2.log
#mv result.log result2.log
#./main -t 8 test.fas
#mv debug.log debug3.log
#mv result.log result3.log

compareResult debug1.log testdata/data1 temp1
compareResult result1.log testdata/data2 temp11
#compareResult result2.log testdata/data2 temp2
#compareResult result3.log testdata/data2 temp3

