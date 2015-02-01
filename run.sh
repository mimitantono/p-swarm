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

rm *.log

make clean all

valgrind --tool=memcheck --dsymutil=yes --leak-check=yes --log-file=leak --show-possibly-lost=no ./test -t 1 test.fas
head -n 1476 result.log > result1.log
valgrind --dsymutil=yes --leak-check=yes --log-file=leak2 --show-possibly-lost=no ./test -t 2 test.fas
head -n 1476 result.log > result2.log

compareResult result1.log testdata/data2 temp1
compareResult result2.log testdata/data2 temp2

