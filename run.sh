#!/bin/bash

make clean
rm *.log

make all
#./main -t 1 test.fas
valgrind --dsymutil=yes --leak-check=yes --log-file=leak --show-leak-kinds=definite ./main -t 1 test.fas

echo ""
echo "Comparing results with oracle:"
if diff debug.log testdata/test.log >temp; then
    echo "**** OK! ****";
    rm temp;
else
    echo "Oh Nooooooo... Failed :(";
    echo "check temp file for diff details";
fi

