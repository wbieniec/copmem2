echo "#############################################################################"
echo "# Please compile on Linux, using make all command                            "
echo "# Please provide valid dataset files                                         "
echo "#############################################################################"

#!/bin/bash
MemFill="./memoryfill 118G"

R="../datasets/hum.all.fa"
Q="../datasets/panTro3.fa"
TST="hp"

for pass in a b c
do
for L in 200 100 80
  do
    for T in 8 4 1
    do
      echo "Test =" $TST "  Pass =" $pass "  L =" $L "  Threads =" $T
      eval $MemFill
      /usr/bin/time --verbose ./copmem2 -t $T -l $L -o $TST-$L-$T-$pass.txt $R $Q
      rm $TST-$L-$T-$pass.txt

      eval $MemFill
      /usr/bin/time --verbose ./copmem2 -mf -t $T -l $L -o $TST-$L-$T-$pass.txt $R $Q
      rm $TST-$L-$T-$pass.txt
    done
  done
done