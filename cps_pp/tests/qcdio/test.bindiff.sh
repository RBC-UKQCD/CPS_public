#include<config.h>
CPS_START_NAMESPACE
for arg in D52C202K3500U007000T* 
do
echo "diffing $arg out.$arg"
diff "$arg" "out.$arg"
done



CPS_END_NAMESPACE
