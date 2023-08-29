# clonesim

## compilation

```
mkdir build
cd build
cmake ../ -DLIBLEMON_ROOT=/scratch/software/lemon/  -DBOOST_ROOT=/scratch/software/boost_1_74_0/
make
```

## usage

```
./generatecnatrees > cnatrees.txt
./simulate -S cnatrees.txt > T.dot
dot -Tpdf T.dot > T.pdf
```