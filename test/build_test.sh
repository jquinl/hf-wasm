clang -Os -Wall -D MY_TESTS=1 -I.  -L. -lm -o out/test.o  ../hf/functions.c ../hf/primitives.c ../hf/integrals_naive.c ../hf/rhf.c ../hf/utils.c tests_hf.c
#clang -Os -Wall -D MY_TESTS=1 -I.  -L. -lm -o out/test.o  ../hf/functions.c ../hf/primitives.c ../hf/integrals_naive.c ../hf/rhf.c ../hf/utils.c tests_fn.c

#gcc -O0 -D MY_TESTS=1  -I.  -L. ../hf/functions.c ../hf/primitives.c ../hf/integrals_naive.c ../hf/rhf.c ../hf/utils.c tests.c -lm  -o out/test.o  