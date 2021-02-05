CC=icpc -fopenmp -O3
CFLAGS=-I${MKLROOT}/include
LDFLAGS=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl


main: main.o ED.o
	$(CC) $^ ${LDFLAGS} -o $@

main.o: main.cpp ED.h
	$(CC) -c $<

ED.o: ED.cpp ED.h
	$(CC) $(CFLAGS) -c $<

.PHONY:clean
clean:
	rm -rf *.o main