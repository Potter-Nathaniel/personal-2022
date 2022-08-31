all: linAlgTools

linAlgTools: driver.o matrix.o
	g++ -o linAlgTools driver.o matrix.o

driver.o: driver.cpp matrix.h
	g++ -c driver.cpp

matrix.o: matrix.cpp matrix.h
	g++ -c matrix.cpp

clean:
	rm *.o linAlgTools; clear

reset:
	make clean; make; clear;


