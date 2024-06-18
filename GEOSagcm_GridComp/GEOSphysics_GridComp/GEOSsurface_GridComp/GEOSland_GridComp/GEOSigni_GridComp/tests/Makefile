all: cffwi-test
	./cffwi-test

cffwi-test:
	gfortran -g -Wall -o cffwi-test ../cffwi.F90 ./cffwi-test.F90

clean:
	rm ./cffwi.mod
	rm ./cffwi-test
