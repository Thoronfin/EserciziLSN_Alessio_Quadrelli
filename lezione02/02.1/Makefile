CC = g++
CFLAGS = -Wall -O3 --std=c++11

esegui : main.exe
	./main.exe
main.exe : main.o random.o Integrator.o Funzioni.o Statistica.o
	$(CC) random.o main.o Integrator.o Funzioni.o Statistica.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
Integrator.o : Integrator.cpp Integrator.h Funzioni.h random.h
	$(CC) -c Integrator.cpp -o Integrator.o $(CFLAGS)
Funzioni.o : Funzioni.cpp Funzioni.h
	$(CC) -c Funzioni.cpp -o Funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Statistica.o : Statistica.cpp Statistica.h
	$(CC) -c Statistica.cpp -o Statistica.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out 
