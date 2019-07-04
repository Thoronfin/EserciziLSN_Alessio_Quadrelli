CC = g++
CFLAGS = -Wall -O3 --std=c++11

esegui : main.exe
	./main.exe
main.exe : main.o random.o Rw.o Statistica.o
	$(CC) random.o main.o Rw.o Statistica.o -o main.exe
main.o : main.cpp random.h Rw.h Statistica.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h 
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Rw.o : Rw.cpp Rw.h random.h
	$(CC) -c Rw.cpp -o Rw.o $(CFLAGS)
Statistica.o : Statistica.cpp Statistica.h 
	$(CC) -c Statistica.cpp -o Statistica.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

