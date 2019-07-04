CC = g++
CFLAGS = -Wall -O3 --std=c++11
esegui : main.exe
	./main.exe
main.exe : main.o random.o Buffon.o Statistica.o
	$(CC) random.o main.o Buffon.o Statistica.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
Buffon.o : Buffon.cpp Buffon.h random.h
	$(CC) -c Buffon.cpp -o Buffon.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Statistica.o : Statistica.cpp Statistica.h
	$(CC) -c Statistica.cpp -o Statistica.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out 
