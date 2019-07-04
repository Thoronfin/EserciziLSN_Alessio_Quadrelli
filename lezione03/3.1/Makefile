CC = g++
CFLAGS = -Wall -O3 --std=c++11
esegui : main.exe
	./main.exe
main.exe : main.o random.o OptionPrice.o Statistica.o
	$(CC) random.o main.o OptionPrice.o Statistica.o -o main.exe
main.o : main.cpp random.h OptionPrice.h Statistica.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h 
	$(CC) -c random.cpp -o random.o $(CFLAGS)
OptionPrice.o : OptionPrice.cpp OptionPrice.h random.h
	$(CC) -c OptionPrice.cpp -o OptionPrice.o $(CFLAGS)
Statistica.o : Statistica.cpp Statistica.h 
	$(CC) -c Statistica.cpp -o Statistica.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

