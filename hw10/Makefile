CC = gcc
CFLAGS = -W -Wall -g
TARGET = na_hw10

$(TARGET) : na_hw10.o nrutil.o gaussj.o mrqmin.o gasdev.o jacobi.o fgauss.o ran1.o mrqcof.o covsrt.o
			 $(CC) $(CFLAGS) -o $@ $^ -lm

na_hw10.o : na_hw10.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

nrutil.o : nrutil.c
		$(CC) $(CFLAGS) -c -o $@ $^

gaussj.o : gaussj.c
		$(CC) $(CFLAGS) -c -o $@ $^

mrqmin.o : mrqmin.c
		$(CC) $(CFLAGS) -c -o $@ $^

gasdev.o : gasdev.c
		$(CC) $(CFLAGS) -c -o $@ $^

jacobi.o : jacobi.c
		$(CC) $(CFLAGS) -c -o $@ $^

mrqcof.o : mrqcof.c
		$(CC) $(CFLAGS) -c -o $@ $^

covsrt.o : covsrt.c
		$(CC) $(CFLAGS) -c -o $@ $^

ran1.o : ran1.c
		$(CC) $(CFLAGS) -c -o $@ $^

fgauss.o : fgauss.c
		$(CC) $(CFLAGS) -c -o $@ $^



clean :
		rm *.o 
		rm $(TARGET)
