CC = gcc
CFLAGS = -W -Wall -g
TARGET = na_hw7

$(TARGET) : na_hw7.o nrutil.o ran1.o gasdev.o jacobi.o eigsrt.o 
			 $(CC) $(CFLAGS) -o $@ $^ -lm

na_hw7.o : na_hw7.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

nrutil.o : nrutil.c
		$(CC) $(CFLAGS) -c -o $@ $^

ran1.o : ran1.c
		$(CC) $(CFLAGS) -c -o $@ $^

gasdev.o : gasdev.c
		$(CC) $(CFLAGS) -c -o $@ $^

jacobi.o : jacobi.c
		$(CC) $(CFLAGS) -c -o $@ $^

eigsrt.o : eigsrt.c
		$(CC) $(CFLAGS) -c -o $@ $^

clean :
		rm *.o 
		rm $(TARGET)
