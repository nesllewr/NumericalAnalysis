CC = gcc
CFLAGS = -W -Wall -g
TARGET = na_hw6

$(TARGET) : na_hw6.o nrutil.o ran1.o gasdev.o
			 $(CC) $(CFLAGS) -o $@ $^ -lm

na_hw6.o : na_hw6.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

nrutil.o : nrutil.c
		$(CC) $(CFLAGS) -c -o $@ $^

ran1.o : ran1.c
		$(CC) $(CFLAGS) -c -o $@ $^

gasdev.o : gasdev.c
		$(CC) $(CFLAGS) -c -o $@ $^

clean :
		rm *.o 
		rm $(TARGET)
