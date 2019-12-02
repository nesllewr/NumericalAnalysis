CC = gcc
CFLAGS = -W -Wall -g
TARGET = na_hw9

$(TARGET) : na_hw9.o nrutil.o gaussj.o
			 $(CC) $(CFLAGS) -o $@ $^ -lm

na_hw9.o : na_hw9.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

nrutil.o : nrutil.c
		$(CC) $(CFLAGS) -c -o $@ $^

gaussj.o : gaussj.c
		$(CC) $(CFLAGS) -c -o $@ $^


clean :
		rm *.o 
		rm $(TARGET)
