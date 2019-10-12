CC = gcc
CFLAGS = -W -Wall -g
TARGET = na_hw5

$(TARGET) : na_hw5.o nrutil.o gaussj.o ludcmp.o svdcmp.o svbksb.o mprove.o pythag.o
			 $(CC) $(CFLAGS) -o $@ $^ -lm

na_hw5.o : na_hw5.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

nrutil.o : nrutil.c
		$(CC) $(CFLAGS) -c -o $@ $^

gaussj.o : gaussj.c
		$(CC) $(CFLAGS) -c -o $@ $^

ludcmp.o : ludcmp.c
		$(CC) $(CFLAGS) -c -o $@ $^

svdcmp.o : svdcmp.c
		$(CC) $(CFLAGS) -c -o $@ $^

svbksb.o : svbksb.c
		$(CC) $(CFLAGS) -c -o $@ $^

mprove.o : mprove.c
		$(CC) $(CFLAGS) -c -o $@ $^

pythag.o : pythag.c
		$(CC) $(CFLAGS) -c -o $@ $^

clean :
		rm *.o 
		rm $(TARGET)
