CC = gcc
CFLAGS = -W -Wall -g
TARGET = na_hw4

$(TARGET) : na_hw4.o nrutil.o rtbis.o rtflsp.o rtnewt.o rtsafe.o rtsec.o
			 $(CC) $(CFLAGS) -o $@ $^ -lm

na_hw4.o : na_hw4.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

nrutil.o : nrutil.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtbis.o : rtbis.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtflsp.o : rtflsp.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtnewt.o : rtnewt.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtsafe.o : rtsafe.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtsec.o : rtsec.c
		$(CC) $(CFLAGS) -c -o $@ $^


clean :
		rm *.o 
		rm $(TARGET)
