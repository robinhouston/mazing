CFLAGS=-Wall -std=c99 -g
LIBS=-lgmp

mazing: main.c fmc.o mazing.o
	$(CC) $(CFLAGS) -o mazing $^ $(LIBS)

.PHONY: clean
clean:
	rm -fr *.o mazing *.dSYM
