CFLAGS=-Wall -std=c99 -O2
LIBS=-lgmp

mazing: main.c fmc.o mazing.o
	$(CC) $(CFLAGS) -o mazing $^ $(LIBS)

.PHONY: clean
clean:
	rm -fr *.o mazing *.dSYM
