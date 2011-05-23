CFLAGS=-std=c99

mazing: main.c mazing.o
	$(CC) -o mazing $^ $(LIBS)
