TARGET = indexCalculus
LIBS = -lgmp -lm -lpthread 
CC = gcc
CFLAGS =  -Wall -Wextra -DB=0 -DP=\"666276812967623946997\" -DROOT=2 -DWANTED=42 -DVERBOSE=1  -O2 -o -lgmp -lpthread -lm

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
	clear