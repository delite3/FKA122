CC = gcc
CFLAGS = -Wall -O2
SRCDIR = src
TARGET = metro

all: $(TARGET)

$(TARGET): $(SRCDIR)/metro.c
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCDIR)/metro.c

clean:
	rm -f $(TARGET)
