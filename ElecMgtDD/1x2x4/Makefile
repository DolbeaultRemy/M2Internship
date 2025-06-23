CC = gcc
CFLAGS = -fPIC -c -Wall -O3
LDFLAGS = -shared
TARGET = liballfuncs.dll

# Sources
DISPATCHER = dispatcher.c
CFUNCTIONS_DIR = Cfunctions
CFUNCTIONS = $(wildcard $(CFUNCTIONS_DIR)/*.c)

# Objects
DISPATCHER_OBJ = dispatcher.o
CFUNCTIONS_OBJ = $(CFUNCTIONS:.c=.o)
OBJ = $(DISPATCHER_OBJ) $(CFUNCTIONS_OBJ)

.PHONY: all clean

all: $(TARGET)

# Compile the dispatcher
dispatcher.o: dispatcher.c
	$(CC) $(CFLAGS) $< -o $@

# Compile all the C files in Cfunctions/
Cfunctions/%.o: Cfunctions/%.c
	$(CC) $(CFLAGS) $< -o $@

# Final link after the compilation of all the objects
$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ @objs.txt

clean:
	rm -f $(OBJ) $(TARGET)
