CC = gcc
CFLAGS = -fPIC -c -Wall -O3
LDFLAGS = -shared
TARGET = libs/liballfuncs_FF_0.04_rep_1.dll

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

# Compile dispatche
dispatcher.o: dispatcher.c
	$(CC) $(CFLAGS) $< -o $@

# Compile tous les fichiers C dans Cfunctions/
Cfunctions/%.o: Cfunctions/%.c
	$(CC) $(CFLAGS) $< -o $@

# Link final après compilation de tous les objets
$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ @objs.txt

clean:
	rm -f $(OBJ) $(TARGET)