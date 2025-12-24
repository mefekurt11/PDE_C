# PDE Edge-Based Image Compression Codec
# Based on Mainberger et al. paper

# Compiler (use gcc for C)
CC = gcc

# Compiler flags
# -I tells the compiler WHERE to find jbig.h (Homebrew path on Apple Silicon)
CFLAGS = -O3 -std=c11 -I/opt/homebrew/include -Wall

# Linker flags
# -L tells the linker WHERE to find libjbig.a
# -ljbig links the JBIG library
# -lm links the math library
LDFLAGS = -L/opt/homebrew/lib -ljbig -lm

TARGET = pde_codec
SRCS = main.c pde_codec.c

all: $(TARGET)

$(TARGET): $(SRCS) pde_codec.h
	$(CC) $(CFLAGS) $(SRCS) $(LDFLAGS) -o $(TARGET)

clean:
	rm -f $(TARGET) compressed.pde output_color.png

.PHONY: all clean