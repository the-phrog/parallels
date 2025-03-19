CC=g++
SRC=main.cpp
TARGET=sin_array_sum

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) -o $(TARGET) $(SRC)

float:
	$(CC) -o $(TARGET) $(SRC)

double:
	$(CC) -o $(TARGET) $(SRC) -DUSE_DOUBLE

clean:
	rm -f $(TARGET)
