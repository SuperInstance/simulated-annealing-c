CC = gcc
CFLAGS = -O2 -Wall -Wextra -lm

test: test.c sa.h
	$(CC) $(CFLAGS) -o test_sa test.c -lm
	./test_sa

clean:
	rm -f test_sa

.PHONY: test clean
