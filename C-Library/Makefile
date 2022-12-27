.PHONY: all clean

CC=gcc
CFLAGS=-Wall -Wextra -Werror
LDFLAGS=$(shell pkg-config --cflags --libs check)
GCOVFLAGS=-fprofile-arcs -ftest-coverage

all: s21_matrix.a

s21_matrix.a:
	$(CC) $(CFLAGS) -c s21_matrix.c
	ar rc s21_matrix.a *.o
	ranlib s21_matrix.a

clean:
	rm -rf *.o *.html *.gcda *.gcno *.css *.a *.gcov *.info *.out *.cfg *.txt lib_tests

test:
	$(CC) $(CFLAGS) $(LDFLAGS) s21_matrix.c tests/*.c -o lib_tests
	./lib_tests
	
gcov_report:
	$(CC) $(CFLAGS) $(LDFLAGS) $(GCOVFLAGS) s21_matrix.c tests/*.c -o gcov_main
	./gcov_main
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory gcov_report
