CC=gcc
CFLAGS=-Wall -Werror -Wextra -std=c11
GCOV_FLAGS=-fprofile-arcs -ftest-coverage -lgcov
SOURCE=s21_matrix.c s21_helpers.c
OBJ = $(SOURCE:.c=.o) s21_matrix_tests.c
OBJST = $(SOURCE:.c=.o)

ifeq (&(OS), Linux)
		CHECK_FLAGS = -lcheck -lm -pthread -lrt -lsubunit
else
		CHECK_FLAGS = -lcheck
endif

all: s21_matrix.a test gcov_report

s21_matrix.a:
	@$(CC) $(CFLAGS) -c $(SOURCE)
	@ar rc s21_matrix.a $(OBJST)
	@rm -rf *.o
test: $(OBJ)
	@$(CC) $(CHECK_FLAGS) $(OBJ) -g -o test
	@rm -rf *.o
	@./test
	@rm test
gcov_report:
	@$(CC) -o test $(SOURCE) s21_matrix_tests.c $(GCOV_FLAGS) $(CHECK_FLAGS)
	@./test
	@lcov -t "gcov_report" -o s21_matrix_tests.info -c -d .
	@genhtml -o report s21_matrix_tests.info
	@rm -rf *.gcno
	@rm -rf *.gcda
	@open report/src/index.html
	@rm test
clean:
	@rm -rf *.info
	@rm -rf *.a
	@rm -rf test.dSYM
	@rm -rf report
	@rm test
rebuild: clean all