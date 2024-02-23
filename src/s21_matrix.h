#ifndef SRC_S21_MATRIX_
#define SRC_S21_MATRIX_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS 1
#define FAILURE 0
#define ERROR 1

#define OK 0
#define INCORRECT_MATRIX 1
#define CALCULATION_ERROR 2

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

double determinant(double **A, int size);
void s21_sub_determinant(double **A, double **tmp, int skip_row, int skip_col,
                         int size);
void complement(matrix_t *A, matrix_t *result);
int correct_matrix(matrix_t *matrix);

#endif  // SRC_S21_MATRIX_
