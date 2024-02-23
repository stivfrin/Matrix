#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = 0;

  if (rows <= 0 || columns <= 0) status = 1;

  if (status == OK) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));

    if (result->matrix)

      for (int i = 0; i < rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));

        if (!result->matrix[i]) s21_remove_matrix(result);
      }

    else
      status = CALCULATION_ERROR;
  }

  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) free(A->matrix[i]);

    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;

  if (A->rows == B->rows || A->columns == B->columns) {
    for (int i = 0; i < A->rows; ++i)

      for (int j = 0; j < A->columns; ++j)

        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-07) {
          status = FAILURE;
          break;
        }
  } else
    status = FAILURE;

  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = (correct_matrix(A) && correct_matrix(B));

  if (status == OK) {
    if (A->rows != B->rows || A->columns != B->columns)
      status = CALCULATION_ERROR;

    else {
      if (s21_create_matrix(A->rows, A->columns, result) != OK)
        status = CALCULATION_ERROR;

      else

        for (int i = 0; i < A->rows; i++)

          for (int j = 0; j < A->columns; j++)
            result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }

  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = (correct_matrix(A) && correct_matrix(B));

  if (status == OK) {
    if (A->rows != B->rows || A->columns != B->columns)
      status = CALCULATION_ERROR;

    if (s21_create_matrix(A->rows, A->columns, result) != OK)
      status = CALCULATION_ERROR;

    else

      for (int i = 0; i < A->rows; i++)

        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
  }

  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = correct_matrix(A);

  if (status == OK) {
    if (s21_create_matrix(A->rows, A->columns, result) != OK)
      status = CALCULATION_ERROR;

    else

      for (int i = 0; i < A->rows; i++)

        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = A->matrix[i][j] * number;
  }

  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = (correct_matrix(A) && correct_matrix(B));

  if (status == OK) {
    if (!(A && A->rows > 0 && A->columns > 0) ||
        !(B && B->rows > 0 && B->columns > 0))
      status = ERROR;

    if ((A->columns != B->rows) || (A->rows != B->columns))
      status = CALCULATION_ERROR;

    if (s21_create_matrix(A->rows, B->columns, result) != OK)
      status = CALCULATION_ERROR;

    else

      for (int i = 0; i < A->rows; i++)

        for (int j = 0; j < B->columns; j++) {
          result->matrix[i][j] = 0;
          for (int k = 0; k < A->columns; k++)
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
  }

  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = correct_matrix(A);

  if (status == OK) {
    if (s21_create_matrix(A->columns, A->rows, result) != OK)
      status = ERROR;

    else

      for (int i = 0; i < A->rows; i++)

        for (int j = 0; j < A->columns; j++)
          result->matrix[j][i] = A->matrix[i][j];
  }

  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = OK;

  if (!(A && A->rows > 0 && A->columns > 0) || !result) status = ERROR;

  if (status == OK) {
    if (A->rows != A->columns) status = CALCULATION_ERROR;

    if (s21_create_matrix(A->columns, A->rows, result) != OK)
      status = CALCULATION_ERROR;

    else
      complement(A, result);
  }
  return status;
}

int s21_determinant(matrix_t *A, double *result) {
  int status = OK;

  if (!(A && A->rows > 0 && A->columns > 0) || !result)
    status = ERROR;

  else {
    if (A->rows != A->columns) status = CALCULATION_ERROR;

    if (A->rows == 1)  // опред-ль = 1му эл. матрицы
      *result = A->matrix[0][0];

    else
      *result = determinant(A->matrix, A->rows);
  }

  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int status = correct_matrix(A);
  double det = 0;
  matrix_t add = {0};
  matrix_t transposed = {0};

  if (status == OK) {
    if (A->rows != A->columns)
      status = CALCULATION_ERROR;

    else {
      int code = s21_determinant(A, &det);

      if (fabs(det) < 1e-7 || code)
        status = CALCULATION_ERROR;
      else {
        if (s21_calc_complements(A, &add))
          status = CALCULATION_ERROR;

        else {
          if (s21_transpose(&add, &transposed) ||
              s21_create_matrix(A->rows, A->columns, result))
            status = CALCULATION_ERROR;
          else

            for (int i = 0; i < A->rows; i++)

              for (int j = 0; j < A->columns; j++)
                result->matrix[i][j] = transposed.matrix[i][j] / det;
        }
      }
    }
  }

  s21_remove_matrix(&transposed);
  s21_remove_matrix(&add);

  return status;
}