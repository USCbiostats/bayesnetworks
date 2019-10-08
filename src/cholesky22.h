// Cholesky Decomposition of a positive definite matrix
// Input values
//    double **x - matrix to be decomposed, assumed to be
//               - symmetric. Only the upper triangular part
//               - is used. This matrix is usually declared
//				 - dynamically.
//    int n      - Dimension of the matrix
//    double **c - Destination matrix. This must already be
//               - declared by the user. This matrix should
//               - also be declared dynamically.
// Return value
//    0: No error
//    1: No input matrix passed
//    2: Dimension of matrix less than 1
//    3: No output matrix specified
//    4: Input matrix not positive definite
//    5: Unknown error (array out of bounds error most likely)

// Note: If using this to sample from a multivariate normal,
//       the matrix passed should be the variance-covariance
//       matrix. The matrix must be of full rank. I.e. knowing
//       the value of any number of the variables cannot uniquely
//       specify the value of any of the remaining variables.

int CholeskyDecomp(double **x, int n, double **c) {
  int i, j, k;
  double sum;

  // Verify input values are valid
  if (x == NULL)
    return 1;
  if (n < 1)
    return 2;
  if (c == NULL)
    return 3;

  try {
    for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
        c[i][j] = x[i][j];
      }
    }
    for (i = 0; i < n; ++i) {
      for (j = i; j < n; ++j) {
        sum = x[i][j];
        for (k = i - 1; k > -1; --k) {
          sum -= c[i][k] * c[j][k];
        }
        if (i == j) {
          // If sum <= 0 then the matrix is not
          // postive definite
          if (sum <= 0.)
            return 4;
          c[i][i] = sqrt(sum);
        } else {
          c[j][i] = sum / c[i][i];
        }
      }
    }
  } catch(...) {
    // An error occurred trying to read the input matrix
    // or writing to the output matrix.
    return 5;
  }
  return 0;
}

// Invert a positive definite symmetric (PDS) matrix
// using lower-upper decomposition with forward and
// back substituion. Decomposition is performed using
// the Cholesky algorithm for positive definite
// symmetric matrices.
//
// Input values
//		double **x - Matrix to invert. The matrix is assumed
//					 to be symmetric and should be declared
//					 dynamically.
//		int n	   - Dimension of the matrix.
//		double **c - Array to hold the inverse matrix. Must
//					 be declared dynamically by user.
//
// Return value
//		0 - No errors
//		1 - No variance-covariance matrix supplied
//		2 - Dimension of matrix less than 1
//		3 - No output array specified
//      4 - Unkown error - Input array overflow most likely
//      5 - Unkown error - Output array overflow most likely
//    >10 - Error in decomposition. Subtract 10 from value
//			to get error number from CholeskyDecomp function
//
int InvertPDS(double **x, int n, double **c) {
  int i, j, k, rc;
  double **s;
  double sum;

  if (x == NULL)
    return 1;
  if (n < 1)
    return 2;
  if (c == NULL)
    return 3;

  rc = 0;
  // Declare matrix to hold lower-upper decompostion
  s = new double *[n];
  for (i = 0; i < n; ++i)
    s[i] = new double[n];

  // Set solution to identity matrix
  try {
    for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
        if (i == j)
          c[i][j] = 1.;
        else
          c[i][j] = 0.;
      }
    }
  } catch (...) {
    // Failure writing to output matrix
    rc = 4;
  }
  // Lower-Upper decomposition using Cholesky algorithm
  if (rc == 0) {
    rc = CholeskyDecomp(x, n, s);
    if (rc == 0) {
      // Solve for each row using back and forward
      // substitution
      try {
        for (i = 0; i < n; ++i) {
          // Forward substitution
          for (j = 0; j < n; ++j) {
            sum = c[i][j];
            for (k = 0; k < j; ++k) {
              sum -= s[j][k] * c[i][k];
            }
            c[i][j] = sum / s[j][j];
          }
          // Back substitution
          for (j = n - 1; j > -1; --j) {
            sum = c[i][j];
            for (k = n - 1; k > j; --k) {
              sum -= s[k][j] * c[i][k];
            }
            c[i][j] = sum / s[j][j];
          }
        }
      } catch (...) {
        // Failure writing to output matrix
        // Should have been caught earlier
        rc = 4;
      }
    } else {
      // Failure in decompostion
      rc += 10;
    }
  }

  // delete lower-upper matrix decomposition
  if (s != NULL) {
    for (i = 0; i < n; ++i) {
      if (s[i] != NULL)
        delete [] s[i];
    }
    delete [] s;
  }

  return rc;
}

// Invert a positive definite symmetric (PDS) matrix
// using lower-upper decomposition with forward and
// back substituion. Decomposition is performed using
// the Cholesky algorithm for positive definite
// symmetric matrices.
//
// Input values
//		double *x  - Matrix to invert. The matrix is assumed
//					 to be symmetric and should be declared
//					 statically.
//		int n	   - Dimension of the matrix.
//		double *c  - Array to hold the inverse matrix. Must
//					 be declared statically.
//
//		Note: Code to properly use this function looks like
//            this.
//					double x[n][n], c[n][n];
//                    . . .
//                  InvertPDS(x[0], n, c[0]);
//
// Return value
//		0 - No errors
//		1 - No variance-covariance matrix supplied
//		2 - Dimension of matrix less than 1
//		3 - No output array specified
//      4 - Unkown error - Input array overflow most likely
//      5 - Unkown error - Output array overflow most likely
//    >10 - Error in decomposition. Subtract 10 from value
//			to get error number from CholeskyDecomp function
//
int InvertPDS(double *x, int n, double *c) {
  int i, rc;
  double **x2, **c2;

  if (x == NULL)
    return 1;
  if (n < 1)
    return 2;
  if (c == NULL)
    return 1;

  rc = 0;
  x2 = new double *[n];
  c2 = new double *[n];


  // Create pointers to rows of matrix
  // Equivalent to dynamic allocation
  try {
    x2[0] = x;
    c2[0] = c;
    for (i = 1; i < n; ++i) {
      x2[i] = x2[i - 1] + n;
      c2[i] = c2[i - 1] + n;
    }
  } catch (...) {
    // Error creating pointers to rows
    rc = 4;
  }

  if (rc == 0)
    // Call routine using dynamically allocated matrices
    rc = InvertPDS(x2, n, c2);

  // Clean up memory
  if (x2 != NULL)
    delete [] x2;
  if (c2 != NULL)
    delete [] c2;
  return rc;
}
