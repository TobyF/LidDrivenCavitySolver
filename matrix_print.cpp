#include "matrix_print.h"
// Print a matrix
void printMatrixCM(double *F, int nx, int ny) {
      for (int j = ny-1; j > -1; --j) {
        for (int i = 0; i < nx; ++i) {
            cout << setw(10) << F[i*ny+j] << " ";
        }
        cout << endl;
    }
}

void printMatrixRM(double *F, int nx, int ny) {
for (int j = 0; j < ny; ++j){
  for (int i = 0; i < nx; ++i) {
            cout << setw(5) << F[i+j*nx] << " ";
        }
        cout << endl;
    }
}
