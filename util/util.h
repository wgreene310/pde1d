#include <Eigen/Core>

#include <nvector/nvector_serial.h>

void print(N_Vector v, const char *title);

Eigen::VectorXd linspace(double a, double b, int n);

template<class T>
void printMat(const T &a, const char *title, const char *format = "%16.9e,") {
  const int m = a.rows();
  const int n = a.cols();
  printf("%s(%d,%d)\n", title, m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf(format, a(i, j));
    }
    printf("\n");
  }
  printf("\n");
}