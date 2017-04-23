#include <Eigen/Core>

#include <nvector/nvector_serial.h>

void print(N_Vector v, const char *title);

Eigen::VectorXd linspace(double a, double b, int n);

template<class T>
void printMat(const T &a, const char *title, const char *format = "%16.9e,") {
  const size_t m = a.rows();
  const size_t n = a.cols();
  printf("%s(%zd,%zd)\n", title, m, n);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      printf(format, a(i, j));
    }
    printf("\n");
  }
  printf("\n");
}