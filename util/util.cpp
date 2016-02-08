#include <stdio.h>

#include "util.h"

void print(N_Vector v, const char *title) {
  double *d = NV_DATA_S(v);
  int len = NV_LENGTH_S(v);
  int count = 0;
  printf("Vector: %s(%d)\n", title, len);
  for (int i = 0; i < len; i++) {
    printf(" %14.9e", d[i]);
    if (++count == 6) {
      printf("\n");
      count = 0;
    }
  }
  if (count)
    printf("\n");
}