#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void main() {
  double c;
  memset(&c, 0, sizeof(double));
  if (c == 0.0) {
    printf("True\n");
  }
  else {
    printf("False\n");
  }
}
