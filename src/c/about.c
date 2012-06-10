#include <stdio.h>

int main(int argc, char *argv[]) {

  printf("%s\n", "==========================================================");
  printf("%s\n", "  Symmetrix (SMTRX) package");
  printf("%s\n", "  Fast matrix operations");
  printf("%s\n", "  See LICENSE.txt for license details.");
  printf("%s%s\n", "  Version: ", SYMTRX_VERSION);
  printf("%s%s\n", "  Build: ", SYMTRX_BUILD);
  printf("%s\n", "==========================================================");

  return 0;

}
