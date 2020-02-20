#include "reader.h"
#include <iostream>
#include <string.h>

// assumes little endian
void printBits(size_t const size, void const *const ptr) {
  unsigned char *b = (unsigned char *)ptr;
  unsigned char byte;
  int i, j;

  for (i = size - 1; i >= 0; i--) {
    for (j = 7; j >= 0; j--) {
      byte = (b[i] >> j) & 1;
      printf("%u", byte);
    }
  }
  puts("");
}
#define BITS(x) printBits(sizeof(x), &x)

int main(int argc, char *argv[])
{
    // unsigned int j = encode(5, 6, 21);
    // printf("%u\n", makeMask(0, 16));
    // printf("%u\n", makeMask(16, 24));
    // printf("%u\n", makeMask(24, 32));
    // // BITS(j);
    // printf("%u\n", j);
    // unsigned int front, back, size;
    // decode(j, &front, &back, &size);
    // printf("Front: %u\tBack: %u\tLaBr: %u\n", front, back, size);
    if (argc != 3){
        std::cout << "Usage: [zinc/si] <output path>\n";
        return -1;
    }
    std::string zinc = "/home/erdos/master/sortering/datafiles/sirius-20180125-101223.data";
    std::string si = "/home/erdos/master/sortering/datafiles/sirius-20180125-091200.data";
    std::string input;
    if (strcmp(argv[1], "zinc") == 0){
        input = zinc;
    } else if (strcmp(argv[1], "si") == 0){
        input = si;
    } else {
        return -1;
    }
    Reader reader(input, argv[2]);
    reader.read();
    reader.summary();
    return 0;
}
