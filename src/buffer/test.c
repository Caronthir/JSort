#include "buffer.h"

int main(int argc, char *argv[])
{
    Buffer* buffer = makeBuffer("duck", "rb");
    printf("Good?: %i\n", buffer->good);
    if (buffer->good){
        bool success = read(buffer);
        printf("Success: %i\n", success);
        printf("Size : %zu\n", buffer->size);
        printf("Current : %zu\n", buffer->i);
    }
    unsigned int val;
    while(good(buffer)){
      val = next(buffer);
      printf("val: %u\n", val);
    }
    close(buffer);
    return 0;
}
