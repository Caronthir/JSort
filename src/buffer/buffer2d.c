#include "buffer2d.h"

Buffer2D* makeBuffer2D(const char *fname, const char *mode){
  Buffer2D *buffer = malloc(sizeof(Buffer2D));
  initBuffer(buffer);
  open(buffer, fname, mode);
  return buffer;
}

Buffer2D *makeBufferW2D(const char *fname) {
    Buffer2D *buffer = makeBuffer2D(fname, "wb");
    return buffer;
}

void close2D(Buffer2D* buffer){
    write2D(buffer);
    free(buffer);
}

bool write2D(Buffer2D* buffer){
  if (buffer->i > buffer->size) {
    printf("Buffer pointer is outside of buffer size");
    return false;
  }
  fwrite(buffer->vals, sizeof(unsigned int), 2*buffer->i, buffer->handle);
  //fflush(buffer->handle);
  buffer->i = 0;
  return true;
}

bool conditionalWrite2D(Buffer2D *buffer) {
  if (buffer->i >= buffer->size) {
    return write2D(buffer);
  } else {
    return false;
  }
}

void fill2D(Buffer2D* buffer, unsigned int x, unsigned int y){
    buffer->vals[buffer->i][0] = x;
    buffer->vals[buffer->i][1] = y;
    increment(buffer);
    conditionalWrite2D(buffer);
}

