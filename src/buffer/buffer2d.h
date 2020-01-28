#ifndef BUFFER2D_H
#define BUFFER2D_H

#include <bits/types/FILE.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "buffer.h"


typedef struct Buffer2D {
  size_t i;
  size_t size;
  bool good;
  FILE *handle;
  unsigned int vals[BUFFER_SIZE][2];
} Buffer2D;

Buffer2D *makeBuffer2D(const char *fname, const char *mode);
Buffer2D *makeBufferW2D(const char *fname);
void fill2D(Buffer2D*, unsigned int x, unsigned int y);
bool write2D(Buffer2D *);
bool conditinalWrite2D(Buffer2D *);
//bool conditionalRead(Buffer *);
//bool read(Buffer2D *);
void close2D(Buffer2D *);
//unsigned int next(Buffer *);

#endif /* BUFFER2D_H */
