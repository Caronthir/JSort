#ifndef BUFFER_H
#define BUFFER_H

#include <bits/types/FILE.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define BUFFER_SIZE 65536/2

typedef struct Buffer
{
  size_t i;
  size_t size;
  bool good;
  FILE *handle;
  unsigned int vals[BUFFER_SIZE];
} Buffer;

Buffer* makeBuffer(const char* fname, const char* mode);
Buffer* makeBufferR(const char* fname);
Buffer* makeBufferW(const char* fname);
void    initBuffer(void *);
bool    open (void*, const char* fname, const char* mode);
bool    write(Buffer*);
bool    conditinalWrite(Buffer*);
bool    conditionalRead(Buffer*);
bool    read (Buffer*);
void    close(Buffer*);
void    fill(Buffer*, unsigned int);
unsigned int  next(Buffer*);
bool    good(void*);
void    increment(void *);

#endif /* BUFFER_H */
