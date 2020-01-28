#include "buffer.h"


Buffer* makeBuffer(const char* fname, const char* mode){
    Buffer* buffer = malloc(sizeof(Buffer));
    initBuffer(buffer);
    open(buffer, fname, mode);
    return buffer;
}

Buffer* makeBufferR(const char* fname){
    Buffer* buffer = makeBuffer(fname, "rb");
    read(buffer);
    return buffer;
}

void initBuffer(void *buffer_) {
    Buffer* buffer = (Buffer*) buffer_;
    buffer->i = 0;
    buffer->size = BUFFER_SIZE;
    buffer->good = false;
}

bool open(void *buffer_, const char *fname, const char *mode) {
    Buffer* buffer = buffer_;
    buffer->handle = fopen(fname, mode);
    if (buffer->handle == NULL){
        printf("Could not read file %s", fname);
        buffer->good = false;
        return false;
    }
    buffer->good = true;
    return true;
}

void close(Buffer* buffer){
    free(buffer);
}

bool read(Buffer* buffer){
    buffer->size = fread(buffer->vals, sizeof(unsigned int), BUFFER_SIZE, buffer->handle);
    buffer->i = 0;
    return buffer->size == BUFFER_SIZE;
}

bool write(Buffer* buffer){
    if(buffer->i > buffer->size){
        printf("Buffer pointer is outside of buffer size");
        return false;
    }
    fwrite(buffer->vals, sizeof(unsigned int), buffer->i, buffer->handle);
    buffer->i = 0;
    return true;
}

bool conditionalRead(Buffer* buffer){
    if (buffer->i >= buffer->size) {
        return read(buffer);
    } else {
        return false;
    }
}

bool conditionalWrite(Buffer* buffer){
    if(buffer->i >= buffer->size){
        return write(buffer);
    } else {
        return false;
    }
}

unsigned int next(Buffer* buffer){
    if (!buffer->good){
        return 0;
    }
    unsigned int val = buffer->vals[buffer->i];
    increment(buffer);
    conditionalRead(buffer);
    return val;
}

inline void increment(void* buffer_){
    Buffer* buffer = (Buffer*)buffer_;
    buffer->i++;
    if (buffer->i >= buffer->size && buffer->size != BUFFER_SIZE){
        buffer->good = false;
    }
}

bool good(void *buffer) { return ((Buffer *)buffer)->good; };

