#ifndef READER_H
#define READER_H

#include <istream>
#include <ostream>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "event.h"

#define BUFFER_SIZE 32'768
const size_t BUFFER_SIZE_OUT = BUFFER_SIZE/sizeof(LaBrEvent);

using cstring = const std::string&;

class Reader
{
public:
    Reader(cstring readpath, cstring savepath);
    virtual ~Reader();

    bool read();
    bool readChunk();
    void write();

    std::string inputpath;

    unsigned int input[BUFFER_SIZE];
    LaBrEvent   output[BUFFER_SIZE_OUT];
    // The entire input buffer can not be used,
    // as a packet could be split at the end.
    std::size_t last_packet;
    std::size_t current_output;

    FILE* source;
    FILE* destination;

    Event event = Event();

    unsigned int num_events = 0;
    unsigned int valid_words = 0;
    unsigned int unreadable_words = 0;
    unsigned int rejected_words = 0;
    unsigned int total_words = 0;

    size_t filesize = 0;
private:
    bool readPacket(size_t& i);
};



#endif /* READER_H */
