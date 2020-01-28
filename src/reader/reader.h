#ifndef READER_H
#define READER_H

#include <istream>
#include <ostream>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "event.h"

//#define BUFFER_SIZE 32'768
const size_t BUFFER_SIZE = 32'768;
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
    std::string summary();

    std::string inputpath;

    unsigned int input[BUFFER_SIZE];
    LaBrEvent   output[BUFFER_SIZE_OUT];
    std::size_t current_output = 0;

    FILE* source;
    FILE* destination;

    Event event = Event();

    unsigned int num_events              = 0;
    unsigned int valid_words             = 0;
    unsigned int unreadable_words        = 0;
    unsigned int rejected_words          = 0;
    unsigned int total_words             = 0;

    unsigned int rejected_ede_correlate  = 0;
    unsigned int rejected_labr_correlate = 0;
    unsigned int rejected_ede_empty      = 0;
    unsigned int rejected_decoding       = 0;
    unsigned int rejected_pileup         = 0;

    size_t filesize = 0;
private:
    bool readPacket(size_t& i);
    void updateCounters();
};



#endif /* READER_H */
