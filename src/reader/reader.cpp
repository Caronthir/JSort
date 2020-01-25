#include "reader.h"
#include "event.h"
#include <filesystem>
#include <stdio.h>

Reader::Reader(cstring readpath, cstring savepath){
    inputpath = std::filesystem::path(readpath).filename();
    source = fopen(readpath.c_str(), "rb");
    filesize = std::filesystem::file_size(readpath);
    //destination.open(savepath, std::ios::out | std::ios::binary | std::ios::trunc);
    destination = fopen(savepath.c_str(), "wb");
}

Reader::~Reader(){
    fclose(source);
    fclose(destination);
}

bool Reader::read(){
    unsigned int word;
    double num_buffers =  filesize/((float)BUFFER_SIZE*sizeof(unsigned int));
    int counter = 0;
    printf("%s\n", inputpath.c_str());
    printf("\e[?25l    ");
    while(readChunk()){
        total_words += BUFFER_SIZE;
        size_t i = 0;
        if(100*counter/num_buffers >= 94){break;}
        while(i < BUFFER_SIZE){
            word = input[i];
            if(isHeader(word)){
                // Fails if the buffer contains an incomplete chunk,
                // in which case the rest will be unreadable
                if(!readPacket(i))
                    break;
            } else if(isEndofBuffer(word)) {
                i++;
                valid_words++;
            } else {
                // Something went wrong. Ignore the rest of the packet
                while(i < BUFFER_SIZE && !isHeader(input[i])){
                    i++;
                    unreadable_words++;
                }
            }
        }
        counter++;
        printf("\r\r\r\r%3.1f%%", 100*counter/num_buffers);
    }
    printf("\e[?25h\n");
    return true;
}

bool Reader::readPacket(size_t& i){
    // Read the packet and move the index to the next word after the packet end.
    unsigned int packet_size = ndw(input[i]);
    if (i + packet_size >= BUFFER_SIZE) {
      std::cerr << "Incomplete chunk. Aborting" << std::endl;
      return false;
    }
    bool success = event.unpack(&input[i + 1], packet_size);
    if (success) {
        // Correlate event and put into output buffer on success.
        //success = event.correlateReject(output[current_output]);
        success = false;
        if (success) {
            num_events++;
            valid_words += packet_size+1;
            current_output++;
            if (current_output >= BUFFER_SIZE_OUT)
                write();
        }
    } else {
      rejected_words++;
    }
    i += packet_size + 1;
    return true;
}

bool Reader::readChunk(){
    size_t words_read = fread(input, sizeof(unsigned int), BUFFER_SIZE, source);
    return words_read == BUFFER_SIZE;
}

void Reader::write(){
    fwrite(&output, sizeof(LaBrEvent), current_output, destination);
    current_output = 0;
}
