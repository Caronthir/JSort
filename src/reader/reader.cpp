#include "reader.h"
#include "event.h"
#include <filesystem>
#include <stdexcept>
#include <stdio.h>
#include <iostream>

//#define DOLOG2
#ifdef DOLOG2
#define LOG(x) std::cout << x << std::endl;
#define DBGS(x) do{x}while(0);
#define HALT() do{char tmp; std::cin >> tmp;}while(0);
#else
#define LOG(x)                                                                 \
  while (0) {                                                                  \
  };
#define DGBS(x)                                                                \
  while (0) {                                                                  \
  };
#define HALT() while(0){};
#endif

Reader::Reader(cstring readpath, cstring savepath){
    inputpath = std::filesystem::path(readpath).filename();
    source = fopen(readpath.c_str(), "rb");
    filesize = std::filesystem::file_size(readpath);
    //destination.open(savepath, std::ios::out | std::ios::binary | std::ios::trunc);
    destination = fopen(savepath.c_str(), "wb");

    for (std::size_t i = 0; i < BUFFER_SIZE_OUT; ++i) {
        output[i] = LaBrEvent();
    }

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
    printf("Going to fill buffer %.2f times\n", num_buffers);
    printf("\e[?25l    ");
    while(readChunk()){
        // Each chunk is BUFFER_SIZE long
        total_words += BUFFER_SIZE;
        size_t i = 0;
        //if(100*counter/num_buffers >= 94){break;}
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
            updateCounters();
        }
        counter++;
        //if (counter == 5064)
        //    break;
        //printf("\r\r\r\r\r\r\r\r\r\r%10i", total_words);
        printf("\r\r\r\r%3.1f%%", 100*counter/num_buffers);
    }
    printf("\e[?25h\n");
    return true;
}

bool Reader::readPacket(size_t& i){
    // Read the packet and move the index to the next word after the packet end.
    const unsigned int packet_size = ndw(input[i]);
    if (i + packet_size >= BUFFER_SIZE) {
      std::cerr << "Incomplete chunk. Aborting" << std::endl;
      return false;
    }
    i++;
    bool success = event.unpack(&input[i], packet_size);
    if (success) {
        // Correlate event and put into output buffer on success.
        success = event.correlateReject(output[current_output]);
        //success = false;
        if (success) {
            num_events++;
            valid_words += packet_size+1;  // +1 for the packet size
            current_output++;
            if (current_output >= BUFFER_SIZE_OUT)
                write();
        }
    }

    if (!success){
      rejected_words += packet_size + 1;
    }
    i += packet_size;
    return true;
}

bool Reader::readChunk(){
    size_t words_read = fread(input, sizeof(unsigned int), BUFFER_SIZE, source);
    LOG("=*= READ " << words_read << " WORDS =*=");
    return words_read == BUFFER_SIZE;
}

void Reader::write(){
    if(current_output > BUFFER_SIZE_OUT)
        throw(std::out_of_range("Buffer pointer outside of buffer"));

    // e, de, front, back, labr_num
    unsigned int header[3];
    unsigned int size;
    ADCTDC* body;
    for (std::size_t i = 0 ; i < current_output; ++i) {
        output[i].serialize(header, &body, &size);
        LOG(header[0] << "\t" << header[1] << "\t"
            << header[2] << "\t" << header[3] << "\t"
            << header[4]);
        HALT();

        // DBGS(for(int i = 0; i < header[4]; i++){};
        //     std::cout << body[i].channel << "\t"
        //               << body[i].adc << "\t" << body[i].tdc << std::endl;
        //     });

        fwrite(header, sizeof(header), 1, destination);
        fwrite(body, sizeof(ADCTDC), size, destination);
    }
    current_output = 0;
}

void Reader::updateCounters(){
    rejected_decoding       += event.rejected_decoding;
    rejected_ede_correlate  += event.rejected_ede_correlate;
    rejected_ede_empty      += event.rejected_ede_empty;
    rejected_labr_correlate += event.rejected_labr_correlate;
}

std::string Reader::summary(){
  double total = (double)total_words;
  double R = (double)rejected_words;
  R = (double)(rejected_decoding+rejected_ede_correlate+rejected_ede_empty+rejected_labr_correlate);
  std::cout.imbue(std::locale(""));
  std::cout << "Total words:      " << total_words << "\n"
            << "Total events:     " << num_events << "\n"
            << "Valid words:      " << valid_words / total * 100 << "%\n"
            << "Unreadable words: " << unreadable_words / total * 100
            << "%\n"
            << "Rejected words:   " << rejected_words / total * 100 << "%\n"
            << "    where\n"
            << "    EΔE Correlate:  " << rejected_ede_correlate/R*100 << "%\n"
            << "    LaBr Correlate: " << rejected_labr_correlate/R*100 << "%\n"
            << "    Decoding:       " << rejected_decoding/R*100 << "%\n"
            << "    EΔE Empty:      " << rejected_ede_empty/R*100 << "%\t( " << rejected_ede_empty <<  " )\n"
            << "    Pileup:         " << rejected_pileup/(float)rejected_ede_empty
            << "% of EΔE Empty  ( " << rejected_pileup  << " ) [ " << rejected_ede_empty - rejected_pileup << " ]\n"
            << std::endl;
  return "";
}
