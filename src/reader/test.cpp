#include "reader.h"
#include <iostream>

int main(int argc, char *argv[])
{
    Reader reader("/home/erdos/master/sortering/datafiles/sirius-20180125-101223.data",
                  "/home/erdos/gits/JSort/src/reader/output.bin");
    reader.read();
    double total = (double)reader.total_words;
    std::cout << "Total words: " << reader.total_words << "\n"
              << "Total events: " << reader.num_events << "\n"
              << "Valid words: " << reader.valid_words/total*100 << "%\n"
              << "Unreadable words: " << reader.unreadable_words/total*100 << "%\n"
              << "Rejected words: " << reader.rejected_words/total*100 << "%" << std::endl;
    return 0;
}
