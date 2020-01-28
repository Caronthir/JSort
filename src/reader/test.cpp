#include "reader.h"
#include <iostream>

int main(int argc, char *argv[])
{
    std::string zinc = "/home/erdos/master/sortering/datafiles/sirius-20180125-101223.data";
    std::string si = "/home/erdos/master/sortering/datafiles/sirius-20180125-091200.data";
    Reader reader(si, "/home/erdos/gits/JSort/src/reader/data/output.bin");
    reader.read();
    reader.summary();
    return 0;
}
