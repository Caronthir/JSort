#ifndef EVENT_H
#define EVENT_H

#include <cstdint>
#include <iostream>


inline unsigned int makeMask(int start, int stop) {
    unsigned int m = 0;
    for (int i = start; i < stop; i++) {
        m = m | 1 << i;
    }
    return m;
}
inline unsigned int decodeMask(unsigned int j, int start, int stop){
    unsigned int mask = makeMask(start, stop);
    j = (j & mask) >> start;
    return j;
}

inline unsigned int encode(unsigned int front, unsigned int back, unsigned int labrsize){
    return front << 24 | back << 16 | labrsize;
}

inline void decode(unsigned int word, unsigned int* __restrict__ front, unsigned int* __restrict__ back, unsigned int* __restrict__ labrsize){
    *front    = decodeMask(word, 24, 32);
    *back     = decodeMask(word, 16, 24);
    *labrsize = decodeMask(word, 0, 16);
}

struct ADCTDC
{
    unsigned int channel;
    int adc;
    int tdc;
};


struct ADC
{
    unsigned int channel;
    unsigned int associated_channel;
    int data;
};

struct LaBrEvent {
public:
    LaBrEvent(){
        for (std::size_t i = 0; i < 32; ++i) {
            labr[i] = ADCTDC();
        }
    };
    virtual ~LaBrEvent(){};

    int e = 0;
    int de = 0;
    unsigned int front = 0;
    unsigned int back = 0;

    ADCTDC labr[32];
    unsigned int labr_num = 0;
    void clear(){
        labr_num=0;
    }
    inline void serialize(unsigned int * __restrict__ header, ADCTDC** body, unsigned int* __restrict__ size) {
        header[0] = e;
        header[1] = de;
        header[2] = encode(front, back, labr_num);
        *size = labr_num;
        // header[2] = front;
        // header[3] = back;
        // header[4] = labr_num;
        *body = labr;
    }
};

class Event
{
public:
    Event();
    virtual ~Event();

    /* I don't know how, but sometimes there are more than
       32 'hits' in a single event. */
    ADC e[64];
    ADC de[64];
    ADCTDC labr[64];

    // A map to correspond TDC and ADC to the same channel
    int labr_indices[64];
    std::size_t labr_num = 0;

    std::size_t e_num = 0;
    std::size_t de_num = 0;

    mutable bool rejected_ede_correlate  = false;
    mutable bool rejected_labr_correlate = false;
    bool rejected_ede_empty              = false;
    bool rejected_decoding               = false;
    bool rejected_pileup                 = false;
    bool pileup_flag                     = false;

    bool unpack(unsigned int* packet, std::size_t size);
    bool decodeWord(unsigned int word);
    void clear();
    void removePileup();
    bool correlateReject(LaBrEvent&) const;
private:
    inline void setADC(unsigned int channel, unsigned int data);
    inline void setTDC(unsigned int channel, unsigned int data);
    inline void addE  (unsigned int channel, unsigned int data);
    inline void addDE (unsigned int channel, unsigned int data);
    bool correlate(LaBrEvent &) const;
};


static const unsigned int EOB = 0x80000000;
// Extract the header
inline unsigned int boe(unsigned int x) { return ((x & 0xC0000000) >> 28); }
// The length of the packet
inline unsigned int ndw(unsigned int x) { return (x & 0x000000ff); }
// Which box fired
inline unsigned int box(unsigned int x) { return ((x & 0x3f800000) >> 23); }
// The channel
inline unsigned int chn(unsigned int x) { return ((x & 0x007f0000) >> 16); }
// The actual data
inline unsigned int dta(unsigned int x) { return (x & 0x0000ffff); }

inline bool isHeader(unsigned int x) { return boe(x) == 0xC; }
inline bool isEndofBuffer(unsigned int x) { return x == EOB; }
inline bool isGuardRing(unsigned int chn) {
  return (chn & 1) == 0 || chn >= 16;
}

inline bool isGood(const ADCTDC& x){return x.adc > 0 && x.tdc > 0;}
#endif /* EVENT_H */
