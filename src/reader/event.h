#ifndef EVENT_H
#define EVENT_H

#include <cstdint>
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

class LaBrEvent {
public:
    LaBrEvent(){};
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
};

class Event
{
public:
    Event();
    virtual ~Event();

    ADC e[32];
    ADC de[32];
    ADCTDC labr[32];

    // A map to correspond TDC and ADC to the same channel
    int labr_indices[32];
    std::size_t labr_num = 0;

    std::size_t e_num = 0;
    std::size_t de_num = 0;

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
