#include "event.h"
#include <iostream>

//#define DOLOG 1
#ifdef DOLOG
#define LOG(x) std::cout << x << std::endl;
#define DBGS(x) x;
#define HALT() char tmp; std::cin >> tmp;
#else
#define LOG(x) while(0){};
#define DGBS(x) while(0){};
#define HALT() while(0){};
#endif


Event::Event() {
    for(std::size_t i = 0; i < 32; i++){
        e[i] = ADC();
        de[i] = ADC();
        labr[i] = ADCTDC();
        labr_indices[i] = -1;
    }
}

Event::~Event() {}

inline void Event::setADC(unsigned int channel, unsigned int data){
    int index = labr_indices[channel];
    if (index < 0) {
        index = labr_indices[channel] = labr_num++;
        labr[index].channel = channel;
        labr[index].tdc = 0;
    }
    labr[index].adc = data;
}

inline void Event::setTDC(unsigned int channel, unsigned int data){
    int index = labr_indices[channel];
    if (index < 0){
        index = labr_indices[channel] = labr_num++;
        labr[index].channel = channel;
        labr[index].adc = 0;
    }
    labr[index].tdc = data;
}

inline void Event::addE(unsigned int channel, unsigned int data){
    ADC& e_ = e[e_num++];
    e_.channel = channel >> 1;  // Back
    e_.data = data;
}

inline void Event::addDE(unsigned int channel, unsigned int data) {
  ADC& de_ = de[de_num++];
  de_.channel = channel % 8;  // Front
  // Note the integer division
  de_.associated_channel = channel / 8; // Must == e back
  de_.data = data;
}

void Event::clear(){
    e_num = 0;
    de_num = 0;
    for(std::size_t i = 0; i < labr_num; i++){
        auto channel = labr_indices[i];
        labr[channel].adc = 0;
        labr[channel].tdc = 0;
    }
    for(std::size_t i = 0; i < 32; i++){
        labr_indices[i] = -1;
    }
    labr_num = 0;
    rejected_labr_correlate = false;
    rejected_ede_correlate  = false;
    rejected_ede_empty      = false;
    rejected_decoding       = false;
    rejected_pileup         = false;
    pileup_flag             = false;
}


bool Event::unpack(unsigned int *packet, std::size_t size){
    clear();
    HALT();
    LOG("======= CLEAR ====== Size: " << size);
    unsigned int word; //nextword
    for (std::size_t i = 0; i < size; ++i) {
        word = packet[i];
        // This can fail V
        //nextword = packet[i + 1];
        // Check for corrupt header
        if (boe(word) != 0) {
            LOG(">> Misplaced header <<");
            continue;
        }

        if (!decodeWord(word)){
            rejected_decoding = true;
            return false;
        }
    }

    // Reject the event if there are more than 1 Δe event in the
    // same detect, i.e. pileup
    if (e_num > 0 && de_num > 1){
        pileup_flag = true;
        removePileup();
    }

    // Reject if the front counter or back counter is empty
    // _after_ pileup is removed
    if (e_num == 0 || de_num == 0){
        rejected_ede_empty = true;
        if(pileup_flag){
            rejected_pileup = true;
            LOG("++++ PILEUP REJECT ++++")
                } else {
            LOG("++++ EDE REJECT ++++")
        }
        return false;
    }

    return true;
}

bool Event::decodeWord(unsigned int word){
    const auto box_num = box(word);
    const auto channel = chn(word);
    const auto data = dta(word);

    switch (box_num) {
    case 0x00: {
        break;
    }
    case 0x01: {
        break;
    }
    case 0x02: {
        break;
    }
    case 0x10: {
        setTDC(channel, data);
        break;
    }
    case 0x20: case 0x24: {
        // Energy of LaBr ch 0-31, MADC ch 0-31
        // Reject if adc < 0
        //if (data < 0){
        //    return true;
        //}
        setADC(channel, data);
        break;
    }
    case 0x21: {
        // Energy of E ch 0-31
        if (isGuardRing(channel))
            return true;
        LOG(" e  "  << data)
        addE(channel, data);
        break;
    }
    case 0x22: {
        // Energy of Δe ch 0-31
        addDE(channel, data);
        LOG("Δe1 " << data)
        break;
    }
    case 0x23: {
        // Energy of Δe2 ch 32-61
        addDE(channel+32, data);
        LOG("Δe2 " << data)
        break;
    }
    default: {
        return false;
    }
    }
    return true;
}

void Event::removePileup(){
    // Remove multiple Δe hitting the save back detector
    bool bad_indices[32] = {false};
    unsigned int bad_elements = 0;
    for (std::size_t i = 0; i < de_num; ++i) {
        for (std::size_t j = i+1; j < de_num; ++j) {
            if (de[i].associated_channel == de[j].associated_channel){
                bad_indices[i] = bad_indices[j] = true;
                bad_elements += 2;
            }
        }
    }

    if (bad_elements == 0)
        return;
    if (bad_elements == de_num){
        // No point in trying to fix the problem, just reject it later.
        de_num = 0;
        return;
    }

    auto length = de_num;
    de_num = 0;

    // Move all the good elements down. A side effect is that the
    // elements not overwritten remains, but this should not be a
    // problem as the de_num points to the last viable element
    for (std::size_t i = 0; i < length; ++i) {
        if(!bad_indices[i]){
            de_num++;
            if(i != de_num){
                de[de_num].channel = de[i].channel;
                de[de_num].data = de[i].data;
                de[de_num].associated_channel = de[i].associated_channel;
            }
        }
    }
}

bool Event::correlate(LaBrEvent& event) const{
    const ADC* de_c;
    const ADC* e_c;
    unsigned int back, front;
    bool found = false;
    LOG("#e " << e_num << "\t#de " << de_num)
    for (std::size_t i = 0; i < de_num; ++i) {
        de_c = &de[i];
        for (std::size_t j = 0; j < e_num; ++j) {
            e_c = &e[j];
            back = e_c->channel;
            LOG("e back " << back << "\t de asoch " << de_c->associated_channel)
            if (de_c->associated_channel != back)
                continue;
            front = de_c->channel;
            LOG("de front: " << front);
            // Reject if there are more than 1 possible correlation
            if (found)
                return false;
            found = true;
            event.front = front+1;
            event.back = back+1;
            event.de = de_c->data;
            event.e = e_c->data;
        }
    }
    return found;
}

bool Event::correlateReject(LaBrEvent& labrevent) const{
    labrevent.clear();
    bool accept = correlate(labrevent);
    if (!accept) {
        LOG("++++ EDE CORR ++++");
      rejected_ede_correlate = true;
      return false;
    }
    accept = false;

    const ADCTDC *labr_c;
    LOG("#labr " << labr_num);
    for (std::size_t i = 0; i < labr_num; i++) {
        labr_c = &labr[i];
        LOG("  (" << i << ")\tchn  " << labr_c->channel
            << "\tadc " << labr_c->adc << "\ttdc  " << labr_c->tdc
            << "\tgood?  " << isGood(*labr_c));
        if (isGood(*labr_c)) {
            accept = true;
            labrevent.labr[labrevent.labr_num] = *labr_c;
            labrevent.labr_num++;
        }
    }
    rejected_labr_correlate = !accept;
    if (rejected_labr_correlate){
        LOG("++++ LABR CORR ++++");
    }

    return accept;
}
