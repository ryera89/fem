#ifndef ELEMENT_H
#define ELEMENT_H

#include "node.h"
#include <array>

enum class ELEMENT_TYPE{TRI,QUAD};

template<ELEMENT_TYPE etype,uint8_t nnum>
struct element{
    static constexpr ELEMENT_TYPE type = etype;
    static constexpr uint8_t nod_number = nnum;
    std::array<node,nnum> nodes;
};


#endif // ELEMENT_H
