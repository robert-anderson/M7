//
// Created by Robert John Anderson on 2020-02-13.
//

#ifndef M7_SHAREDTABLEFILLER_H
#define M7_SHAREDTABLEFILLER_H

#include "Table.h"

/*
class StripingFiller {
public:
    StripingFiller(){

    }
};


enum SharingStrategy{
    private_temp,
    shared_stripe
};

template<enum SharingStrategy strategy>
class SharedTableFiller {
    using filler_t = typename std::conditional<strategy==private_temp, Table, StripingFiller>::type;
    filler_t filler;

public:
    SharedTableFiller(Table &table):m_table(table){};

};
*/


#endif //M7_SHAREDTABLEFILLER_H
