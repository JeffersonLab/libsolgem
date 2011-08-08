#ifndef __GEMC_TYPES_H
#define __GEMC_TYPES_H

////////////////////////////////////////////////////////
//  Data for extracting things from GEMC
//
//  we'll hardcode them here, but it would be nice to
//  maybe get them into a database
//  I guess this could also be done through a mysql
//  interface, but I think that makes it more complicated
//  and breakable

#define NBANKS 2

// Tag numbers associated in the GEMC banks
#define __GEM_TAG  110
#define __FLUX_TAG 800


int __gemc_types_datasize[NBANKS] = {22,20};
int __gemc_types_detidnum[NBANKS] = {23,21};

int data_size(int tag){
    return __gemc_types_datasize[tag];
}

int data_detid(int tag){
    return __gemc_types_detidnum[tag];
}

#endif//__GEMC_TYPES_H
