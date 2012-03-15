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
#define __GENERATED_TAG  10
#define __GEM_TAG  110
#define __FLUX_TAG 800

#define __GENERATED_SIZE 7

#define __GEM_DRIFT_ID 6
#define __GEM_STRIP_ID 19

int __gemc_types_datasize[NBANKS] = {23,21};
int __gemc_types_detidnum[NBANKS] = {23,21};

// FIXME:  Need to do this better,  map?
int data_size(int tag){
    if( tag == __GEM_TAG ) {
	return __gemc_types_datasize[0];
    }

    if( tag == __FLUX_TAG ) {
	return __gemc_types_datasize[1];
    }

    return 0;
}

int data_detid(int tag){
    if( tag == __GEM_TAG ) {
	return __gemc_types_detidnum[0];
    } 

    if( tag == __FLUX_TAG ) {
	return __gemc_types_detidnum[1];
    }

    return -1;
}

#endif//__GEMC_TYPES_H
