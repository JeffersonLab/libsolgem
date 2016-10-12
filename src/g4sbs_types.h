#ifndef __G4SBS_TYPES_H
#define __G4SBS_TYPES_H

////////////////////////////////////////////////////////
//  Data for extracting things from G4SBS
//
//  we'll hardcode them here, but it would be nice to
//  maybe get them into a database
//  I guess this could also be done through a mysql
//  interface, but I think that makes it more complicated
//  and breakable


#define NBANKS 2  // => ???

// Tag numbers associated in the G4SBS banks => ???
#define __GENERATED_TAG  10
#define __GEM_TAG  110
#define __FLUX_TAG 800

#define __GENERATED_SIZE 7

// hits only come from drift (so far) in g4sbs
//#define __GEM_DRIFT_ID 5
//#define __GEM_COPPER_FRONT_ID 3
//#define __GEM_COPPER_BACK_ID 6 //?
//#define __GEM_STRIP_ID 18

// warning: will have, at some point, to correspond to TSolGEMData....
int __g4sbs_types_datasize[NBANKS] = {23, 21};
int __g4sbs_types_detidnum[NBANKS] = {23, 21};


// FIXME:  Need to do this better,  map?
int g4sbs_data_size(int tag){
    if( tag == __GEM_TAG ) {
	return __g4sbs_types_datasize[0];
    }
    if( tag == __FLUX_TAG ) {
    	return __g4sbs_types_datasize[1];
    }
    return 0;
}

int g4sbs_data_detid(int tag){
    if( tag == __GEM_TAG ) {
	return __g4sbs_types_detidnum[0];
    } 
    if( tag == __FLUX_TAG ) {
    	return __g4sbs_types_detidnum[1];
    }
    return -1;
}

#endif//__G4SBS_TYPES_H
