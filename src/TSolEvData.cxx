#include "TSolEvData.h"
#include "TSolEVIOFile.h"

TSolEvData::TSolEvData(){
}

TSolEvData::~TSolEvData(){
}

void TSolEvData::Clear(){
    THaEvData::Clear();
}

int TSolEvData::LoadEvent( const int *, THaCrateMap * ){
    // Just not do anything here.  We're never going to 
    // call it in the general analyzer
    return 0;
}

int TSolEvData::LoadEvent( TSolEVIOFile *f, THaCrateMap *map ){
    // This is where we load up the data structure

    unsigned int i, chan, hit;
    Int_t crate, slot;

    // Initialize cratemaps
    if (first_decode) {
	fMap = map;
	init_cmap();
	if (init_slotdata(fMap) == HED_ERR) return HED_ERR;
	first_decode = false;
    }

    // Clear out previous data
    Clear();
    for( i = 0; i < (unsigned int) fNSlotClear; i++ ){
	crateslot[fSlotClear[i]]->clearEvent();
    }

    // Set basic information about our event
    // We just have normal "physics" events
    evscaler = 0;
    event_length = 0;
    event_type = 1;
    event_num = f->GetEvNum();

    // Copy data from our file into the THaEvData
    // structure

    f->ReadNextEvent();

    // We have several hitdata objects
    // they are all guaranteed to have unique
    // crate/slot numbers

    hitdata *data;
    std::vector <double> eviodata;
    Int_t raw;
    Float_t fraw;

    /*
    for( i = 0; i < f->GetNData(); i++ ){
	data  = f->GetHitData(i);
	crate = data->GetCrate();
	slot  = data->GetSlot();

	// Empty or malformed data just get ignored
	if( data->GetNHits() < 1 ) continue;

	for( chan = 0; chan < data->GetSize(); chan++ ){
	    eviodata = data->GetDataArray(chan);

	    // We want be SURE to do this without a cast
	    // but it needs to become a float first so we can pack it
	    // into an integer for transport
	    // This is a dumb kluge and I feel bad and I'm sorry - SPR
	    for( hit = 0; hit < eviodata.size(); hit++ ){
		fraw = (Float_t) eviodata[hit];
		memcpy( &raw, &fraw, sizeof(Float_t) );
		if (crateslot[idx(crate,slot)]->loadData("evio",chan,raw,raw)
			== SD_ERR) return HED_ERR;
	    }
	}
    }
    */

    // Et, viola!

    return HED_OK;
}
