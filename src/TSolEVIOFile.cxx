#include "TSolEVIOFile.h"

#include "evioUtil.hxx"
#include "evioFileChannel.hxx"

#include "gemc_types.h"

TSolEVIOFile::TSolEVIOFile(){
    fFilename[0] = '\0';
    fChan = NULL;

    return;
}

TSolEVIOFile::TSolEVIOFile(const char *f){
    SetFilename(f);
    return;
}

TSolEVIOFile::~TSolEVIOFile(){ return; }

void TSolEVIOFile::SetFilename( const char *f ){
   strcpy( fFilename, f );
   return;
}

Int_t TSolEVIOFile::Open(){
    // Return 0 on fail, 1 on success
    if( fFilename[0] == '\0' ){ return 0; }

    try {
	fChan = new evio::evioFileChannel(fFilename, "r", 0 );

	if( !fChan ){ return 0; }

	fChan->open();
    } catch (evio::evioException e) {
	// Problem opening
	fprintf(stderr, "%s\n", e.toString().data() );
    }

    fEvNum = -1;

    return 1;
}

Int_t TSolEVIOFile::Close(){
    // Return 0 on fail, 1 on success
    try {
	if( !fChan ){ return 0; }

	fChan->close();
    } catch (evio::evioException e) {
	// Problem opening
	fprintf(stderr, "%s\n", e.toString().data() );
    }

    return 1;
}

Int_t TSolEVIOFile::ReadNextEvent(){
    // Channel not open
    if( !fChan ){ return 0; }

    Clear();

    if( !fChan->read() ){
	return 0;
    }

    fEvNum++;

    evio::evioDOMNodeList::iterator iter;

    evio::evioDOMTree EDT(fChan);

    // Read in header info
    // Not much interesting here
    evio::evioDOMNodeListP eHdrNodeList = EDT.getNodeList(evio::tagNumEquals(1, 1));
    
    // Extract Detector IDs
    evio::evioDOMNodeListP eDigNodeList = EDT.getNodeList(evio::tagNumEquals(__GEM_TAG, 100));
    for( iter = eDigNodeList->begin(); iter != eDigNodeList->end(); iter++ ){
	if( (*iter)->isContainer() ){
	    ExtractDetIDs( (*iter)->getChildList(), __GEM_TAG  );
	}
    }

    // Raw events
    evio::evioDOMNodeListP eRawNodeList = EDT.getNodeList(evio::tagNumEquals(__GEM_TAG, 200));

    unsigned int i = 0;

    for( iter = eRawNodeList->begin(); iter != eRawNodeList->end(); iter++ ){
	// Make sure that we have a number of events equal to the number of
	// detectors with hits as we found before when extracting the detector IDs
	if(i != fDetID.size() ){
	    fprintf(stderr, "%s number of detector IDs found does not equal number of data entries\n",
		    __PRETTY_FUNCTION__ );
	    return 0;
	}

	/////////////////////////////////////////////////////////////////////
	// Add events to data arrays
	if( (*iter)->isContainer() && i < fDetID.size() ){
	    BuildData( (*iter)->getChildList(), fDetID[i++], __GEM_TAG );
	}
    }

    return 0;
}

void TSolEVIOFile::ExtractDetIDs( evio::evioDOMNodeList *hits, int tag  ){
    evio::evioDOMNodeList::const_iterator iter;
    int vnum;

    // Loop through everything with this tag
    for( iter = hits->begin(); iter != hits->end(); iter++ ){
	// Extract hit node
	const evio::evioDOMNodeP v = *iter;
	vnum = v->num;

	// Make sure this child leaf is actually data
	if( !v->isLeaf() ) continue;
	// Make sure this is the size we expect
	if( vnum != data_size(tag) ) continue;

	vector<int> *vec = v->getVector<int>();

	// Extract detector ID and save

	fDetID.push_back( (*vec)[data_detid(tag)-1] );
    }

    return;
}

void TSolEVIOFile::BuildData( evio::evioDOMNodeList *hits, int slot, int tag  ){
    evio::evioDOMNodeList::const_iterator iter;
    int vnum, i;

    // Loop through everything with this tag
    for( iter = hits->begin(); iter != hits->end(); iter++ ){
	// Extract hit node

	const evio::evioDOMNodeP v = *iter;
	vnum = v->num;

	// Make sure this child leaf is actually data
	if( !v->isLeaf() ) continue;
	// Make sure this is the size we expect
	if( vnum != data_size(tag) ) continue;

	vector<double> *vec = v->getVector<double>();

	// Build CODA event using extracted detector IDs as slots
	// Fix crate number by tag, (multihit) channels are for each variable
	for(  i = 0; i < data_size(tag); i++ ){
	    AddDatum(tag, slot, i, (*vec)[i] );
	}
    }

    return;
}

void TSolEVIOFile::AddDatum( int crate, int slot, int chan, double datum ){
    // Look to see if this slot is already defined

    unsigned int i;
    for( i = 0; i < fHitData.size() && (fHitData[i]->GetSlot() != slot)
	   && (fHitData[i]->GetCrate() != crate); i++ ){}

    // Didn't find it... add to end
    if( i == fHitData.size() ){
	fHitData.push_back( new hitdata(crate, slot, data_size(crate) ) );
    }
}

void TSolEVIOFile::Clear(){
    fDetID.clear();

    // Clear out hit data

    unsigned int i;
    for( i = 0; i < fHitData.size(); i++ ){
	fHitData[i]->Clear();
    }

    return;
}

///////////////////////////////////////////////////////////////

hitdata::hitdata(int crate, int slot, unsigned int size){
    fCrate = crate;
    fSlot  = slot;
    fData  = new vector<double>[size];
}

void hitdata::AddDatum(int chan, double data ){
    fData[chan].push_back(data);
    return;
}

bool hitdata::IsFilled(){
    /*  Make sure that all the data arrays are the
     *  same size.  If they're not, we have a problem!
     */
    unsigned int first = fData[0].size();
    unsigned int i;
    for( i = 1; i < fSize && first == fData[i].size(); i++ ){}
    if( i == fSize ) return true;

    return false;
}

unsigned int hitdata::GetNHits(){
    if( !IsFilled() ){
	return 0;
    } else {
	return fData[0].size();
    }

    return 0;
}

void hitdata::Clear(){
    unsigned int i;
    for( i = 0; i < fSize; i++ ){
	fData[i].clear();
    }
    return;
}

hitdata::~hitdata(){
    delete fData;
}













