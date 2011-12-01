#include "TSolEVIOFile.h"

#include "evioUtil.hxx"
#include "evioFileChannel.hxx"

#include "gemc_types.h"

#ifndef __CINT__


TSolEVIOFile::TSolEVIOFile(){
    fFilename[0] = '\0';
    fChan = NULL;

    printf("Hi, I'm an EVIO file!\n");

    return;
}

TSolEVIOFile::TSolEVIOFile(const char *f){
    printf("Hi, I'm an EVIO file for %s!\n", f);
    SetFilename(f);
    return;
}

TSolEVIOFile::~TSolEVIOFile(){
    printf("Goodbye, cruel world\n");
    return; 
}

void TSolEVIOFile::SetFilename( const char *f ){
   strcpy( fFilename, f );
   return;
}

Int_t TSolEVIOFile::Open(){
    // Return 0 on fail, 1 on success
    if( fFilename[0] == '\0' ){ return 0; }

    try {
	// SPR - Unclear to me what a good value to use for the 
	// buffer size is.  I picked something that works...
	fChan = new evio::evioFileChannel(fFilename, "r", 1<<24 );

	if( !fChan ){ 
	    fprintf(stderr, "%s: evioFileChannel could not be made\n",__PRETTY_FUNCTION__);
	    return 0; 
	}

	fChan->open();
    } catch (evio::evioException e) {
	// Problem opening
	fprintf(stderr, "%s %s line %d:  %s\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, e.toString().data() );
	return 0;
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
    // Return 1 on success

    // Channel not open
    if( !fChan ){ 
	fprintf(stderr, "%s %s line %d Channel not open\n",
		__FILE__,__PRETTY_FUNCTION__,__LINE__ );
	return 0; 
    }

    Clear();

    bool res;

    try {
	res = fChan->read(); 
    } catch (evio::evioException e) {
	// Problem opening
	fprintf(stderr, "%s:  %s\n",__PRETTY_FUNCTION__, e.toString().data() );
    }


    if( !res ){
	// Don't need to print this out.  Not really an error
#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Channel read return is false...  probably end of file\n",
		__FILE__, __FUNCTION__, __LINE__ );
#endif//DEBUG
	return 0;
    }

    fEvNum++;

    evio::evioDOMNodeList::iterator iter;

    evio::evioDOMTree EDT(fChan);

    // Read in header info
    // Not much interesting here
    evio::evioDOMNodeListP eHdrNodeList = EDT.getNodeList(evio::tagNumEquals(1, 1));

    //  Loop over these to get number of events
    for( iter = eHdrNodeList->begin(); iter != eHdrNodeList->end(); iter++ ){
	const evio::evioDOMNodeP ev = *iter;
	const vector<int> *vec = ev->getVector<int>();
	// Event number should be (*vec)[0];
    }


    
    // Extract Detector IDs
    evio::evioDOMNodeListP eDigNodeList = EDT.getNodeList(evio::tagNumEquals(__GEM_TAG, 100));
    for( iter = eDigNodeList->begin(); iter != eDigNodeList->end(); iter++ ){
#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Processing digitized (integer data) hits\n",
		__FILE__, __FUNCTION__, __LINE__ );
#endif//DEBUG

	if( (*iter)->isContainer() ){
	    ExtractDetIDs( (*iter)->getChildList(), __GEM_TAG  );
	}
    }


    // Raw events
    evio::evioDOMNodeListP eRawNodeList = EDT.getNodeList(evio::tagNumEquals(__GEM_TAG, 200));

    unsigned int i = 0;

    for( iter = eRawNodeList->begin(); iter != eRawNodeList->end(); iter++ ){
#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Processing raw (float data) hits\n",
		__FILE__, __FUNCTION__, __LINE__ );
#endif//DEBUG

	/////////////////////////////////////////////////////////////////////
	// Add events to data arrays
	if( (*iter)->isContainer() ){
	    BuildData( (*iter)->getChildList() );
	}
    }

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Completed\n",
		__FILE__, __FUNCTION__, __LINE__ );
#endif//DEBUG

    return 1;
}

void TSolEVIOFile::ExtractDetIDs( evio::evioDOMNodeList *hits, int tag  ){
    evio::evioDOMNodeList::const_iterator iter;
    int vnum;
    unsigned int i;

#ifdef  DEBUG
    fprintf(stderr, "%s %s line %d: Extracting IDs from 0x%08x\n",
	    __FILE__, __FUNCTION__, __LINE__, (unsigned int) hits);
#endif//DEBUG

    // Loop through everything with this tag
    for( iter = hits->begin(); iter != hits->end(); iter++ ){
	// Extract hit node
	const evio::evioDOMNodeP v = *iter;
	vnum = v->num;

	// Make sure this child leaf is actually data
	if( !v->isLeaf() ){
	    continue;
	}

	// num is our identifier for the data stored
	// Make sure this agrees with what we say it should be
	if( vnum != data_detid(tag) ){
	    continue;
	}

	vector<int> *vec = v->getVector<int>();

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Number of hits found %d\n",
		__FILE__, __FUNCTION__, __LINE__, vec->size());
#endif//DEBUG

	for( i = 0; i < vec->size(); i++ ){
	    fHitData.push_back( new hitdata((*vec)[i], data_size(tag)) );
	}
    }

    return;
}

void TSolEVIOFile::BuildData( evio::evioDOMNodeList *hits  ){
    evio::evioDOMNodeList::const_iterator iter;
    int vnum;
    unsigned int i;

    // Loop through everything with this tag
    for( iter = hits->begin(); iter != hits->end(); iter++ ){
	// Extract hit node

	const evio::evioDOMNodeP v = *iter;

	// Make sure this child leaf is actually data
	if( !v->isLeaf() ){
	    continue;
	}

	// vnum is the variable number
	vnum = v->num;
	vector<double> *vec = v->getVector<double>();

	// Build hit event using extracted detector IDs 
	for(  i = 0; i < vec->size(); i++ ){
	   fHitData[i]->SetData(vnum, (*vec)[i]); 
	}
    }

    return;
}

void TSolEVIOFile::Clear(){
    // Clear out hit data

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Deleting hits\n",
		__FILE__, __FUNCTION__, __LINE__);
#endif//DEBUG

    unsigned int i;
    for( i = 0; i < fHitData.size(); i++ ){
	delete fHitData[i];
    }

    fHitData.clear();

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Hits deleted\n",
		__FILE__, __FUNCTION__, __LINE__);
#endif//DEBUG

    return;
}

TSolGEMData *TSolEVIOFile::GetGEMData(){
    // Pack data into TSolGEMData
   
    unsigned int i;

    hitdata *h;

    TSolGEMData *gd = new TSolGEMData(GetNData());
    gd->SetEvent(fEvNum);
    gd->SetRun(0);

    for( i = 0; i < GetNData(); i++ ){
	h = GetHitData(i);

	// Vector information
	TVector3 p(h->GetData(20), h->GetData(21), h->GetData(22));
	gd->SetMomentum(i, p);

	TVector3 li(h->GetData(5), h->GetData(6), h->GetData(7));
	gd->SetHitEntrance(i, li);

	TVector3 lo(h->GetData(9), h->GetData(10), h->GetData(11));
	gd->SetHitExit(i, lo);

	TVector3 lr(h->GetData(2), h->GetData(3), h->GetData(4));
	gd->SetHitReadout(i, lr);

	///////////////////////////////

	gd->SetHitEnergy( i, h->GetData(1)*1e6 ); // Gives eV
	gd->SetHitChamber( i, h->GetDetID() );
	gd->SetParticleID( i, (UInt_t) h->GetData(18) );
	gd->SetParticleType( i, (UInt_t) h->GetData(13) );
    }

    return gd;
}


///////////////////////////////////////////////////////////////
// hitdata classes

hitdata::hitdata(int detid, int size ){
    fDetID = detid;
    fData  = new double[size];
    fSize  = size;
    fFillbits = 0;

    if( size > sizeof( long long int )*8 ){
	fprintf(stderr, "%s %s line %d:  Error:  Event size too long for bit pattern storage (requested %d, have %d)\n",
	      __FILE__, __PRETTY_FUNCTION__, __LINE__, size, 
	      sizeof(long long int)*8);
	exit(1);
    }

    // There is no value indexed at 0, so we'll just set it to 0 for
    // sanity's sake and not deal with crazy offsets all over

    fFillbits |= 1;
    fData[0] = 3.1415927;
}

void hitdata::SetData(int idx, double data ){
    if( idx < 0 || idx >= fSize ){
	fprintf(stderr, "%s %s line %d:  Error:  index out of range (%d oor of size %d)\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fSize);
	return;

    }

    fFillbits |= (1<<idx);

    fData[idx] = data;
    return;
}

double hitdata::GetData(int idx){
    if( idx < 0 || idx >= fSize ){
	fprintf(stderr, "%s %s line %d:  Error:  index out of range (%d oor of size %d)\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fSize);
	return 1e9;
    }

    if( !(fFillbits & (1<<idx)) ){
	fprintf(stderr, "%s %s line %d:  Error:  Accessing unset data (idx %d of 0x%08x) val: %f\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, (int) this, fData[idx] );
	return 1e9;
    }

    return fData[idx];
}

bool hitdata::IsFilled(){
    if( fFillbits == ((1<<fSize) - 1) ){
	return true;
    }

    return false;
}

hitdata::~hitdata(){
    delete fData;
}


#endif//__CINT__









