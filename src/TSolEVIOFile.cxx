#include "TSolEVIOFile.h"

#include "evioUtil.hxx"
#include "evioFileChannel.hxx"

#include "gemc_types.h"

#ifndef __CINT__


TSolEVIOFile::TSolEVIOFile(){
    fFilename[0] = '\0';
    fChan = NULL;

    return;
}

TSolEVIOFile::TSolEVIOFile(const char *f){
    SetFilename(f);
    return;
}

TSolEVIOFile::~TSolEVIOFile(){
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

    bool res = false;

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
	//const evio::evioDOMNodeP ev = *iter;
	//const vector<int> *vec = ev->getVector<int>();
	// Event number should be (*vec)[0];
    }

    // Extract generated data
    evio::evioDOMNodeListP eGenNodeList = EDT.getNodeList(evio::tagNumEquals(__GENERATED_TAG, 200));

    for( iter = eGenNodeList->begin(); iter != eGenNodeList->end(); iter++ ){
#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Processing raw (float data) generated tracks\n",
		__FILE__, __FUNCTION__, __LINE__ );
#endif//DEBUG

	/////////////////////////////////////////////////////////////////////
	// Add events to data arrays
	if( (*iter)->isContainer() ){
	    BuildGenerated( (*iter)->getChildList() );
	}
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

void TSolEVIOFile::BuildGenerated( evio::evioDOMNodeList *hits  ){
    evio::evioDOMNodeList::const_iterator iter;
    int vnum;
    unsigned int i;

    // Loop through everything with this tag
    for( iter = hits->begin(); iter != hits->end(); iter++ ){
	// Extract node

	const evio::evioDOMNodeP v = *iter;

	// Make sure this child leaf is actually data
	if( !v->isLeaf() ){
	    continue;
	}

	// vnum is the variable number
	vector<double> *vec = v->getVector<double>();

	if( fGenData.size() < vec->size() ){
	    for(  i = 0; i < vec->size(); i++ ){
		fGenData.push_back( new gendata() );
	    }
	}

	vnum = v->num;
	// Build hit event 
	for(  i = 0; i < vec->size(); i++ ){
	    // 4,5,6 are vertex positions.  GEMC spits them out in cm
	    // which we will convert to mm for consistancy
	    if( 4 <= vnum && vnum <= 6 ){
		fGenData[i]->SetData(vnum/10-1, (*vec)[i]*10.0); 
	    } else {
		fGenData[i]->SetData(vnum/10-1, (*vec)[i]); 
	    }
	}
    }

    return;
}

void TSolEVIOFile::BuildData( evio::evioDOMNodeList *hits ){
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

	// Build hit event 
	for(  i = 0; i < vec->size(); i++ ){
	    fHitData[i]->SetData(vnum, (*vec)[i]); 
	}

	// Extract weight data which is stored in each hit 
	// and put it into generated data.  This is a dumb kludge
	// and it makes me feel bad.  SPR 4/5/2012
	if( vnum == 19 ){
	    for( i = 0; i < fGenData.size() && vec->size()>0; i++ ){
		fGenData[i]->SetData(7, (*vec)[0]);
	    }
	}
    }

    return;
}

void TSolEVIOFile::Clear(){
    // Clear out hit and generated data

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Deleting hits\n",
		__FILE__, __FUNCTION__, __LINE__);
#endif//DEBUG

    unsigned int i;
    for( i = 0; i < fHitData.size(); i++ ){
	delete fHitData[i];
    }

    for( i = 0; i < fGenData.size(); i++ ){
	delete fGenData[i];
    }

    fHitData.clear();
    fGenData.clear();

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Hits deleted\n",
		__FILE__, __FUNCTION__, __LINE__);
#endif//DEBUG

    return;
}

TSolGEMData *TSolEVIOFile::GetGEMData(){
    // Pack data into TSolGEMData
   
    unsigned int i,j;

    hitdata *h, *hs;

    bool matchedstrip;

    // Upper limit of what we need.  Probably only a 
    // factor of 2 high
    TSolGEMData *gd = new TSolGEMData(GetNData());
    gd->SetEvent(fEvNum);
    gd->SetRun(0);

//    printf("NEXT EVENT ---------------------------\n");

    if (GetNData() == 0){
	gd->SetNHit (0);
    }
    else
    {
	int ngdata = 0;
	for( i = 0; i < GetNData(); i++ ){
	    h = GetHitData(i);

	    // Chamber IDs are tagged as 
	    // xxxyy  where xxx is the chamber num and yy is the 
	    // plane num  we find the drift planes and then
	    // match them to the corresponding readout hits

	    if( h->GetDetID()%100 == __GEM_DRIFT_ID &&  h->GetData(1)>0.0 ){
		// Vector information
		TVector3 p(h->GetData(20), h->GetData(21), h->GetData(22));
		gd->SetMomentum(ngdata, p);

		TVector3 li(h->GetData(5), h->GetData(6), h->GetData(7));
		gd->SetHitEntrance(ngdata, li);

		TVector3 lo(h->GetData(9), h->GetData(10), h->GetData(11));
		gd->SetHitExit(ngdata, lo);

		//	    printf("%d %f %f\n", h->GetDetID()/100, li.X(), li.Y()  );

		gd->SetHitEnergy(ngdata, h->GetData(1)*1e6 ); // Gives eV
		gd->SetParticleID(ngdata, (UInt_t) h->GetData(18) );
		gd->SetParticleType(ngdata, (UInt_t) h->GetData(13) );

		// Chamber ID starts indexing a 0 whereas we start conventionally
		// at 1 
		gd->SetHitChamber(ngdata, h->GetDetID()/100 - 1 );

		////////////////////////////////////////////
		// Search for entrance/exit hits in surrounding Cu plane

		for( j = 0; j < GetNData(); j++ ){
		    hs = GetHitData(j);

		    if( hs->GetDetID()%100 == __GEM_COPPER_FRONT_ID &&    // is prior Cu plane
			    hs->GetDetID()/100 == h->GetDetID()/100 && // same detector
			    ((UInt_t) hs->GetData(18)) == ((UInt_t) h->GetData(18))  // same particle
		      ){
			// Found matching hit, replace entrance data
			li = TVector3(hs->GetData(5), hs->GetData(6), hs->GetData(7));
			gd->SetHitEntrance(ngdata, li);
			break;
		    }
		}

		for( j = 0; j < GetNData(); j++ ){
		    hs = GetHitData(j);

		    if( hs->GetDetID()%100 == __GEM_COPPER_BACK_ID &&    // is subsequent Cu plane
			    hs->GetDetID()/100 == h->GetDetID()/100 && // same detector
			    ((UInt_t) hs->GetData(18)) == ((UInt_t) h->GetData(18))  // same particle
		      ){
			// Found matching hit, replace exit data
			lo = TVector3(hs->GetData(9), hs->GetData(10), hs->GetData(11));
			gd->SetHitExit(ngdata, lo);
			break;
		    }
		}

		////////////////////////////////////////////
		// Search other hits for the corresponding 
		// hit on the strip


		matchedstrip = false;
		for( j = 0; j < GetNData(); j++ ){
		    hs = GetHitData(j);

		    if( hs->GetDetID()%100 == __GEM_STRIP_ID &&    // is strip plane
			    hs->GetDetID()/100 == h->GetDetID()/100 && // same detector
			    ((UInt_t) hs->GetData(18)) == ((UInt_t) h->GetData(18))  // same particle
		      ){
			if( !matchedstrip ){
			    // This is the truth information
			    TVector3 lr(hs->GetData(2), hs->GetData(3), hs->GetData(4));
			    gd->SetHitReadout(ngdata, lr);
			    matchedstrip = true;
			} else {
			    fprintf(stderr, "%s %s line %d: Found multiple readout plane hits matching drift hit.  Truth information may be inaccurate\n",
				    __FILE__, __FUNCTION__, __LINE__);
			}
		    } 
		}

		if( !matchedstrip && (UInt_t) h->GetData(18) == 1 ){
		    // FIXME
		    // SPR 12/2/2011
		    // This isn't the greatest way to do this but we're usually intersted in 
		    // Particle 1 when we're looking at doing tracking.  Not all things depositing
		    // energy leave a corresponding hit in the cathode plane.  Maybe we can look at
		    // the mother IDs or something later.  

		    fprintf(stderr, "%s %s line %d: Did not find readout plane hit corresponding to drift hits.  No truth information\n",
			    __FILE__, __FUNCTION__, __LINE__);
		    TVector3 lr(-1e9, -1e9, -1e9);
		    gd->SetHitReadout(ngdata, lr);
		}

		ngdata++;
	    }
	}
	gd->SetNHit(ngdata);
    }

    return gd;
}


///////////////////////////////////////////////////////////////
// hitdata classes

hitdata::hitdata(int detid, unsigned int size ){
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

void hitdata::SetData(unsigned int idx, double data ){
    if( idx >= fSize ){
	fprintf(stderr, "%s %s line %d:  Error:  index out of range (%d oor of size %d)\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fSize);
	return;

    }

    fFillbits |= (1<<idx);

    fData[idx] = data;
    return;
}

double hitdata::GetData(unsigned int idx) const {
    if( idx >= fSize ){
	fprintf(stderr, "%s %s line %d:  Error:  index out of range (%d oor of size %d)\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fSize);
	return 1e9;
    }

    if( !(fFillbits & (1<<idx)) ){
	fprintf(stderr, "%s %s line %d:  Error:  Accessing unset data (idx %d) val: %f\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fData[idx] );
	return 1e9;
    }

    return fData[idx];
}

bool hitdata::IsFilled() const {
    if( fFillbits == ((1<<fSize) - 1) ){
	return true;
    }

    return false;
}

hitdata::~hitdata(){
    delete fData;
}

///////////////////////////////////////////////////////////////
// gendata classes

// Size is 1 bigger because we are also including the weight
// Set that default to 1
gendata::gendata():hitdata(-1, __GENERATED_SIZE+1){
    SetData(7,1.0);
}

#endif//__CINT__
