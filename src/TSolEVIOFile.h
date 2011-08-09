#ifndef __TSOLEVIOFILE_H
#define __TSOLEVIOFILE_H

// Put prototypes here first so that it doens't freak out
// over the hidden code

#ifdef  __CINT__
namespace evio {
    class evioFileChannel;
    class evioDOMNodeList;
}
#endif//__CINT__

// Hide these from the ROOT interpreter
// we don't need them anyways
#ifndef __CINT__
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"
#endif//__CINT__

#include "TROOT.h"

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing hit data
//

class hitdata {
    public:
	hitdata(int crate, int slot, unsigned int size);
	virtual ~hitdata();

	void AddDatum( int chan, double data );
	vector <double> GetDataArray( int chan ){ return fData[chan];}

	bool IsFilled();
	int GetSlot(){ return fSlot; }
	int GetCrate(){ return fCrate; }
	unsigned int GetSize(){ return fSize; }

	void Clear();

	unsigned int GetNHits();
    private:
	unsigned int fSize;
	int fSlot;
	int fCrate;

	vector <double> *fData;
};



///////////////////////////////////////////////////////////////////////////////


class TSolEVIOFile {

    public:
	TSolEVIOFile();
	TSolEVIOFile( const char *name );
	virtual ~TSolEVIOFile();

	void  SetFilename( const char *name );
	void  Clear();
	Int_t Open();
	Int_t Close();

	const char* GetFileName() { return fFilename; }

	Int_t ReadNextEvent();
	void  ExtractDetIDs( evio::evioDOMNodeList *, int );
	void  BuildData( evio::evioDOMNodeList *, int, int );
	void  AddDatum(int crate, int slot, int chan, double data );

	UInt_t GetNData(){ return fHitData.size(); }

	UInt_t GetEvNum(){ return fEvNum; }

	hitdata *GetHitData(Int_t i){ return fHitData[i]; }
    private:
	char  fFilename[255];
	evio::evioFileChannel *fChan;

	vector<int> fDetID;
	vector<hitdata *> fHitData;

	unsigned int fEvNum;
};



#endif//__TSOLEVIOFILE_H
