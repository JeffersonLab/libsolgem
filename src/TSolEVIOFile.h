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

#include "TSolGEMData.h"

#define __DEFAULT_DATA_SIZE 32

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing hit data
//
// Stores an arbitrary double data in dynamically allocated
// arrays.  Allows us to add in data as we get it and then check
// to make sure all entries in the array are filled

class hitdata {
    public:
	hitdata( int detid, unsigned int size = __DEFAULT_DATA_SIZE );
	~hitdata();

	int     GetDetID(){ return fDetID;}

	void    SetData( unsigned int, double );
	double  GetData( unsigned int );
	double *GetData(){ return fData; }

	bool    IsFilled();

    private:
	int     fDetID;
	unsigned int     fSize;
	long long int fFillbits;
	double *fData;
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
	void  BuildData( evio::evioDOMNodeList * );
	void  AddDatum(int crate, int slot, int chan, double data );

	UInt_t GetNData(){ return fHitData.size(); }

	UInt_t GetEvNum(){ return fEvNum; }

	hitdata *GetHitData(Int_t i){ return fHitData[i]; }

	TSolGEMData *GetGEMData();

    private:
	char  fFilename[255];
	evio::evioFileChannel *fChan;

	vector<hitdata *> fHitData;

	unsigned int fEvNum;
};



#endif//__TSOLEVIOFILE_H
