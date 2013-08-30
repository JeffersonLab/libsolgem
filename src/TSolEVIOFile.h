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
#include "evioFileChannel.hxx"
#include "evioUtil.hxx"
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
	virtual ~hitdata();

	int     GetDetID() const { return fDetID;}

	void    SetData( unsigned int, double );
	double  GetData( unsigned int ) const ;
	double *GetData(){ return fData; }

	bool    IsFilled() const ;

    protected:
	int     fDetID;
	unsigned int     fSize;
	long long int fFillbits;
	double *fData;
};

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing generated track data
//

class gendata : public hitdata {
    public:
	gendata();
	~gendata(){;}

	int	GetPID() const { return IsFilled()? (int) fData[0] : -1e9; }
	double  GetWeight() const { return fData[7]; }
	TVector3 GetP() const { return IsFilled()? TVector3(fData[1], fData[2], fData[3]) : TVector3(-1e9, -1e9, -1e9 ); }
	TVector3 GetV() const { return IsFilled()? TVector3(fData[4], fData[5], fData[6]) : TVector3(-1e9, -1e9, -1e9 ); }
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

	const char* GetFileName() const { return fFilename; }

	Int_t ReadNextEvent();
	void  ExtractDetIDs( evio::evioDOMNodeList *, int );
	void  BuildData( evio::evioDOMNodeList * );
	void  BuildGenerated( evio::evioDOMNodeList * );
  //	void  AddDatum(int crate, int slot, int chan, double data );

	UInt_t GetNData() const { return fHitData.size(); }
	UInt_t GetNGen() const { return fGenData.size(); }

	UInt_t GetEvNum() const { return fEvNum; }

	hitdata *GetHitData(Int_t i) const { return fHitData[i]; }
	gendata *GetGenData(Int_t i) const { return fGenData[i]; }

	TSolGEMData *GetGEMData();
        void GetGEMData(TSolGEMData* gd, Int_t source = 0 );

    private:
	char  fFilename[255];
	evio::evioFileChannel *fChan;

	vector<hitdata *> fHitData;
	vector<gendata *> fGenData;

	unsigned int fEvNum;
};



#endif//__TSOLEVIOFILE_H
