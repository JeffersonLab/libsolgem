#ifndef __TSOLSIMG4SBSFILE_H
#define __TSOLSIMG4SBSFILE_H

// Put prototypes here first so that it doens't freak out
// over the hidden code


/* #ifdef  __CINT__ */
/* namespace evio { */
/*     class evioFileChannel; */
/*     class evioDOMNodeList; */
/* } */
/* #endif//__CINT__ */

/* // Hide these from the ROOT interpreter */
/* // we don't need them anyways */
#ifndef __CINT__
//#include "evioFileChannel.hxx"
#include "evioUtil.hxx"
#endif//__CINT__

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "g4sbs_gep_tree_with_spin.h"

#include "TSolGEMData.h"


#define __DEFAULT_DATA_SIZE 32

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing hit data
//
// Stores an arbitrary double data in dynamically allocated
// arrays.  Allows us to add in data as we get it and then check
// to make sure all entries in the array are filled

class g4sbshitdata {
    public:
	g4sbshitdata( int detid, unsigned int size = __DEFAULT_DATA_SIZE );
	virtual ~g4sbshitdata();

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

class g4sbsgendata : public g4sbshitdata {
    public:
	g4sbsgendata();
	~g4sbsgendata(){;}
	
	int	GetPID() const { return IsFilled()? (int) fData[0] : -1e9; }
	double  GetWeight() const { return fData[7]; }
	TVector3 GetP() const { return IsFilled()? TVector3(fData[1], fData[2], fData[3]) : TVector3(-1e9, -1e9, -1e9 ); }
	TVector3 GetV() const { return IsFilled()? TVector3(fData[4], fData[5], fData[6]) : TVector3(-1e9, -1e9, -1e9 ); }
	
	/* int	GetPID() const { return IsFilled()? (int) fData[0] : -1e9; } */
	/* int	GetMID() const { return IsFilled()? (int) fData[1] : -1e9; } */
	/* int	GetTRID() const { return IsFilled()? (int) fData[2] : -1e9; } */
	/* //TVector3 GetP() const { return IsFilled()? TVector3(fData[1], fData[2], fData[3]) : TVector3(-1e9, -1e9, -1e9 ); } */
	/* double GetP() const { return IsFilled()? fData[3] : -1e9; } */
	/* double GetBeta() const { return IsFilled()? fData[4] : -1e9; } */
	/* TVector2 GetDXDY() const { return IsFilled()? TVector2(fData[5], fData[6]) : TVector2(-1e9, -1e9 ); } */
	/* TVector2 GetTrueXY() const { return IsFilled()? TVector2(fData[7], fData[8]) : TVector2(-1e9, -1e9 ); } */
	/* TVector3 GetV() const { return IsFilled()? TVector3(fData[9], fData[10], fData[11]) : TVector3(-1e9, -1e9, -1e9 ); } */
	/* TVector3 GetS() const { return IsFilled()? TVector3(fData[12], fData[13], fData[14]) : TVector3(-1e9, -1e9, -1e9 ); } */
	/* double  GetWeight() const { return fData[15]; } */
};


///////////////////////////////////////////////////////////////////////////////


class TSolSimG4SBSFile {

    public:
	TSolSimG4SBSFile();
	TSolSimG4SBSFile( const char *name );
	virtual ~TSolSimG4SBSFile();

	void  SetFilename( const char *name );
        void  SetSource( Int_t i ) { fSource = i; }
	void  Clear();
	Int_t Open();
	Int_t Close();

	const char* GetFileName() const { return fFilename; }
        Int_t GetSource() const { return fSource; }
	Int_t ReadNextEvent();
	/* void  ExtractDetIDs( evio::evioDOMNodeList *, int );// to be replaced */
	/* void  BuildData( evio::evioDOMNodeList * );// to be replaced */
	/* void  BuildGenerated( evio::evioDOMNodeList * );// to be replaced */
	//	void  AddDatum(int crate, int slot, int chan, double data );

	UInt_t GetNData() const { return fg4sbsHitData.size(); }
	UInt_t GetNGen() const { return fg4sbsGenData.size(); }

	UInt_t GetEvNum() const { return fEvNum; }
	
	g4sbshitdata *GetHitData(Int_t i) const { return fg4sbsHitData[i]; }
	g4sbsgendata *GetGenData(Int_t i) const { return fg4sbsGenData[i]; }

	TSolGEMData *GetGEMData();
        void GetGEMData(TSolGEMData* gd);

    private:
	char  fFilename[255];
	TFile *fFile;
	//TChain *fChain;
	g4sbs_gep_tree_with_spin *fTree;
	Int_t fSource;   // User-defined source ID (e.g. MC run number)

	vector<g4sbshitdata *> fg4sbsHitData;
	vector<g4sbsgendata *> fg4sbsGenData;

	unsigned int fEvNum;
};



#endif//__TSOLSIMG4SBSFILE_H
