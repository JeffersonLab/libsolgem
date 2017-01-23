#ifndef __TSBSGEANT4FILE_H
#define __TSBSGEANT4FILE_H

// Put prototypes here first so that it doens't freak out
#ifndef __CINT__
#include "evioUtil.hxx"
#endif//__CINT__

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "g4sbs_tree.h"
#include "TSolGEMData.h"
#include "TSolDBManager.h"

#define __DEFAULT_DATA_SIZE 32

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing hit data
//
// Stores an arbitrary double data in dynamically allocated
// arrays.  Allows us to add in data as we get it and then check
// to make sure all entries in the array are filled
// // ___________________________________________________________ //
// // hit_data: {GEM plane, E deposited, X_RO_x, X_RO_y, X_RO_z, 
// //            X_in_x, X_in_y, X_in_z, tmin,
// //            X_out_x, X_out_y, X_out_z, tmax,
// //            type (prim, second), x_vtx, y_vtx, z_vtx,
// //            (???), Particle ID, (???), 
// //            px, py, pz}
// // the strucutre of the data array is identical to the structure 
// // of the hitdata array defined in TSolEVIOFile class

class g4sbshitdata {
    public:
        //Default constructor. The size may depend on the data we examine (GEM, CDET, ECal)
        //TODO: at some point, include CDET and ECal hits
  	g4sbshitdata( int detid, unsigned int size = __DEFAULT_DATA_SIZE );
	virtual ~g4sbshitdata();
	
	//Get detector ID
	int     GetDetID() const { return fDetID;}

	// Get/set one specific element of the data for this hit
	void    SetData( unsigned int, double );
	double  GetData( unsigned int ) const ;
	double *GetData(){ return fData; }//Get all data array 
	
	bool    IsFilled() const ;
	
    protected:
	int     fDetID;//detector ID
	unsigned int     fSize;//data array size;
	long long int fFillbits;
	double *fData;//data array: See in .cxx the sequence of this data array for g4sbs GEMs
};

////////////////////////////////////////////////////////////////////////////
// Auxilliary class for storing generated track data
// // ___________________________________________________________ //
// //gendata: {PID, px, py, pz, x_vtx, y_vtx, z_vtx, weight};
// // the strucutre of this data array is identical to the structure 
// // of the gendata array defined in TSolEVIOFile class

class g4sbsgendata : public g4sbshitdata {
    public:
	g4sbsgendata();
	~g4sbsgendata(){;}
	
	int	GetPID() const { return IsFilled()? (int) fData[0] : -1e9; }//G4 particle ID
	double  GetWeight() const { return fData[7]; }//cross section
	TVector3 GetP() const { return IsFilled()? TVector3(fData[1], fData[2], fData[3]) : TVector3(-1e9, -1e9, -1e9 ); }//Track momentum 3-vector
	TVector3 GetV() const { return IsFilled()? TVector3(fData[4], fData[5], fData[6]) : TVector3(-1e9, -1e9, -1e9 ); }//Track vtx 3-vector
};


///////////////////////////////////////////////////////////////////////////////


class TSBSGeant4File {

 public:
  //constructor may be inputed a data file to input some of the paramaters used by this class
  //NB: if the second file path does not select a valid file, default parameters will be used.
  TSBSGeant4File();// Default constructor
  TSBSGeant4File( const char *name);// Constructor with input file name: recommanded
  virtual ~TSBSGeant4File();// Default destructor

  void ReadGasData(const char* filename); // NB: See comment lines 128-129 
  
  // Standard getters and setters
  void  SetFilename( const char *name );
  void  SetSource( Int_t i ) { fSource = i; }
  void  Clear();
  Int_t Open();
  Int_t Close();
  
  const char* GetFileName() const { return fFilename; }
  Int_t GetSource() const { return fSource; }
  
  // This is actually where the data is read: 
  Int_t ReadNextEvent();
  
  //return the size of the hit arrays
  UInt_t GetNData() const { return fg4sbsHitData.size(); }
  UInt_t GetNGen() const { return fg4sbsGenData.size(); }
  
  UInt_t GetEvNum() const { return fEvNum; }
  
  //get one hit from the hit data arrays
  g4sbshitdata *GetHitData(Int_t i) const { return fg4sbsHitData[i]; }
  g4sbsgendata *GetGenData(Int_t i) const { return fg4sbsGenData[i]; }
  
  //get GEM data
  TSolGEMData *GetGEMData();
  void GetGEMData(TSolGEMData* gd);
  
  double FindGasRange(double p); // NB: See comment lines 128-129 
  
 private:
  // Members
  char  fFilename[255];
  TFile *fFile;
  g4sbs_tree *fTree;// needed to easily unfold root file data
  Int_t fSource;   // User-defined source ID (e.g. MC run number)  // Temp: Do we use that ?
  //double fZSpecOffset; // Offset with which the GEM hits are registered in g4sbs for GEP.
  
  // NB: 2017/01/16: The use of electron range in ionized gas 2017/01/18: Not anymore...
  // is now deprecated due to the addition of X_in and X_out in the g4sbs data on my side
  /* // These two parameters are used to calculate the range in the gas  */
  /* // for very low momentum particles (electrons). */
  /* // This avoids to calculate stupid values for the particle position  */
  /* // at the exit of the GEM ionizable gas (since it is not included in g4sbs output). */
  /* // -> Shall it be ? */
  // This is not necessary for TSolEVIOFile as the hit exit is included in evio files.
  //char fgasdatafile[255];
  vector<double> feMom;
  vector<double> fgasErange;
  
  //hit data arrays
  vector<g4sbshitdata *> fg4sbsHitData;
  vector<g4sbsgendata *> fg4sbsGenData;
  
  unsigned int fEvNum;// global event incrementer

  TSolDBManager *fManager;
};



#endif//__TSOLSIMG4SBSFILE_H
