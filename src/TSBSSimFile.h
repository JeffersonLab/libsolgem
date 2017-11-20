/////////////////////////////////////////////////////////////////////
//
//   TSBSSimFile
//
//   Interface to an input file with simulated SoLID data
//
//   Ole Hansen (ole@jlab.org)  December 2011
//
/////////////////////////////////////////////////////////////////////

#ifndef __TSBSSimFile_h
#define __TSBSSimFile_h

#include "THaRunBase.h"
#include "TString.h"
#include "ha_compiledata.h"

class TFile;
class TTree;
class TBranch;
class TSBSSimEvent;

class TSBSSimFile : public THaRunBase {
 public:
  TSBSSimFile(const char* filename, const char* description = "");
  TSBSSimFile(const TSBSSimFile &run);
  virtual ~TSBSSimFile();
  virtual TSBSSimFile &operator=(const THaRunBase &rhs);
  // for ROOT RTTI
  TSBSSimFile() : fROOTFile(0), fTree(0), fEvent(0), fEntry(0) {}

  virtual void  Print( Option_t* opt="" ) const;

  Int_t         Close();
  virtual Int_t Compare( const TObject* obj ) const;
#if ANALYZER_VERSION_CODE >= 67072  // ANALYZER_VERSION(1,6,0)
  const UInt_t* GetEvBuffer() const;
#else
  const  Int_t* GetEvBuffer() const;
#endif
  Int_t         Init();
  const char*   GetFileName() const { return fROOTFileName.Data(); }
  Int_t         Open();
  Int_t         ReadEvent();
  void          SetFileName( const char* name ) { fROOTFileName = name; }

 protected:
  virtual Int_t ReadDatabase() {return 0;}

  TString fROOTFileName;  //  Name of input file
  TFile* fROOTFile;       //! Input ROOT file
  TTree* fTree;           //! Input Tree with simulation data
  TSBSSimEvent* fEvent;   //! Current event

  ULong64_t fNEntries;    //! Number of entries in tree
  ULong64_t fEntry;       //! Current entry number

  ClassDef(TSBSSimFile,1) // Interface to input file with simulated SoLID data
};

#endif
