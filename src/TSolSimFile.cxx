//*-- Author :    Ole Hansen (ole@jlab.org)    9-Dec-2011

/////////////////////////////////////////////////////////////////////
//
//   TSolSimFile
//
//   Interface to an input file with simulated SoLID spectrometer data
//
//   Takes raw digitized simulation data from ROOT input file and
//   uses them to fill a TSolSimEvent object. A pointer to the event
//   object is available via GetEvBuffer() for use by the decoder.
//
/////////////////////////////////////////////////////////////////////

#include "TSolSimFile.h"
#include "TSolSimEvent.h"
#include "evio.h"     // for S_SUCCESS

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"

#include <cstring>
#include <libgen.h>    // for POSIX basename()
#include <cstdlib>
#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------
TSolSimFile::TSolSimFile(const char* filename, const char* description) :
  THaRunBase(description), fROOTFileName(filename), fROOTFile(0), fTree(0), 
  fEvent(0), fNEntries(0), fEntry(0)
{
  // Constructor

  // Use default if no file name given
  if( fROOTFileName.IsNull() ) {
    Info( __FUNCTION__, "Using default input file MCdata.root" );
    fROOTFileName = "MCdata.root";
  }
}

//-----------------------------------------------------------------------------
TSolSimFile::TSolSimFile(const TSolSimFile &run)
  : THaRunBase(run), fROOTFileName(run.fROOTFileName), 
    fROOTFile(0), fTree(0), fEvent(0), fNEntries(0), fEntry(0)
{
}

//-----------------------------------------------------------------------------
TSolSimFile& TSolSimFile::operator=(const THaRunBase& rhs)
{
  if (this != &rhs) {
    THaRunBase::operator=(rhs);
    if( rhs.InheritsFrom("TSolSimFile") )
      fROOTFileName = static_cast<const TSolSimFile&>(rhs).fROOTFileName;
    fROOTFile = 0;
    fTree = 0;
    fEvent = 0;
    fNEntries = fEntry = 0;
  }
  return *this;
}

//_____________________________________________________________________________
Int_t TSolSimFile::Compare( const TObject* obj ) const
{
  // Compare two TSolSimFiles. They are different if either their number
  // or their input file name differs.

  if (this == obj) return 0;
  const THaRunBase* rhs = dynamic_cast<const THaRunBase*>(obj);
  if( !rhs ) return -1;
  // operator< compares fNumber
  if( *this < *rhs )       return -1;
  else if( *rhs < *this )  return  1;
  const TSolSimFile* rhsr = dynamic_cast<const TSolSimFile*>(rhs);
  if( !rhsr ) return 0;
  if( fROOTFileName < rhsr->fROOTFileName ) return -1;
  else if( rhsr->fROOTFileName < fROOTFileName ) return 1;
  return 0;
}

//-----------------------------------------------------------------------------
Int_t TSolSimFile::Init()
{
  // Initialize the run. Sets run date, reads run database etc.

  // TODO: get date from MC production file?
  fDate.Set(2012,1,1,0,0,0);
  fAssumeDate = kTRUE;
  fDataSet |= kDate;

  Int_t ret = THaRunBase::Init();
  if( !ret ) {
    char* s = strdup(fROOTFileName);
    fName = basename(s);
    free(s);
    fNumber = 1;
  }
  return ret;
}

//-----------------------------------------------------------------------------
Int_t TSolSimFile::Open()
{
  // Open ROOT input file

  fROOTFile = new TFile(fROOTFileName, "READ", "SoLID MC data");
  if( !fROOTFile || fROOTFile->IsZombie() ) {
    Error( __FUNCTION__, "Cannot open input file %s", fROOTFileName.Data() );
    Close();
    return -1;
  }

  fTree = static_cast<TTree*>( fROOTFile->Get(treeName) );
  if( !fTree ) {
    Error( __FUNCTION__, "Tree %s does not exist in the input file.",
	   treeName );
    Close();
    return -2;
  }

  fTree->SetBranchStatus("*", kFALSE);

  // Set up reading of the event branch
  delete fEvent; fEvent = 0;

  UInt_t found = 0;
  fTree->SetBranchStatus( eventBranchName, kTRUE, &found );
  if( found > 0 ) {
    fTree->SetBranchAddress( eventBranchName, &fEvent );
  }
  else {
    Error( __FUNCTION__, "No event branch \"%s\" in the input tree.",
	   eventBranchName );
    Close();
    return -3;
  }

  fNEntries = fTree->GetEntries();
  fEntry = 0;

  fOpened = kTRUE;
  return 0;
}

//-----------------------------------------------------------------------------
Int_t TSolSimFile::Close()
{
  delete fTree; fTree = 0;
  if (fROOTFile) {
    fROOTFile->Close();
    delete fROOTFile; fROOTFile = 0;
  }
  delete fEvent; fEvent = 0;
  fOpened = kFALSE;
  return 0;
}

//-----------------------------------------------------------------------------
Int_t TSolSimFile::ReadEvent()
{
  // Read next event from ROOT file

  if( fEntry >= fNEntries )
    return EOF;

  Int_t ret;
  if( !IsOpen() ) {
    ret = Open();
    if( ret ) return ret;
  }

  // Clear the event to get rid of anything still hanging around
  if( fEvent ) fEvent->Clear();

  // Read input file
  ret = fTree->GetEntry(fEntry++);
  if( ret == 0 )
    return EOF;
  if( ret < 0 )
    return -128;  // CODA_ERR

  return S_SUCCESS;
}

//-----------------------------------------------------------------------------
const Int_t *TSolSimFile::GetEvBuffer() const
{
  if( !IsOpen() ) return 0;

  return reinterpret_cast<Int_t*>(fEvent);
}

//-----------------------------------------------------------------------------
TSolSimFile::~TSolSimFile()
{
  if( IsOpen() )
    Close();
}

//_____________________________________________________________________________
void TSolSimFile::Print( Option_t* opt ) const
{
  // Print run info and status

  TString sopt(opt);
  bool do_header = sopt.Contains("start", TString::kIgnoreCase);
  if( sopt.IsNull() || do_header ) {
    THaRunBase::Print("STARTINFO");
  }

  if( sopt.IsNull() || sopt.Contains("status", TString::kIgnoreCase) ) {
    cout << "Analyzed events:       " << fNumAnalyzed << endl;
    cout << "Initialized:           " << fIsInit << endl;
    cout << "Opened:                " << fOpened << endl;
  }

}

//-----------------------------------------------------------------------------
ClassImp(TSolSimFile)
