#include "TSolDBManager.h"
#include "TSolSimDecoder.h"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector2.h"

TSolDBManager * TSolDBManager::fManager = NULL;

TSolDBManager::TSolDBManager() 
: fErrID(-999), fErrVal(-999.)
{
}
//______________________________________________________________
TSolDBManager::~TSolDBManager()
{
}
//______________________________________________________________
void TSolDBManager::LoadGeneralInfo(const string& fileName)
{
    ifstream input(fileName.c_str());
    if (!input.is_open()){
        cout<<"cannot find general information file "<<fileName
            <<". Exiting the program"<<endl;
        exit(0);
    }
    const string prefix = "generalinfo.";
    DBRequest request[] = {
        {"do_map_sector",       &fDoMapSector         , kInt,    0, 1},
        {"self_define_sector",  &fDoSelfDefinedSector , kInt,    0, 1},
        {"sector_mapped",       &fMappedSector        , kInt,    0, 1},
        {"ntracker",            &fNTracker            , kInt,    0, 1},
        {"nsector",             &fNSector             , kInt,    0, 1},
        {"nreadout",            &fNReadOut            , kInt,    0, 1},
        {"gem_drift_id",        &fGEMDriftID          , kInt,    0, 1},
        {"gem_copper_front_id", &fGEMCopperFrontID    , kInt,    0, 1},
        {"gem_copper_back_id",  &fGEMCopperBackID     , kInt,    0, 1},
        {"gem_strip_id",        &fGEMStripID          , kInt,    0, 1},
        {"faec_id",             &fFAECID              , kInt,    0, 1},
        {"laec_id",             &fLAECID              , kInt,    0, 1},
        {"nsignal",             &fNSigParticle        , kInt,    0, 1},
        {"chan_per_slot",       &fChanPerSlot         , kInt,    0, 1},
        {"modules_per_readout", &fModulesPerReadOut   , kInt,    0, 1},
        {"self_defined_sector", &fDoSelfDefinedSector , kInt,    0, 1},
        {"self_defined_sector", &fDoSelfDefinedSector , kInt,    0, 1},
        {"z_g4sbsspecoffset",   &fZg4sbsSpecOffset    , kDouble, 0, 1},
        {"gasdatafilename",     &fGasDataFilename     , kString, 0, 1},
	{"z0",                  &fgZ0                 , kDouble, 0, 1},
	{"calo_z",              &fgCaloZ              , kDouble, 0, 1},
	{"calo_res",            &fgCaloRes            , kDouble, 0, 1},
	{"docalo",              &fgDoCalo             , kInt,    0, 1},
        { 0 }
    };
    int pid, tid;
    DBRequest signalRequest[] = {
        {"pid",                 &pid,                   kInt, 0, 1},
        {"tid",                 &tid,                   kInt, 0, 1},
        { 0 }
    };
    int err = LoadDB( input, request,  prefix);
    
    if( err ) exit(2); 
    
    for (int i=0; i<fNSigParticle; i++){
        ostringstream signal_prefix(prefix, ios_base::ate);
        signal_prefix<<"signal"<<i+1<<".";
        
        err = LoadDB(input, signalRequest, signal_prefix.str());
        
        fSigPID.push_back(pid);
        fSigTID.push_back(tid);
        
        if( err ) exit(2); 
    }
    
    /*
    for (int i=0; i<fNTracker; i++){
         vector<GeoInfo> thisInfo;
         thisInfo.clear();
         fGeoInfo[i] = thisInfo;
    }
    */
    
    fModulesPerChamber = fModulesPerReadOut * fNReadOut;
    
    fChambersPerCrate  = 
    (TSolSimDecoder::GetMAXSLOT()/fModulesPerChamber/fNTracker) * fNTracker;
}
/*
//______________________________________________________________
void TSolDBManager::LoadGeoInfo(const string& fileName)
{
    ifstream input(fileName.c_str());
    if (!input.is_open()){
        cout<<"cannot find general information file "<<fileName
            <<". Exiting the program"<<endl;
        exit(0);
    }
    
    const string prefix = "gemc.";
    
    GeoInfo thisGeo;
    
    DBRequest request[] = {
        { "r0",               &thisGeo.r0,             kDouble, 0, 1},
        { "r1",               &thisGeo.r1,             kDouble, 0, 1},
        { "phi0",             &thisGeo.phi0,           kDouble, 0, 1},
        { "dphi",             &thisGeo.dphi,           kDouble, 0, 1},
        { "z0",                &thisGeo.z,              kDouble, 0, 1},
        { "depth",            &thisGeo.depth,          kDouble, 0, 1},
        { 0 }
    };
    
    DBRequest plane_request[] = {
        { "x.stripangle",     &thisGeo.stripangle_u,   kDouble, 0, 1},
        { "x.pitch",          &thisGeo.pitch_u,        kDouble, 0, 1},
        { "y.stripangle",     &thisGeo.stripangle_v,   kDouble, 0, 1},
        { "y.pitch",          &thisGeo.pitch_v,        kDouble, 0, 1},
        { 0 }
    };
    
    for (int i=0; i<fNTracker; i++){
        map<int, vector<GeoInfo> >::iterator it = fGeoInfo.find(i);
        if (it == fGeoInfo.end()) { cout<<"unexpected tracker id "<<i<<endl; }
    
        for (int j=0; j<fNSector; j++){
            ostringstream sector_prefix(prefix, ios_base::ate);
            int idx = i*fNSector + j;
            sector_prefix<<"gem"<<idx+1<<".";
            
            int err = LoadDB(input, request, sector_prefix.str());
            if( err ) exit(2);
            
            sector_prefix<<"gem"<<idx+1;
            err = LoadDB(input, plane_request, sector_prefix.str());
            if (err) exit(2);
            
            fGeoInfo[i].push_back(thisGeo);
        }
    }
}
*/
//______________________________________________________________
string TSolDBManager::FindKey( ifstream& inp, const string& key )
{
  static const string empty("");
  string line;
  string::size_type keylen = key.size();
  inp.seekg(0); // could probably be more efficient, but it's fast enough
  while( getline(inp,line) ) {
    if( line.size() <= keylen )
      continue;
    if( line.compare(0,keylen,key) == 0 ) {
      if( keylen < line.size() ) {
	string::size_type pos = line.find_first_not_of(" \t=", keylen);
	if( pos != string::npos )
	  return line.substr(pos);
      }
      break;
    }
  }
  return empty;
}
//_________________________________________________________________________
bool TSolDBManager::CheckIndex(int i, int j, int k)
{
    if (i >= fNTracker || i < 0){
        cout<<"invalid tracker ID requested: "<<i<<endl;
        return false;
    }
    else if (j >= fNSector || j < 0){
        cout<<"invalid sector id requested: "<<j<<endl;
        return false;
    }
    else if (k >= fNReadOut || k < 0){
        cout<<"invalid readout id requested: "<<k<<endl;
    }
    return true;
}
//_________________________________________________________________
int TSolDBManager::LoadDB( ifstream& inp, DBRequest* request, const string& prefix )
{
  DBRequest* item = request;
  while( item->name ) {
    ostringstream sn(prefix, ios_base::ate);
    sn << item->name;
    const string& key = sn.str();
    string val = FindKey(inp,key);
    if( !val.empty() ) {
      istringstream sv(val);
      switch(item->type){
        case kDouble:
          sv >> *((double*)item->var);
          break;
        case kInt:
          sv >> *((Int_t*)item->var);
          break;
        default:
          return 1;
        break;
      }
      if( !sv ) {
	cerr << "Error converting key/value = " << key << "/" << val << endl;
	return 1;
      }
    } else {
      cerr << "key \"" << key << "\" not found" << endl;
      return 2;
    }
    ++item;
  }
  return 0;
}
//_____________________________________________________________________
const int & TSolDBManager::GetSigPID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigPID[i];
}
//______________________________________________________________________
const int & TSolDBManager::GetSigTID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigTID[i];
}
/*
//______________________________________________________________________
const double & TSolDBManager::GetSectorZ(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[i].at(j).z;
}
//_______________________________________________________________________
const double & TSolDBManager::GetSectorRMin(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[i].at(j).r0;
}
//________________________________________________________________________
const double & TSolDBManager::GetSectorRMax(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[i].at(j).r1;
}
//_________________________________________________________________________
const double & TSolDBManager::GetSectorPhiStart(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[i].at(j).phi0;
}
//_________________________________________________________________________
const double & TSolDBManager::GetSectorPhiCover(int i, int j)
{
    if (!CheckIndex(i, j)) return fErrVal;
    return fGeoInfo[i].at(j).dphi;
}
//_________________________________________________________________________
const double & TSolDBManager::GetSectorStripAngle(int i, int j, int k)
{
    if (!CheckIndex(i, j, k)) return fErrVal;
    if (k == 0) return fGeoInfo[i].at(j).stripangle_u;
    else return fGeoInfo[i].at(j).stripangle_u;
}
//_________________________________________________________________________
const double & TSolDBManager::GetSectorPitch(int i, int j, int k)
{
    if (!CheckIndex(i, j, k)) return fErrVal;
    if (k == 0) return fGeoInfo[i].at(j).pitch_u;
    else return fGeoInfo[i].at(j).pitch_u;
}
//__________________________________________________________________________
int TSolDBManager::GetSectorIDFromPos(double& x, double& y, int& itracker)
{
    if (!CheckIndex(itracker)) return fErrVal;
    double thisPhi = atan2(y, x);
    for (unsigned int i=0; i<fGeoInfo[itracker].size(); i++){
        double phiCenter = (fGeoInfo[itracker].at(i).phi0 + 
                            fGeoInfo[itracker].at(i).dphi / 2.) /180. *TMath::Pi();
        phiCenter = TVector2::Phi_mpi_pi(phiCenter);
        double deltaPhi = fabs(TVector2::Phi_mpi_pi(thisPhi - phiCenter));
        if (deltaPhi < fGeoInfo[itracker].at(i).dphi / 2. / 180. *TMath::Pi()){
            return (int)i;
            break;
        }
    }
    return -1;
}
*/







