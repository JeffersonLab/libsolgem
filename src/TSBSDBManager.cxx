#include "TSBSDBManager.h"
#include "TSBSSimDecoder.h"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector2.h"
#include "TRandom3.h"

TSBSDBManager * TSBSDBManager::fManager = NULL;

TSBSDBManager::TSBSDBManager() 
: fErrID(-999), fErrVal(-999.)
{
}
//______________________________________________________________
TSBSDBManager::~TSBSDBManager()
{
}
//______________________________________________________________
void TSBSDBManager::LoadGeneralInfo(const string& fileName)
{  
  //The "sector-plane" concept is not suitable for SBS GEMTrackers. 
  //Instead, "Plane-Module" is introduced. "Plane" means tracking plane and 
  //"Module" means a independent GEM module which is a sub division of the "Plane"
  
    ifstream input(fileName.c_str());
    if (!input.is_open()){
        cout<<"cannot find general information file "<<fileName
            <<". Exiting the program"<<endl;
        exit(0);
    }
    const string prefix = "generalinfo.";

    std::vector<Int_t>* NModule = 0;
    NModule = new vector<Int_t>;
    DBRequest request[] = {
        {"do_map_sector",       &fDoMapSector         , kInt,    0, 1},
        {"self_define_sector",  &fDoSelfDefinedSector , kInt,    0, 1},
        {"sector_mapped",       &fMappedSector        , kInt,    0, 1},
	{"nchamber",            &fNChamber            , kInt,    0, 1},
        {"nsector",             &fNSector             , kInt,    0, 1},
	{"nplane",              &fNGEMPlane           , kInt,    0, 1},
	{"nmodule",             NModule               , kIntV        },
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
	{"g4sbs_detectortype",  &fg4sbsDetectorType   , kInt,    0, 1},
	{"g4sbs_z_specoffset",  &fg4sbsZSpecOffset    , kDouble, 0, 1},
	{"calothreshold",       &fCaloThr             , kDouble, 0, 1},
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
    if( err ) {cout<<"Load DB error"<<endl;exit(2);} 

    
    if(fNGEMPlane!=NModule->size()){cout<<"Check consistency of number of GEM Planes"<<endl;exit(2);}
    int nGEMtot=0;
    for(int i=0;i<fNGEMPlane;i++)
      {
	int nmodule = NModule->at(i);
	fNModule.push_back(nmodule);
	for(int j=0;j<nmodule;j++)
	  {
	    fmPMtoIgem[i][j]=nGEMtot;
	    fmIgemtoPlane[nGEMtot]=i;
	    fmIgemtoModule[nGEMtot]=j;
	    nGEMtot++;
	  }
      }
   
    
    for (int i=0; i<fNSigParticle; i++){
        ostringstream signal_prefix(prefix, ios_base::ate);
        signal_prefix<<"signal"<<i+1<<".";
        
        err = LoadDB(input, signalRequest, signal_prefix.str());
        
        fSigPID.push_back(pid);
        fSigTID.push_back(tid);
	
	if( err ) exit(2); 
    }
        
    for (int i=0; i<GetNGEMPlane(); i++){
      vector<GeoInfo> thisInfo;
      thisInfo.clear();
      fPMGeoInfo[i] = thisInfo;
    }

    delete NModule;

    //fModulesPerChamber = fModulesPerReadOut * fNReadOut;
    
    // fChambersPerCrate = 
    // (TSBSSimDecoder::GetMAXSLOT()/fModulesPerChamber/fNChamber) * fNChamber;
}

void TSBSDBManager::LoadGeoInfo(const string& prefix)
{
  const string& fileName = "db_"+prefix+".dat";
    
  ifstream input(fileName.c_str());
  if (!input.is_open()){
    cout<<"cannot find geometry file "<<fileName
	<<". Exiting the program"<<endl;
    exit(0);
  }
      
  GeoInfo thisGeo;
  
  DBRequest request[] = {
    {"dmag",        &thisGeo.dmag,         kDouble, 0, 1},
    {"d0",          &thisGeo.d0,           kDouble, 0, 1},
    {"xoffset",     &thisGeo.xoffset,      kDouble, 0, 1},
    {"dx",          &thisGeo.dx,           kDouble, 0, 1},
    {"dy",          &thisGeo.dy,           kDouble, 0, 1},
    //{"thetaH",      &thisGeo.thetaH,       kDouble, 0, 1},
    {"thetaV",      &thisGeo.thetaV,       kDouble, 0, 1},
    {"depth",       &thisGeo.depth,        kDouble, 0, 1},
    { 0 }
  };
  
  DBRequest plane_request[] = {
    { "x.stripangle",     &thisGeo.stripangle_u,   kDouble, 0, 1},
    { "x.pitch",          &thisGeo.pitch_u,        kDouble, 0, 1},
    { "y.stripangle",     &thisGeo.stripangle_v,   kDouble, 0, 1},
    { "y.pitch",          &thisGeo.pitch_v,        kDouble, 0, 1},
    { 0 }
  };
  for (int i=0; i<fNGEMPlane; i++){
    map<int, vector<GeoInfo> >::iterator it = fPMGeoInfo.find(i);
    if (it == fPMGeoInfo.end()) { cout<<"unexpected GEM Plane "<<i<<endl; }
    
    for (int j=0; j<fNModule[i]; j++){
      ostringstream plane_prefix(prefix, ios_base::ate);
      int idx = j;
      plane_prefix<<".plane"<<i<<".module"<<j<<".";
      
      int err = LoadDB(input, request, plane_prefix.str());
      if( err ) exit(2);
     
      err = LoadDB(input, plane_request, plane_prefix.str());
      if (err) exit(2);
      
      fPMGeoInfo[i].push_back(thisGeo);
    }
  }
}



//______________________________________________________________
string TSBSDBManager::FindKey( ifstream& inp, const string& key )
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
bool TSBSDBManager::CheckIndex(int i, int j, int k)//(plane, module, readoutAxis)
{
    if (i >= fNChamber || i < 0){
        cout<<"invalid chamber ID requested: "<<i<<endl;
        return false;
    }
    else if(j>=fNModule[i]|| j<0){
      cout<<"invalid module id requested: "<<j<<endl;
      return false;
    }
    else if (k >= fNReadOut || k < 0){
        cout<<"invalid readout id requested: "<<k<<endl;
    }
    return true;
}
//_________________________________________________________________
int TSBSDBManager::LoadDB( ifstream& inp, DBRequest* request, const string& prefix )
{
  DBRequest* item = request;
  while( item->name ) {
    ostringstream sn(prefix, ios_base::ate);
    sn << item->name;
    const string& key = sn.str();
    string val = FindKey(inp,key);
    Int_t tempval;
    if( !val.empty() ) {
      istringstream sv(val);
      switch(item->type){
        case kDouble:
          sv >> *((double*)item->var);
          break;
        case kInt:
          sv >> *((Int_t*)item->var);
          break;
        case kIntV:
	  while(1){
	    if(!sv.good()) break;
	    sv >> tempval;
            ((std::vector<Int_t>*)item->var)->push_back(tempval);
	  }
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
const int & TSBSDBManager::GetSigPID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigPID[i];
}
//______________________________________________________________________
const int & TSBSDBManager::GetSigTID(unsigned int i)
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigTID[i];
}

//______________________________________________________________________
const double & TSBSDBManager::GetDMag(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).dmag;
}
//______________________________________________________________________
const double & TSBSDBManager::GetD0(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).d0;
}
//______________________________________________________________________
const double & TSBSDBManager::GetXOffset(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).xoffset;
}
//______________________________________________________________________
const double & TSBSDBManager::GetDX(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).dx;
}
//______________________________________________________________________
const double & TSBSDBManager::GetDY(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).dy;
}
//______________________________________________________________________
// const double & TSBSDBManager::GetThetaH(int i, int j)
// {
//     if (!CheckIndex(i, j)) return fErrVal;
//     return fGeoInfo[j].at(i).thetaH;
// }
//______________________________________________________________________
const double & TSBSDBManager::GetThetaV(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).thetaV;
}
//_________________________________________________________________________
const double & TSBSDBManager::GetStripAngle(int i, int j, int k)
{
  if (!CheckIndex(i, j, k)) return fErrVal;
  if (k == 0) return fPMGeoInfo[i].at(j).stripangle_u;
  else return fPMGeoInfo[i].at(j).stripangle_u;
}
//_________________________________________________________________________
//const double & TSBSDBManager::GetPitch(int i, int j, int k)
//{
//    if (!CheckIndex(i, j, k)) return fErrVal;
//    if (k == 0) return fGeoInfo[j].at(i).pitch_u;
//    else return fGeoInfo[j].at(i).pitch_u;
//}



int TSBSDBManager::GetModuleIDFromPos(int iplane, double x, double y)
{
  if (!CheckIndex(iplane)) return fErrVal;
  
  int module = -1;
  for(int k = 0; k<fPMGeoInfo[iplane].size(); k++){
    if(fPMGeoInfo[iplane].at(k).xoffset-fPMGeoInfo[iplane].at(k).dx/2.0<=x && 
       x<=fPMGeoInfo[iplane].at(k).xoffset+fPMGeoInfo[iplane].at(k).dx/2.0)
      {
        module = k;
      }
  }

  return module;
}

//__________________________________________________________________________

double TSBSDBManager::GetPosFromModuleStrip(int iproj, int iplane, int imodule, int istrip)
{
  if (!CheckIndex(iplane, imodule)) return fErrVal;

  double pos = fErrVal;
  if(iproj==0){
    pos = fPMGeoInfo[iplane].at(imodule).pitch_u*istrip
         -fPMGeoInfo[iplane].at(imodule).dx/2.0
         +fPMGeoInfo[iplane].at(imodule).xoffset;
    }
  
  if(iproj==1){
    pos = fPMGeoInfo[iplane].at(imodule).pitch_v*istrip
      -fPMGeoInfo[iplane].at(imodule).dy/2.0;
  }
  
  //cout << " " << pos << endl;
  return( pos );
}







