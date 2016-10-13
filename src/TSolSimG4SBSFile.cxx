#include "TSolSimG4SBSFile.h"

//#include "evioUtil.hxx"
//#include "evioFileChannel.hxx"

#include "g4sbs_types.h"
//#include "g4sbs_gep_tree_with_spin.h"

#ifndef __CINT__


TSolSimG4SBSFile::TSolSimG4SBSFile() : fFile(0), fSource(0) {
  fFilename[0] = '\0';
}

TSolSimG4SBSFile::TSolSimG4SBSFile(const char *f) : fFile(0), fSource(0) {
  SetFilename(f);
}

TSolSimG4SBSFile::~TSolSimG4SBSFile() {
  Clear();
  delete fFile;
}

void TSolSimG4SBSFile::SetFilename( const char *f ){
  if( !f ) return;
  strcpy( fFilename, f );
}

Int_t TSolSimG4SBSFile::Open(){
    // Return 0 on fail, 1 on success
    if( fFilename[0] == '\0' ){ return 0; }

    //try {
	// SPR - Unclear to me what a good value to use for the 
	// buffer size is.  I picked something that works...
    delete fFile;
    fFile = new TFile(fFilename);//evio::evioFileChannel(fFilename, "r", 1<<24 );
    
    if( !fFile->IsOpen() ){ 
      fprintf(stderr, "%s: File could not be made\n",__PRETTY_FUNCTION__);
      return 0; 
    }
    
    TChain* C1 = (TChain*)fFile->Get("T");

    fTree = new g4sbs_gep_tree_with_spin(C1, false);

    //fFile->open();
    //}  catch (evio::evioException e) {
    // 	// Problem opening
    // 	fprintf(stderr, "%s %s line %d:  %s\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, e.toString().data() );
    // 	delete fFile; fFile = 0;
    // 	return 0;
    // }
    
    fEvNum = -1;
 
    return 1;
}

Int_t TSolSimG4SBSFile::Close(){
    // Return 0 on fail, 1 on success
    Int_t ret = 1;
    //try {
    if( !fFile->IsOpen() ){ return 0; }
    
    fFile->Close();
    // } catch (evio::evioException e) {
    // 	// Problem closing
    // 	fprintf(stderr, "%s\n", e.toString().data() );
    // 	ret = 0;
    // }
    
    delete fFile; fFile = 0;
    return ret;
}

Int_t TSolSimG4SBSFile::ReadNextEvent(){
    // Return 1 on success

    // Channel not open
    if( !fFile->IsOpen() ){ 
	fprintf(stderr, "%s %s line %d Channel not open\n",
		__FILE__,__PRETTY_FUNCTION__,__LINE__ );
	return 0; 
    }

    Clear();
    
    int n_hits = 0;
    bool res = false;
    
    //cout << fEvNum << endl;
    
    fEvNum++;
    
    //cout << fEvNum << endl;
    
    res = fTree->GetEntry(fEvNum);
    
    if( !res ){
      // Don't need to print this out.  Not really an error
#ifdef  DEBUG
      fprintf(stderr, "%s %s line %d: Channel read return is false...  probably end of file\n",
	      __FILE__, __FUNCTION__, __LINE__ );
#endif//DEBUG
      return 0;
    }
    
    double weight = fTree->ev_solang*fTree->ev_sigma;
    int det_id;//0: FT, 1: FPPs
    
    int pid;
    int type;
    int plane;
    double edep;
    double tmin;
    double tmax;
    
    TVector3 Mom;
    double pz;

    TVector3 X_in;
    TVector3 X_out;
    TVector3 X_RO;
    
    TVector3 Vtx;
    
    double hit_data_temp[23];
    double gen_data_temp[8];
    
    for(int i = 0; i<fTree->Harm_FT_hit_nhits; i++){
      det_id = 0;
      
      pid = fTree->Harm_FT_hit_pid->at(i);
      type = fTree->Harm_FT_hit_mid->at(i)+1;
      plane = fTree->Harm_FT_hit_plane->at(i);
      edep = fTree->Harm_FT_hit_edep->at(i)*1.0e3;
      tmin = fTree->Harm_FT_hit_tmin->at(i);
      tmax = fTree->Harm_FT_hit_tmax->at(i);
      
      pz = sqrt( pow(fTree->Harm_FT_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FT_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FT_hit_typ->at(i), 2) + 1.0) );
      
      Mom = TVector3(fTree->Harm_FT_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FT_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      X_in = TVector3(fTree->Harm_FT_hit_tx->at(i)*1.0e3, // in mm
		      fTree->Harm_FT_hit_ty->at(i)*1.0e3, // in mm
		      fTree->Harm_FT_hit_z->at(i) *1.0e3);// in mm
      
      X_out = TVector3(fTree->Harm_FT_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_txp->at(i), // in mm 
		       fTree->Harm_FT_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FT_hit_typ->at(i), // in mm
		       fTree->Harm_FT_hit_z->at(i) *1.0e3+3.0);// in mm
      
      X_RO = TVector3(fTree->Harm_FT_hit_tx->at(i)*1.0e3+9.185*fTree->Harm_FT_hit_txp->at(i), // in mm 
		      fTree->Harm_FT_hit_ty->at(i)*1.0e3+9.185*fTree->Harm_FT_hit_typ->at(i), // in mm 
		      fTree->Harm_FT_hit_z->at(i) *1.0e3+9.185);// in mm
      
      Vtx = TVector3(fTree->Harm_FT_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FT_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FT_hit_vz->at(i)*1.0e3);// in mm
      
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = type;
      hit_data_temp[17] = -1.0e-9;
      hit_data_temp[18] = pid;
      hit_data_temp[19] = -1.0e-9;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id,  g4sbs_data_size(__GEM_TAG)));

      for(int j = 0; j<23; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      
      gen_data_temp[0] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+1] = Mom[k];
	gen_data_temp[k+4] = Vtx[k];
      }
      gen_data_temp[7] = weight;
      
      fg4sbsGenData.push_back(new g4sbsgendata());
      
      for(int j = 0; j<8; j++){
	fg4sbsGenData[n_hits]->SetData(j, gen_data_temp[j]);
      }
      n_hits++;
      
#ifdef  DEBUG
      cout << "detector ID: " << det_id << ", plane: " << plane << endl
	   << "particle ID: " << pid << ", type (1, primary, >1 secondary): " << type << endl
	   << "energy deposit (eV): " << edep << endl;
      cout << "Momentum (MeV): ";
      for(int k = 0; k<3; k++){
	cout << Mom[k] << ", ";
      }
      cout << " norm " << fTree->Harm_FT_hit_p->at(i) << endl
	   << "dpx/dpz = " << fTree->Harm_FT_hit_txp->at(i)
	   << ", dpx/dpz = " << fTree->Harm_FT_hit_typ->at(i)
	   << endl;
      cout << "hit position at drift entrance (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_in[k] << ", ";
      }
      cout << " time : " << tmin << endl;
      cout << "hit position at drift exit (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_out[k] << " ";
      }
      cout << " time : " << tmax << endl;
      cout << "hit position at readout (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_RO[k] << ", ";
      }
      cout << endl;
      cout << "Vertex position (mm): ";
      for(int k = 0; k<3; k++){
	cout << Vtx[k] << ", ";
      }
      cout << endl;
#endif//DEBUG
    }
    
    for(int i = 0; i<fTree->Harm_FPP1_hit_nhits; i++){
      det_id = 1;
      
      pid = fTree->Harm_FPP1_hit_pid->at(i);
      type = fTree->Harm_FPP1_hit_mid->at(i)+1;
      plane = fTree->Harm_FPP1_hit_plane->at(i);
      edep = fTree->Harm_FPP1_hit_edep->at(i)*1.0e3;
      tmin = fTree->Harm_FPP1_hit_tmin->at(i);
      tmax = fTree->Harm_FPP1_hit_tmax->at(i);
      
      pz = sqrt( pow(fTree->Harm_FPP1_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FPP1_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FPP1_hit_typ->at(i), 2) + 1.0) );
      
      Mom = TVector3(fTree->Harm_FPP1_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FPP1_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      X_in = TVector3(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3, // in mm
		      fTree->Harm_FPP1_hit_ty->at(i)*1.0e3, // in mm
		      fTree->Harm_FPP1_hit_z->at(i) *1.0e3);// in mm
      
      X_out = TVector3(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_txp->at(i), 
		       fTree->Harm_FPP1_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP1_hit_typ->at(i), 
		       fTree->Harm_FPP1_hit_z->at(i) *1.0e3+3.0);// in mm
      
      X_RO = TVector3(fTree->Harm_FPP1_hit_tx->at(i)*1.0e3*9.185*fTree->Harm_FPP1_hit_txp->at(i), 
		      fTree->Harm_FPP1_hit_ty->at(i)*1.0e3*9.185*fTree->Harm_FPP1_hit_typ->at(i), 
		      fTree->Harm_FPP1_hit_z->at(i) *1.0e3+9.185);// in mm
      
      Vtx = TVector3(fTree->Harm_FPP1_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP1_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP1_hit_vz->at(i)*1.0e3);// in mm
      
      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = type;
      hit_data_temp[17] = -1.0e-9;
      hit_data_temp[18] = pid;
      hit_data_temp[19] = -1.0e-9;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id,  g4sbs_data_size(__GEM_TAG)));

      for(int j = 0; j<23; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      
      gen_data_temp[0] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+1] = Mom[k];
	gen_data_temp[k+4] = Vtx[k];
      }
      gen_data_temp[7] = weight;
      
      fg4sbsGenData.push_back(new g4sbsgendata());
      
      for(int j = 0; j<8; j++){
	fg4sbsGenData[n_hits]->SetData(j, gen_data_temp[j]);
      }
      n_hits++;
      
#ifdef DEBUG
      cout << "detector ID: " << det_id << ", plane: " << plane << endl
	   << "particle ID: " << pid << ", type (1, primary, >1 secondary): " << type << endl
	   << "energy deposit (MeV): " << edep << endl;
      cout << "Momentum (MeV): ";
      for(int k = 0; k<3; k++){
	cout << Mom[k] << ", ";
      }
      cout << " norm " << p << endl;
      cout << "hit position at drift entrance (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_in[k] << ", ";
      }
      cout << " time : " << tmin << endl;
      cout << "hit position at drift exit (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_out[k] << " ";
      }
      cout << " time : " << tmax << endl;
      cout << "hit position at readout (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_RO[k] << ", ";
      }
      cout << endl;
      cout << "Vertex position (mm): ";
      for(int k = 0; k<3; k++){
	cout << Vtx[k] << ", ";
      }
      cout << endl;

      cout << "detector ID: " << det_id << ", plane: " << plane << endl
	   << "particle ID: " << pid << ", type (1, primary, >1 secondary: " << type << endl
	   << "energy deposit (GeV): " << edep << endl;
      cout << "Momentum (GeV): ";
      Mom.Print();
      cout << endl;
      cout << "hit position at drift entrance (mm): ";
      X_in.Print();
      cout << " time : " << tmin << endl;
      cout << "hit position at drift exit (mm): ";
      X_out.Print();
      cout << " time : " << tmax << endl;
      cout << "hit position at Readout (mm): ";
      X_RO.Print();
      cout << endl;
      cout << "Vertex position (mm): ";
      Vtx.Print();
      cout << endl;
#endif//DEBUG      
    }
    
    for(int i = 0; i<fTree->Harm_FPP2_hit_nhits; i++){
      det_id = 1;
      
      pid = fTree->Harm_FPP2_hit_pid->at(i);
      type = fTree->Harm_FPP2_hit_mid->at(i)+1;
      plane = 5+fTree->Harm_FPP2_hit_plane->at(i);
      edep = fTree->Harm_FPP2_hit_edep->at(i)*1.0e3;
      tmin = fTree->Harm_FPP2_hit_tmin->at(i);
      tmax = fTree->Harm_FPP2_hit_tmax->at(i);
      
      pz = sqrt( pow(fTree->Harm_FPP2_hit_p->at(i), 2)/
		 ( pow(fTree->Harm_FPP2_hit_txp->at(i), 2) + 
		   pow(fTree->Harm_FPP2_hit_typ->at(i), 2) + 1.0) );
      
      Mom = TVector3(fTree->Harm_FPP2_hit_txp->at(i)*pz*1.0e3, // in MeV
		     fTree->Harm_FPP2_hit_typ->at(i)*pz*1.0e3, // in MeV
		     pz*1.0e3);// in MeV
      
      X_in = TVector3(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3, // in mm
		      fTree->Harm_FPP2_hit_ty->at(i)*1.0e3, // in mm
		      fTree->Harm_FPP2_hit_z->at(i) *1.0e3);// in mm
      
      X_out = TVector3(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_txp->at(i), // in mm 
		       fTree->Harm_FPP2_hit_ty->at(i)*1.0e3+3.0*fTree->Harm_FPP2_hit_typ->at(i), // in mm
		       fTree->Harm_FPP2_hit_z->at(i) *1.0e3+3.0);// in mm
      
      X_RO = TVector3(fTree->Harm_FPP2_hit_tx->at(i)*1.0e3+9.185*fTree->Harm_FPP2_hit_txp->at(i), // in mm
		      fTree->Harm_FPP2_hit_ty->at(i)*1.0e3+9.185*fTree->Harm_FPP2_hit_typ->at(i), // in mm
		      fTree->Harm_FPP2_hit_z->at(i) *1.0e3+9.185);// in mm
      
      Vtx = TVector3(fTree->Harm_FPP2_hit_vx->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP2_hit_vy->at(i)*1.0e3, // in mm
		     fTree->Harm_FPP2_hit_vz->at(i)*1.0e3);// in mm

      hit_data_temp[0] = (double)plane;
      hit_data_temp[1] = edep;
      hit_data_temp[8] = tmin;
      hit_data_temp[12] = tmax;
      hit_data_temp[13] = type;
      hit_data_temp[17] = -1.0e-9;
      hit_data_temp[18] = pid;
      hit_data_temp[19] = -1.0e-9;
      for(int k = 0; k<3; k++){
	hit_data_temp[k+2] = X_RO[k];
	hit_data_temp[k+5] = X_in[k];
	hit_data_temp[k+9] = X_out[k];
	hit_data_temp[k+14] = Vtx[k];
	hit_data_temp[k+20] = Mom[k];
      }
      
      fg4sbsHitData.push_back(new g4sbshitdata(det_id,  g4sbs_data_size(__GEM_TAG)));
      
      for(int j = 0; j<23; j++){
	fg4sbsHitData[n_hits]->SetData(j, hit_data_temp[j]);
      }
      
      gen_data_temp[0] = pid;
      for(int k = 0; k<3; k++){
	gen_data_temp[k+1] = Mom[k];
	gen_data_temp[k+4] = Vtx[k];
      }
      gen_data_temp[7] = weight;
      
      fg4sbsGenData.push_back(new g4sbsgendata());
      
      for(int j = 0; j<8; j++){
	fg4sbsGenData[n_hits]->SetData(j, gen_data_temp[j]);
      }
      n_hits++;
      
#ifdef  DEBUG
      cout << "detector ID: " << det_id << ", plane: " << plane << endl
	   << "particle ID: " << pid << ", type (1, primary, >1 secondary): " << type << endl
	   << "energy deposit (eV): " << edep << endl;
      cout << "Momentum (MeV): ";
      for(int k = 0; k<3; k++){
	cout << Mom[k] << ", ";
      }
      cout << " norm " << p << endl;
      cout << "hit position at drift entrance (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_in[k] << ", ";
      }
      cout << " time : " << tmin << endl;
      cout << "hit position at drift exit (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_out[k] << " ";
      }
      cout << " time : " << tmax << endl;
      cout << "hit position at readout (mm): ";
      for(int k = 0; k<3; k++){
	cout << X_RO[k] << ", ";
      }
      cout << endl;
      cout << "Vertex position (mm): ";
      for(int k = 0; k<3; k++){
	cout << Vtx[k] << ", ";
      }
      cout << endl;
#endif//DEBUG          
    }
    
    return 1;
}


void TSolSimG4SBSFile::Clear(){
    // Clear out hit and generated data

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Deleting hits\n",
		__FILE__, __FUNCTION__, __LINE__);
#endif//DEBUG

    unsigned int i;
    for( i = 0; i < fg4sbsHitData.size(); i++ ){
	delete fg4sbsHitData[i];
    }

    for( i = 0; i < fg4sbsGenData.size(); i++ ){
	delete fg4sbsGenData[i];
    }

    fg4sbsHitData.clear();
    fg4sbsGenData.clear();

#ifdef  DEBUG
	fprintf(stderr, "%s %s line %d: Hits deleted\n",
		__FILE__, __FUNCTION__, __LINE__);
#endif//DEBUG

    return;
}


TSolGEMData* TSolSimG4SBSFile::GetGEMData()
{
  // Return TSolGEMData object filled with GEM data of present event.
  // The returned object pointer must be deleted by the caller!

  TSolGEMData* gd = new TSolGEMData();

  GetGEMData(gd);
  return gd;
}

void TSolSimG4SBSFile::GetGEMData(TSolGEMData* gd)
{
  // Pack data into TSolGEMData
   
//    printf("NEXT EVENT ---------------------------\n");

    if( !gd ) return;
    gd->ClearEvent();
    gd->SetSource(fSource);
    gd->SetEvent(fEvNum);

    if (GetNData() == 0) {
      return;
    }
    gd->InitEvent(GetNData());

    g4sbshitdata *h;//, *hs;
    //bool matchedstrip;
    unsigned int i, ngdata = 0;// j,
    for(i=0; i<GetNData(); i++){
      h = GetHitData(i);

      if( h->GetData(1)>0.0 ){
	
	// Chamber IDs are tagged as 
	// xy  where x is the det id (0 for FT, 1 for FPPs) y the plane number, 
	// labeled from 0 to 9 instead of 1 to 10, for convenience reasons:
	// FT: chambers 0-5, FPPs: chambers 10-19
	
	//if( h->GetDetID()%100 == __GEM_DRIFT_ID &&  h->GetData(1)>0.0 ){
	// Vector information
	TVector3 p(h->GetData(20), h->GetData(21), h->GetData(22));
	gd->SetMomentum(ngdata, p);
	
	TVector3 li(h->GetData(5), h->GetData(6), h->GetData(7));
	gd->SetHitEntrance(ngdata, li);
	
	TVector3 lo(h->GetData(9), h->GetData(10), h->GetData(11));
	gd->SetHitExit(ngdata, lo);
	
	// Average over entrance and exit time
	gd->SetHitTime(ngdata, (h->GetData(8)+h->GetData(12))/2.0);
	
	TVector3 vert(h->GetData(14), h->GetData(15), h->GetData(16));
	gd->SetVertex(ngdata, vert);
	
	TVector3 lr(h->GetData(2), h->GetData(3), h->GetData(4));
	gd->SetHitReadout(ngdata, lr);
	// printf("%d %f %f\n", h->GetDetID()/100, li.X(), li.Y()  );
	
	gd->SetHitEnergy(ngdata, h->GetData(1)*1e6 ); // Gives eV
	gd->SetParticleID(ngdata, (UInt_t) h->GetData(18) );
	gd->SetParticleType(ngdata, (UInt_t) h->GetData(13) );
	
	gd->SetHitChamber(ngdata,  h->GetDetID()*10+h->GetData(0)-1);
	
	ngdata++;
      }
    }
    gd->SetNHit(ngdata);
}


///////////////////////////////////////////////////////////////
// hitdata classes 

g4sbshitdata::g4sbshitdata(int detid, unsigned int size ){
  fDetID = detid;
  fData  = new double[size];
  fSize  = size;
  fFillbits = 0;
  
  if( size > sizeof( long long int )*8 ){
    fprintf(stderr, "%s %s line %d:  Error:  Event size too long for bit pattern storage (requested %d, have %ld)\n",
	    __FILE__, __PRETTY_FUNCTION__, __LINE__, size, 
	    sizeof(long long int)*8);
    exit(1);
  }
  
  // There is no value indexed at 0, so we'll just set it to 0 for
  // sanity's sake and not deal with crazy offsets all over
  
  fFillbits |= 1;
  fData[0] = 3.1415927;
}

void g4sbshitdata::SetData(unsigned int idx, double data ){
    if( idx >= fSize ){
	fprintf(stderr, "%s %s line %d:  Error:  index out of range (%d oor of size %d)\n",__FILE__, __PRETTY_FUNCTION__, __LINE__, idx, fSize);
	return;

    }

    fFillbits |= (1<<idx);

    fData[idx] = data;
    return;
}

double g4sbshitdata::GetData(unsigned int idx) const {
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

bool g4sbshitdata::IsFilled() const {
    if( fFillbits == ((1<<fSize) - 1) ){
	return true;
    }

    return false;
}

g4sbshitdata::~g4sbshitdata(){
    delete fData;
}

///////////////////////////////////////////////////////////////
// gendata classes

// Size is 1 bigger because we are also including the weight
// Set that default to 1
g4sbsgendata::g4sbsgendata():g4sbshitdata(-1, __GENERATED_SIZE+1){
    SetData(7,1.0);
}

#endif//__CINT__
