TSolSpec* dds;
TSolGEMChamber* ddy;
TSolSimGEMDigitization* ddd;
TSolGEMData* ddgd;

void DigDemo()
{
  dds = new TSolSpec ("spectrometer", "SOLiD spectrometer");
  dds->Init();

  ddy = new TSolGEMChamber ("testchamber","Test chamber");
  ddy->SetName("testchamber");
  ddy->Init();
  

  dds->AddGEM (*ddy);

  ddd = new TSolSimGEMDigitization (*dds);
  
  ddgd = new TSolGEMData (1); // 1 hit wonder
  ddgd->SetRun (1000);
  ddgd->SetEvent (0);
  ddgd->SetMomentum (0, TVector3 (0.01, 0.02, 3.00));
  ddgd->SetHitEntrance (0, TVector3 (-0.01, 0.003, 1.55) * 1000.0); // mm
  ddgd->SetHitExit (0, TVector3 (-0.02, 0.001, 1.58) * 1000.0);
  ddgd->SetHitReadout (0, TVector3 (-0.022, 0.0008, 1.59) * 1000.0);
  ddgd->SetHitEnergy (0, 1e3); // eV
  ddgd->SetHitChamber (0, 0);
  ddgd->SetParticleID (0, 1);
  ddgd->SetParticleType (0, 0);
  //  ddgd->Print();
  //  ddgd->PrintHit (0);

  ddd->InitTree (*dds, "digdemo.root");
  ddd->Digitize (*ddgd, *dds);
  ddd->PrintCharges();
  ddd->WriteTree();
  ddd->CloseTree();
}
