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
  dds->Print();

  ddd = new TSolSimGEMDigitization (*dds);
  ddd->Print();

  ddgd = new TSolGEMData (1); // 1 hit wonder
  ddgd->SetRun (1000);
  ddgd->SetEvent (0);
  ddgd->SetMomentum (0, TVector3 (0.01, 0.02, 3.00));
  ddgd->SetHitEntrance (0, TVector3 (-0.01, 0.003, 1.55));
  ddgd->SetHitExit (0, TVector3 (-0.02, 0.001, 1.58));
  ddgd->SetHitEntrance (0, TVector3 (-0.022, 0.0008, 1.59));
  ddgd->SetHitEnergy (0, 3.1);
  ddgd->SetHitChamber (0, 0);
  ddgd->SetParticleID (0, 1);
  ddgd->SetParticleType (0, 0);
  ddgd->Print();
  ddgd->PrintHit (0);

  ddd->Digitize (*ddgd);
}
