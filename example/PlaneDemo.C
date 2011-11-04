void PlaneDemo()
{
  TSolGEMChamber y("testchamber","Test chamber");
  y.SetName("testchamber");
  y.Init();

  y.GetPlane(0).Print();
  y.GetPlane(1).Print();
}
