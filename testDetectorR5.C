
// You have to load the class before ... ;-)
// .L DetectorK.cxx++

//void standardPlots() {
R5Detector* its = 0;
R5Detector* CreateDetector();
void AddESDTrack(AliESDEvent* ev, AliESDtrack* trc);

AliESDEvent* esdEv = 0;

void testDetectorR5(int nev = 10, int ntracks=10) {

  its = CreateDetector();
  esdEv = new AliESDEvent();
  esdEv->CreateStdContent();

  // create dummy vertex for beam diamond
  const double kBeamSig = 50e-4;
  const double beamPos[3] = {0.,0.,0.};
  const double beamSig[3] = {kBeamSig, kBeamSig, 6. };
  AliESDVertex diamond(beamPos,beamSig,"diamond");
  
  
  for (int iev=0;iev<nev;iev++) {
    esdEv->Reset();
    esdEv->SetMagneticField(its->GetBField()*10.);
    esdEv->SetDiamond(&diamond);
    // generate random vertex
    double vx = gRandom->Gaus(diamond.GetX(), diamond.GetXRes());
    double vy = gRandom->Gaus(diamond.GetY(), diamond.GetYRes());
    double vz = gRandom->Gaus(diamond.GetZ(), diamond.GetZRes());
    
    for (int i=0;i<ntracks;i++) {
      //
      double eta = 2.*gRandom->Rndm()-1.; // random eta between -1:1
      double pT = 2.*gRandom->Rndm()+0.2; // random pT between 0.2 and 2.2
      double phi = gRandom->Rndm()*TMath::Pi()*2;
      int charge = gRandom->Rndm()>0.5 ? 1:-1;
      printf("inp %e %e %d %e\n",pT,eta,charge,phi);
      Bool_t res = its->ProcessTrack(pT, eta, 0.139, charge, phi, vx, vy, vz);    
      if (!res) continue; // track is not reconstructed

      // transfer the track to ESDevent
      AliESDtrack* esdTr = (AliESDtrack*)its->GetProbeTrackInwardAsESDTrack();
      // if needed, add fake TPC flags, though this may create problems in AOD filtering
      esdTr->SetStatus(AliESDtrack::kTPCin|AliESDtrack::kTPCout|AliESDtrack::kTPCrefit);     
      // add extra info if needed, for instance, MC label
      esdTr->SetLabel(i);
      its->AddESDTrackToEvent(esdEv, esdTr); // Call this for proper transfer of info to ESD event
      
    }
    // fit vertexTracks
    AliESDUtils::RefitESDVertexTracks(esdEv);
    int ntr = esdEv->GetNumberOfTracks();
    printf("Event %d: ntr = %d\n",iev,ntr);
    printf("Generated vertex: %+e %+e %+e\n",vx,vy,vz);
    esdEv->GetPrimaryVertexTracks()->Print();
    for (int itr=0;itr<ntr;itr++) {
      AliESDtrack* esdtr = esdEv->GetTrack(itr);
      printf("Tr#%2d Pt: %5.2f Eta: %+4.1f Phi: %+5.2f",itr,esdtr->Pt(), esdtr->Eta(), esdtr->Phi());
      // since esd track has can account for 7 ITS layers pattern at mose, I stored
      // the hits pattern info in the TPC data...
      TBits &hits = esdtr->GetTPCClusterMap();
      //TBits &fakes = esdtr->GetTPCSharedMap(); // here we store fakes, but at the moment they are not set
      printf(" Hits: ");
      for (int ilr=0;ilr<its->GetNActiveLayers();ilr++) {
	printf("%c", hits.TestBitNumber(ilr) ? '+':'-');
      }
      printf("\n");
    }
  }
}



R5Detector* CreateDetector()
{
  AliESDtrack::OnlineMode(kTRUE); // to avoid friend track creation
  R5Detector* det = new R5Detector("ALICE","ITS");

  det->SetPropagateToOrigin(kTRUE); // if we want all tracks to be propagated to DCA to 0/0/0.
  
  det->SetBField(1.);
  // new ideal Pixel properties?
  Double_t x0IB     = 0.001;
  Double_t x0OB     = 0.005;
  Double_t xRho     = 0.;
  Double_t resRPhiIB     = 0.0001;
  Double_t resZIB        = 0.0001;
  Double_t resRPhiOB     = 0.0005;
  Double_t resZOB        = 0.0005;
  Double_t eff           = 0.98;
  //
  // select Z span in such a way to have +-1 unit eta coverage for vertex at 2sigmaZ (~12cm) from nominal IP
  // i.e. Zmax >= 12 + R/tan( 2*atan(exp(-1.)) ) = 12 + R/0.851
  
  det->AddLayer((char*)"vertex",  0.0,  0.1, 0, 0); // dummy vertex for matrix calculation

  det->AddLayer((char*)"bpipe",   1.6,  200., 0.0022);
  det->AddLayer((char*)"ddd1",    1.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff); 
  det->AddLayer((char*)"ddd2",    2.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff); 
  det->AddLayer((char*)"ddd3",    3.8,  21.0, x0IB, xRho, resRPhiIB, resZIB,eff);
  det->AddLayer((char*)"ddd3a",   8.0,  21.0, x0IB, xRho, resRPhiOB, resZOB,eff); 
  det->AddLayer((char*)"ddd4",   20.0,  42.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  det->AddLayer((char*)"ddd5",   25.0,  42.0, x0OB, xRho, resRPhiOB, resZOB,eff); 

  //det->AddLayer((char*)"ddd6",  35.0, 80.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  det->AddLayer((char*)"ddd7",   40.0,  80.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  det->AddLayer((char*)"ddd8",   55.0,  80.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  
  //  det->AddLayer((char*)"dddZ",  90., 130., x0OB, xRho, resRPhiOB, resZOB,eff); 
  det->AddLayer((char*)"dddY",   80.0, 130.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  det->AddLayer((char*)"dddX",  100.0, 130.0, x0OB, xRho, resRPhiOB, resZOB,eff); 
  return det;
}

