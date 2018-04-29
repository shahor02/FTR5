#ifndef DETECTORK_H
#define DETECTORK_H

#include <TNamed.h>
#include "AliLog.h"
#include <TList.h>
#include <TClonesArray.h>
#include <TArrayD.h>
#include "AliExternalTrackParam.h"

//-------------------------------------------------------------------------
// Current support and development: Ruben Shahoyan (Ruben.Shahoyan@cern.ch) 
//-------------------------------------------------------------------------

class KMCLayer;
class KMCCluster;
class TH2F;
class TH1F;
class TGraph;
class TArrayI;

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCProbe : public AliExternalTrackParam {
 public:
  enum {kBitKilled=BIT(14)};
  enum {kNDOF=5};
  enum {kY2=0,kZ2=2,kSnp2=5,kTgl2=9,kPtI2=14};
  enum {kY,kZ,kSnp,kTgl,kPtI};
  //
  KMCProbe();
  KMCProbe(KMCProbe& src);
  KMCProbe& operator=(const KMCProbe& src);
  virtual ~KMCProbe() {}
  virtual   void   Print(Option_t* option = "") const;
  virtual   void   Reset();
  void      ResetCovMat();
  //
  void      Kill(Bool_t v=kTRUE)                          {SetBit(kBitKilled,v);}
  Bool_t    IsKilled()                              const {return TestBit(kBitKilled);}
  //
  Bool_t    CorrectForMeanMaterial(const KMCLayer* lr, Bool_t inward=kTRUE);
  Bool_t    GetXatLabR(Double_t r,Double_t &x, Double_t bz, Int_t dir=0) const;
  Bool_t    PropagateToR(double r, double b, int dir=0);
  //
  void      SetMass(double m=0.14)                      {fMass = m;}
  Double_t  GetMass()                             const {return fMass;}
  Double_t  GetChi2()                             const {return fChi2;}
  void      SetChi2(double chi2)                        {fChi2 = chi2;}
  void      AddChi2(double chi2)                        {fChi2 += chi2;}
  void      SetInnerChecked(int n)                      {fInnerChecked = n;}
  void      SetOuterChecked(int n)                      {fOuterChecked = n;}
  int       GetInnerChecked()                     const {return fInnerChecked;}
  int       GetOuterChecked()                     const {return fOuterChecked;}
  UShort_t  GetNHits()                            const {return fNHits;}
  UShort_t  GetNFakeHits()                        const {return fNHitsFake;}
  UInt_t&   GetHitsPatt()                               {return fHits;}
  UInt_t&   GetFakesPatt()                              {return fFakes;}
  Double_t  GetNormChi2(Bool_t penalize=kFALSE)   const;
  void      AddHit(Int_t lr, double chi2, int clID=-1);
  void      ResetHit(Int_t lr);
  Bool_t    IsHit(Int_t lr)                       const {return (lr<fgNLayers) ? IsWBit(fHits,lr)  : kFALSE;}
  Bool_t    IsHitFake(Int_t lr)                   const {return (lr<fgNLayers) ? IsWBit(fFakes,lr) : kFALSE;}
  Bool_t    PropagateToCluster(KMCCluster* cl, Double_t b);
  // protected: 
  static void   SetWBit(UInt_t &patt,UInt_t bit)               {patt |= 0x1<<bit;}
  static void   ResetWBit(UInt_t &patt,UInt_t bit)             {patt &= ~(0x1<<bit);}
  static Bool_t IsWBit(const UInt_t &patt,const UInt_t bit)    {return patt&(0x1<<bit);}
  static void   SetNLayers(Int_t n)                            {fgNLayers = n;}
  static int    GetNLayers()                                   {return fgNLayers;}
  //
  static Double_t GetMissingHitPenalty()                        {return fgMissingHitPenalty;}
  static void     SetMissingHitPenalty(double p=2.)             {fgMissingHitPenalty = p;}
  //
 public:
  enum {kMaxITSLr=12};
  Float_t fMass;   // mass
  Float_t fChi2;   // total chi2
  UInt_t  fHits;   // pattern on hits (max 32!)
  UInt_t  fFakes;  // pattern of fakes among hits
  UShort_t fNHits;    // total hits
  UShort_t fNHitsFake; // number of fake ITS hits
  UShort_t fInnerChecked; // innermost layer checked
  UShort_t fOuterChecked; // innermost layer checked
  //
  static Int_t    fgNLayers;
  static Double_t fgMissingHitPenalty;
  ClassDef(KMCProbe,1);  
};

//_______________________________________
inline Double_t KMCProbe::GetNormChi2(Bool_t penalize) const
{
  // normalized chi2, penilized for missing hits
  double chi2 = fChi2;
  if (penalize) {
    int nMiss = (fOuterChecked-fInnerChecked+1) - fNHits;
    chi2 = fChi2 + nMiss*fgMissingHitPenalty;
  }
  return chi2/( (fNHits<<1)-kNDOF);
}

//_______________________________________
inline void KMCProbe::AddHit(Int_t lr, double chi2, Int_t clID) {
  // note: lr is active layer ID
  if (lr<0) return;
  fNHits++;
  fChi2 += chi2;
  SetWBit(fHits,lr); 
  if (clID>-1) {
    SetWBit(fFakes,lr);
    fNHitsFake++;
  }
}

//_______________________________________
inline void KMCProbe::ResetHit(Int_t lr) {
  // note: lr is active layer ID
  if (IsWBit(fHits,lr))  {fNHits--;     ResetWBit(fHits,lr);}
  if (IsWBit(fFakes,lr)) {fNHitsFake--; ResetWBit(fFakes,lr);}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCCluster : public TObject {
 public:
  //
  enum {kBitKilled=BIT(14), kDummy = -999};
  KMCCluster(Float_t y=0, Float_t z=0, Float_t x=0, Float_t phi=0, Int_t id=-1) : fY(y),fZ(z),fX(x),fPhi(phi),fID(id) {}
  KMCCluster(KMCCluster &src);
  KMCCluster& operator=(const KMCCluster& src);
  virtual ~KMCCluster() {}
  void Reset() { SetID(kDummy); TObject::Clear();}
  //
  Double_t GetY()    const   {return fY;}
  Double_t GetX()    const   {return fX;}
  Double_t GetZ()    const   {return fZ;}
  Double_t GetPhi()  const   {return fPhi;}
  Int_t    GetID()   const   {return fID;}
  void     SetID(int id)     {fID = id;}
  //
  void    Kill(Bool_t v=kTRUE)          {SetBit(kBitKilled,v);}
  Bool_t  IsKilled()              const {return TestBit(kBitKilled);}
  Float_t fY; 
  Float_t fZ; 
  Float_t fX;
  Float_t fPhi;
  Int_t   fID;
  void Set(Float_t y, Float_t z, Float_t x, Float_t phi) {fY=y; fZ=z; fX=x; fPhi=phi; ResetBit(kBitKilled);}
  virtual void Print(Option_t * = 0) const;
  ClassDef(KMCCluster,1);
};

//____________________________________________________________________________
inline Bool_t KMCProbe::PropagateToCluster(KMCCluster* cl, double b)
{
  // propagate track to cluster frame
  if (!Rotate(cl->GetPhi()) || !PropagateTo(cl->GetX(),b)) {
    AliDebugF(2,"Failed to propager track to cluster at phi=%.3f X=%.3f",cl->GetPhi(),cl->GetX());
    if (AliLog::GetGlobalDebugLevel()>1) Print();
    return kFALSE;
  }
  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCLayer : public TNamed {
public:
  enum {kBitVertex=BIT(15)};
  KMCLayer(char *name);
  Float_t GetRadius()     const {return fR;}
  Float_t GetRadL()       const {return fx2X0;}
  Float_t GetXTimesRho()  const {return fXRho;}
  Float_t GetPhiRes()     const {return fPhiRes;}
  Float_t GetZRes()       const {return fZRes;}
  Float_t GetLayerEff()   const {return fEff;}
  Int_t   GetActiveID()   const {return fActiveID;}
  void    SetActiveID(int id) {fActiveID = id;}
  virtual void  Print(Option_t* option = "") const;
  //
  Bool_t IsActive()     const {return !fIsDead;}
  Bool_t IsDead()       const {return fIsDead;}
  Bool_t IsVertex()     const {return TestBit(kBitVertex);}
  //
  KMCCluster* GetMCCluster()        {return (KMCCluster*)&fClMC;}
  //
  void Reset() {
    if (IsActive()) {fClMC.Reset();}
  }
  //
  //
  Float_t fR; 
  Float_t fx2X0;
  Float_t fXRho;    // x*density
  Float_t fPhiRes; 
  Float_t fZRes;   
  Float_t fEff;
  Bool_t  fIsDead;
  Int_t   fActiveID;   // active layer id
  //
  KMCCluster   fClMC;       // MC cluster (from MS scattered track)
  //
  ClassDef(KMCLayer,1);
};


//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class KMCDetector : public TNamed {
 public:
  enum {kUtilHisto=BIT(14)};
  KMCDetector();
  KMCDetector(char *name,char *title);
  virtual ~KMCDetector();

  void AddLayer(char *name, Float_t radius, Float_t radL, Float_t xrho=0., Float_t phiRes=-1, Float_t zRes=-1, Float_t eff=-1);
  Int_t GetLayerID(Int_t actID) const;

  void PrintLayout(); 
  //  void PlotLayout(Int_t plotDead = kTRUE);
  
  void SetBField(Float_t bfield) {fBFieldG = bfield*10; }
  Float_t GetBField() const {return fBFieldG/10; }

  void SetIntegrationTime(Float_t integrationTime) {fIntegrationTime = integrationTime; }
  Float_t GetIntegrationTime() const { return fIntegrationTime; }

  void SetdNdEtaCent(Int_t dNdEtaCent ) {fdNdEtaCent = dNdEtaCent; }
  Float_t GetdNdEtaCent() const { return fdNdEtaCent; }

  Int_t GetNLayers()          const {return fLayers.GetEntries(); }
  Int_t GetNActiveLayers()    const {return fNActiveLayers; }

  Int_t GetFirstActiveLayer() const {return fFirstActiveLayer;}
  Int_t GetLastActiveLayer() const {return fLastActiveLayer;}

  Int_t GetFirstActiveLayerTracked() const {return fFirstActiveLayerTracked;}
  Int_t GetLastActiveLayerTracked() const {return fLastActiveLayerTracked;}

  // Helper functions
  Double_t ThetaMCS                 ( Double_t mass, Double_t RadLength, Double_t momentum ) const;
  Double_t ProbGoodHit              ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ); 
  Double_t ProbGoodChiSqHit         ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ); 
  Double_t ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ, double confL); 
  Double_t ProbNullChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ, double confL); 

  // Howard W. hit distribution and convolution integral
  Double_t HitDensity        ( Double_t radius )   ;
  Double_t UpcHitDensity     ( Double_t radius )   ;
  Double_t OneEventHitDensity    ( Double_t multiplicity, Double_t radius ) const   ;
  void     CalcDensFactorEta(double eta);
  
  void   ApplyMS(KMCProbe* trc,  double x2x0) const;

  // method to extend AliExternalTrackParam functionality
  Bool_t IsZero(double val, double tol=1e-9) const {return TMath::Abs(val)<tol;}
  TList *GetLayers()                   const {return (TList*)&fLayers;}
  KMCLayer* GetLayer(Int_t i)          const {return (KMCLayer*) fLayers.At(i);}
  KMCLayer* GetActiveLayer(Int_t actID)    const {int pid=GetLayerID(actID); return pid<0 ? 0:GetLayer(pid);}
  KMCLayer* GetLayer(const char* name) const {return (KMCLayer*) fLayers.FindObject(name);}
  KMCProbe* GetProbeTrack()       const {return (KMCProbe*)&fProbe;}
  void   ClassifyLayers();
  void   Reset() { for (int i=fNLayers;i--;) GetLayer(i)->Reset(); }                  

  /*  
  Bool_t SolveSingleTrackViaKalman(Double_t mass, Double_t pt, Double_t eta);
  Bool_t SolveSingleTrackViaKalmanMC(int offset=6);
  Bool_t SolveSingleTrack(Double_t mass, Double_t pt, Double_t eta, TObjArray* sumArr=0, int nMC=10000,int offset=6);
  KMCProbe* KalmanSmooth(int actLr, int actMin,int actMax) const;
  KMCProbe* KalmanSmoothFull(int actLr, int actMin,int actMax) const; //TBD
  void   EliminateUnrelated();
  //

  KMCProbe* PrepareKalmanTrack(double pt, double lambda, double mass, int charge, double phi=0,double x=0,double y=0,double z=0);
  int TransportKalmanTrackWithMS(KMCProbe *probTr);
  Bool_t PropagateToLayer(KMCProbe* trc, KMCLayer* lr, int dir) const;
  Bool_t UpdateTrack(KMCProbe* trc, KMCLayer* lr, KMCCluster* cl, Bool_t goToCluster=kTRUE) const;
 
  //
  Bool_t   GetUseBackground()               const {return fUseBackground;}
  void     SetUseBackground(Bool_t v=kTRUE)       {fUseBackground = v;}
  void     CheckTrackProlongations(KMCProbe *probe, KMCLayer* lr, KMCLayer* lrP);
  void     ResetSearchLimits() {fBgYMin=fBgZMin=1e6; fBgYMax=fBgZMax=-1e6; fNBgLimits=0;}
  void     UpdateSearchLimits(KMCProbe* probe, KMCLayer* lr);
  Int_t    GenBgClusters(KMCLayer* lr);
  Bool_t   NeedToKill(KMCProbe* probe) const;
  Double_t PropagateBack(KMCProbe* trc);
  //
  // definition of reconstructable track
  void     RequireMaxChi2Cl(double cut=25.)           {fMaxChi2Cl = cut>0 ? cut:9; fMaxChi2ClSQ = TMath::Sqrt(fMaxChi2Cl);}
  void     RequireMinITSHits(Int_t n=4)               {fMinITSHits = n;}
  void     RequireMaxNormChi2NDF(double cut=5.)       {fMaxNormChi2NDF = cut>0 ? cut:9;}
  void     RequirePattern(UInt_t *patt, int groups);
  //
  Double_t GetMaxChi2Cl()                      const {return fMaxChi2Cl;}
  Double_t GetMaxNormChi2NDFusterKMC()              const {return fMaxNormChi2NDF;}
  Int_t    GetMinITSHits()                     const {return fMinITSHits;}
  //
  Double_t GetUpdCalls()                       const {return fUpdCalls;}
  TH2F*    GetHMCLrResidRPhi()                 const {return fHMCLrResidRPhi;}
  TH2F*    GetHMCLrResidZ()                    const {return fHMCLrResidZ;}
  TH2F*    GetHMCLrChi2()                      const {return fHMCLrChi2;}
  //
  void     PrintITS(Option_t* opt="") const {for (int i=0;i<=fLastActiveITSLayer;i++) if (!GetLayer(i)->IsDead()) GetLayer(i)->Print(opt);}
  static void SetVtxConstraint(double d=-1, double z=-1) {fgVtxConstraint[0]=d; fgVtxConstraint[1]=z;}
  //
  void CalcHardSearchLimits(double dzv);
  void SetMaxSeedToPropagate(Int_t n=300) {fMaxSeedToPropagate = n;}
  */
 protected:
 
  Int_t fNLayers;        // total number of layers in the model
  Int_t fNActiveLayers;  // number of active layers in the model
  Int_t fFirstActiveLayer;      // id of 1st active layer
  Int_t fLastActiveLayer;       // id of last active layer
  Int_t fFirstActiveLayerTracked;    // id of first active layer really used for tracking of given pt
  Int_t fLastActiveLayerTracked;    // id of last active layer really used for tracking of given pt
  TList fLayers;                // List of layer pointers
  Float_t fBFieldG;             // Magnetic Field in Gauss (set in Tesla)
  Float_t fIntegrationTime;     // electronics integration time

  Int_t fdNdEtaCent;       // Multiplicity
  Double_t fDensFactorEta;                             // density scaling for non-0 eta
  //
  // reconstruction settings
  Double_t fMaxChi2Cl;   // max cluster-track chi2 
  Double_t fMaxNormChi2NDF;// max chi2/NDF to accept
  Int_t    fMinITSHits;  // min ITS hits in track to accept
  //
  KMCProbe fProbe;
  //
  static Double_t fgVtxConstraint[2];  // if both positive, the vertex is used as constraint (accounted in chi2 but not in update)
  ClassDef(KMCDetector,1);
};


#endif
