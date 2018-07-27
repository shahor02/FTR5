#ifndef DETECTORK_H
#define DETECTORK_H

#include <TRandom3.h>
#include <TNamed.h>
#include "AliLog.h"
#include <TList.h>
#include <TClonesArray.h>
#include <TArrayD.h>
#include <TBits.h>
#include "AliExternalTrackParam.h"
#include "AliKalmanTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDTOFCluster.h"
#include "AliESDTOFHit.h"


//-------------------------------------------------------------------------
// Current support and development: Ruben Shahoyan (Ruben.Shahoyan@cern.ch) 
//-------------------------------------------------------------------------

class R5Layer;
class R5Cluster;
class TH2F;
class TH1F;
class TGraph;
class TArrayI;
class AliCluster;
class TParticle;

//////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------
class R5Cluster : public TObject {
 public:
  //
  enum {kBitKilled=BIT(14), kDummy = -999};
  R5Cluster(Double_t y=0, Double_t z=0, Double_t x=0, Double_t phi=0, Int_t id=-1) : fY(y),fZ(z),fX(x),fPhi(phi),fID(id) {}
  R5Cluster(R5Cluster &src);
  R5Cluster& operator=(const R5Cluster& src);
  virtual ~R5Cluster() {}
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
  Bool_t  IsValid()               const {return fID!=kDummy && !IsKilled();}
  Double_t fY; 
  Double_t fZ; 
  Double_t fX;
  Double_t fPhi;
  Int_t   fID;
  void Set(Double_t y, Double_t z, Double_t x, Double_t phi, int id) {fY=y; fZ=z; fX=x; fPhi=phi; ResetBit(kBitKilled); fID = id;}
  virtual void Print(Option_t * = 0) const;
  ClassDef(R5Cluster,1);
};


//--------------------------------------------------------------------------------------------
class R5Probe : public AliKalmanTrack {
 public:
  enum {kBitKilled=BIT(14)};
  enum {kNDOF=5};
  enum {kY2=0,kZ2=2,kSnp2=5,kTgl2=9,kPtI2=14};
  enum {kY,kZ,kSnp,kTgl,kPtI};
  //
  R5Probe();
  R5Probe(R5Probe& src);
  R5Probe& operator=(const R5Probe& src);
  virtual ~R5Probe() {}
  virtual   void   Print(Option_t* option = "") const;
  virtual   void   Reset();
  void      ResetCovMat();
  void      ResetTrackingInfo();
  //
  void      Kill(Bool_t v=kTRUE)                          {SetBit(kBitKilled,v);}
  Bool_t    IsKilled()                              const {return TestBit(kBitKilled);}
  //
  Bool_t    CorrectForMeanMaterial(const R5Layer* lr, Bool_t inward=kTRUE);
  Bool_t    GetXatLabR(Double_t r,Double_t &x, Double_t bz, Int_t dir=0) const;
  Bool_t    PropagateToR(Double_t r, Double_t b, int dir=0, Bool_t tof = kFALSE);
  //
  void      SetInnerChecked(int n)                      {fInnerChecked = n;}
  void      SetOuterChecked(int n)                      {fOuterChecked = n;}
  int       GetInnerChecked()                     const {return fInnerChecked;}
  int       GetOuterChecked()                     const {return fOuterChecked;}
  UShort_t  GetNHits()                            const {return GetNumberOfClusters();}
  UShort_t  GetNFakeHits()                        const {return fNHitsFake;}
  UInt_t&   GetHitsPatt()                               {return fHits;}
  UInt_t&   GetFakesPatt()                              {return fFakes;}
  Double_t  GetNormChi2(Bool_t penalize=kFALSE)   const;
  void      AddHit(Int_t lr, Double_t chi2, int clID=-1);
  void      ResetHit(Int_t lr);
  Bool_t    IsHit(Int_t lr)                       const {return (lr<fgNLayers) ? IsWBit(fHits,lr)  : kFALSE;}
  Bool_t    IsHitFake(Int_t lr)                   const {return (lr<fgNLayers) ? IsWBit(fFakes,lr) : kFALSE;}
  Bool_t    PropagateToCluster(R5Cluster* cl, Double_t b, bool tof = kFALSE);
  // protected: 
  static void   SetWBit(UInt_t &patt,UInt_t bit)               {patt |= 0x1<<bit;}
  static void   ResetWBit(UInt_t &patt,UInt_t bit)             {patt &= ~(0x1<<bit);}
  static Bool_t IsWBit(const UInt_t &patt,const UInt_t bit)    {return patt&(0x1<<bit);}
  static void   SetNLayers(Int_t n)                            {fgNLayers = n;}
  static int    GetNLayers()                                   {return fgNLayers;}
  //
  static Double_t GetMissingHitPenalty()                        {return fgMissingHitPenalty;}
  static void     SetMissingHitPenalty(Double_t p=2.)             {fgMissingHitPenalty = p;}
  //
  virtual Double_t GetPIDsignal() const { return GetTOFsignal(); }
  void             SetTOFsignal(float s) { fTOFsignal = s; }
  Double_t         GetTOFsignal() const { return fTOFsignal; }
  
 private:
  // dummy methods
  virtual Double_t GetPredictedChi2(const AliCluster *c) const {return -1;}
  virtual Bool_t PropagateTo(Double_t xr, Double_t x0, Double_t rho) {return kFALSE;}
  virtual Bool_t Update(const AliCluster* c, Double_t chi2, Int_t index) {return kFALSE;}
  
 public:
  UInt_t  fHits;   // pattern on hits (max 32!)
  UInt_t  fFakes;  // pattern of fakes among hits
  Short_t fNHitsFake; // number of fake ITS hits
  Short_t fInnerChecked; // innermost layer checked
  Short_t fOuterChecked; // innermost layer checked
  Float_t fTOFsignal;    // TOF signal in ps
  //
  static Int_t    fgNLayers;
  static Double_t fgMissingHitPenalty;
  ClassDef(R5Probe,1);  
};

//_______________________________________
inline Double_t R5Probe::GetNormChi2(Bool_t penalize) const
{
  // normalized chi2, penilized for missing hits
  Double_t chi2 = GetChi2();
  int nh = GetNumberOfClusters();
  if (penalize) {
    int nMiss = (fOuterChecked-fInnerChecked+1) - nh;
    chi2 += nMiss*fgMissingHitPenalty;
  }
  return chi2/( (nh<<1)-kNDOF);
}

//_______________________________________
inline void R5Probe::AddHit(Int_t lr, Double_t chi2, Int_t clID) {
  // note: lr is active layer ID
  if (lr<0) return;
  SetNumberOfClusters( GetNumberOfClusters() + 1);
  SetChi2( GetChi2() + chi2 );
  SetWBit(fHits,lr); 
  if (clID>-1) {
    SetWBit(fFakes,lr);
    fNHitsFake++;
  }
}

//_______________________________________
inline void R5Probe::ResetHit(Int_t lr) {
  // note: lr is active layer ID
  if (IsWBit(fHits,lr))  {SetNumberOfClusters( GetNumberOfClusters() - 1); ResetWBit(fHits,lr);}
  if (IsWBit(fFakes,lr)) {fNHitsFake--; ResetWBit(fFakes,lr);}
}

//////////////////////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
inline Bool_t R5Probe::PropagateToCluster(R5Cluster* cl, Double_t b, bool tof)
{
  // propagate track to cluster frame
  double xyz0[3],xyz1[3];
  if (tof) {
    GetXYZ(xyz0);
  }
  if ( (TMath::Abs(cl->GetPhi() - GetAlpha())>1e-4 &&  !Rotate(cl->GetPhi())) ||
       (TMath::Abs(cl->GetX() - GetX())>1e-4 && !AliExternalTrackParam::PropagateTo(cl->GetX(),b)) ) {
    AliDebugF(2,"Failed to propagate track to cluster at phi=%.3f X=%.3f",cl->GetPhi(),cl->GetX());
    if (AliLog::GetGlobalDebugLevel()>1) Print();
    return kFALSE;
  }
  if (tof) {
    GetXYZ(xyz1);
    double sgn = xyz0[0]*xyz0[0]+xyz0[1]*xyz0[1] < xyz1[0]*xyz1[0]+xyz1[1]*xyz1[1] ? 1 : -1;
    double dx = xyz1[0]-xyz0[0], dy = xyz1[1]-xyz0[1], dz = xyz1[2]-xyz0[2];
    AddTimeStep( sgn*TMath::Sqrt(dx*dx+dy*dy+dz*dz) );
  }
  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class R5Layer : public TNamed {
public:
  enum {kBitVertex=BIT(15)};
  enum {kPosY,kPosZ,kSigY2,kSigZY,kSigZ2,kNPointParam};
  enum {kDiagErr0,kDiagErr1,kDiagTheta};
  R5Layer(char *name);
  Double_t GetRadius()     const {return fR;}
  Double_t GetZMax()       const {return fZMax;}  
  Double_t GetRadL()       const {return fx2X0;}
  Double_t GetXTimesRho()  const {return fXRho;}
  Double_t GetPhiRes()     const {return fPhiRes;}
  Double_t GetZRes()       const {return fZRes;}
  Double_t GetLayerEff()   const {return fEff;}
  Int_t   GetActiveID()   const {return fActiveID;}
  void    SetActiveID(int id) {fActiveID = id;}
  virtual void  Print(Option_t* option = "") const;
  //
  Bool_t IsActive()     const {return !fIsDead;}
  Bool_t IsDead()       const {return fIsDead;}
  Bool_t IsVertex()     const {return TestBit(kBitVertex);}
  //
  Bool_t InZAcceptane(Double_t z) const {return TMath::Abs(z)<fZMax;}

  R5Cluster* GetMCCluster()        {return (R5Cluster*)&fClMC;}
  //
  void Reset() {
    if (IsActive()) {fClMC.Reset();}
    for (int i=kNPointParam;i--;) {
      fExtInward[i] = fExtOutward[i] = fExtComb[i] = -1.;
    }
    fDiagErr[0]=fDiagErr[1]=fDiagErr[2] = -1.;
  }
  //
  void SetExtInward(const R5Probe* probe) {
    fExtInward[kPosY]  = probe->GetY();
    fExtInward[kPosZ]  = probe->GetZ();
    fExtInward[kSigY2] = probe->GetSigmaY2();
    fExtInward[kSigZY] = probe->GetSigmaZY();
    fExtInward[kSigZ2] = probe->GetSigmaZ2();
  }

  void SetExtOutward(const R5Probe* probe) {
    fExtOutward[kPosY]  = probe->GetY();
    fExtOutward[kPosZ]  = probe->GetZ();
    fExtOutward[kSigY2] = probe->GetSigmaY2();
    fExtOutward[kSigZY] = probe->GetSigmaZY();
    fExtOutward[kSigZ2] = probe->GetSigmaZ2();
  }
  //
  Double_t* GetExtInward() const {return (Double_t*)&fExtInward[0];}
  Double_t* GetExtOutward() const {return (Double_t*)&fExtOutward[0];}
  Double_t* GetExtComb()    const {return (Double_t*)&fExtComb[0];}
  static void Diagonalize2x2Matrix(Double_t sigAA, Double_t sigBA, Double_t sigBB, Double_t &sig11, Double_t &sig22, Double_t &theta);
  void CalcExtComb();
  //
  Double_t fR;
  Double_t fZMax;
  Double_t fx2X0;
  Double_t fXRho;    // x*density
  Double_t fPhiRes; 
  Double_t fZRes;   
  Double_t fEff;
  Double_t fExtInward[5]; // estimate from inward propagation
  Double_t fExtOutward[5]; // estimate from outward propagation
  Double_t fExtComb[5];    // combined estimate
  Double_t fDiagErr[3];    // diagonalized errors
  Bool_t  fIsDead;
  Int_t   fActiveID;   // active layer id
  //
  R5Cluster   fClMC;       // MC cluster (from MS scattered track)
  //
  ClassDef(R5Layer,1);
};


//////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------
class R5Detector : public TNamed {
 public:
  enum {kUtilHisto=BIT(14)};
  enum {kInward = -1, kOutward = 1};
  R5Detector();
  R5Detector(char *name,char *title);
  virtual ~R5Detector();

  // main method to check single track
  Bool_t ProcessTrack(Double_t pt, Double_t eta, Double_t mass, int charge, Double_t phi, Double_t x=0.,Double_t y=0., Double_t z=0., Double_t t=0.);
  Bool_t ProcessTrack(const TParticle* part);
  static void AddESDTrackToEvent(AliESDEvent* ev, const AliESDtrack* trc);
  
  void AddLayer(char *name, Double_t radius, Double_t zmax, Double_t radL, Double_t xrho=0., Double_t phiRes=-1, Double_t zRes=-1, Double_t eff=-1);
  Int_t GetLayerID(Int_t actID) const;

  virtual  void Print(const Option_t* opt) const; 
  //  void PlotLayout(Int_t plotDead = kTRUE);
  
  void SetBField(Double_t bfield) {fBFieldG = bfield*10; }
  Double_t GetBField() const {return fBFieldG/10; }

  void SetIntegrationTime(Double_t integrationTime) {fIntegrationTime = integrationTime; }
  Double_t GetIntegrationTime() const { return fIntegrationTime; }

  void SetdNdEtaCent(Int_t dNdEtaCent ) {fdNdEtaCent = dNdEtaCent; }
  Double_t GetdNdEtaCent() const { return fdNdEtaCent; }

  Int_t GetNLayers()          const {return fLayers.GetEntries(); }
  Int_t GetNActiveLayers()    const {return fNActiveLayers; }

  Int_t GetFirstActiveLayer() const {return fFirstActiveLayer;}
  Int_t GetLastActiveLayer() const {return fLastActiveLayer;}

  Int_t GetFirstActiveLayerTracked() const {return fFirstActiveLayerTracked;}
  Int_t GetLastActiveLayerTracked() const {return fLastActiveLayerTracked;}

  void SetTOFResolution(double r=20) { fTOFResolutionPS = r;}
  Double_t GetTOFResolution() const { return fTOFResolutionPS;}

  // Helper functions
  Double_t ThetaMCS                 ( Double_t mass, Double_t RadLength, Double_t momentum ) const;
  Double_t ProbGoodHit              ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ); 
  Double_t ProbGoodChiSqHit         ( Double_t radius, Double_t searchRadiusRPhi, Double_t searchRadiusZ ); 
  Double_t ProbGoodChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ, Double_t confL); 
  Double_t ProbNullChiSqPlusConfHit ( Double_t radius, Double_t leff, Double_t searchRadiusRPhi, Double_t searchRadiusZ, Double_t confL); 

  // Howard W. hit distribution and convolution integral
  Double_t HitDensity        ( Double_t radius )   ;
  Double_t UpcHitDensity     ( Double_t radius )   ;
  Double_t OneEventHitDensity    ( Double_t multiplicity, Double_t radius ) const   ;
  void     CalcDensFactorEta(Double_t eta);
  
  void   ApplyMS(R5Probe* trc,  Double_t x2x0) const;

  // method to extend AliExternalTrackParam functionality
  Bool_t    IsZero(Double_t val, Double_t tol=1e-9) const {return TMath::Abs(val)<tol;}
  TList*    GetLayers()                   const {return (TList*)&fLayers;}
  R5Layer* GetLayer(Int_t i)          const {return (R5Layer*) fLayers.At(i);}
  R5Layer* GetActiveLayer(Int_t actID)    const {int pid=GetLayerID(actID); return pid<0 ? 0:GetLayer(pid);}
  R5Layer* GetLayer(const char* name) const {return (R5Layer*) fLayers.FindObject(name);}
  R5Probe* GetProbeGen()         const {return (R5Probe*)&fProbeInMC0;}
  R5Probe* GetProbeTrackInward() const {return (R5Probe*)&fProbeInward;}
  R5Probe* GetProbeOutMC()         const {return (R5Probe*)&fProbeOutMC;}
  const AliESDtrack* GetProbeTrackInwardAsESDTrack();
  void      ClassifyLayers();
  void      Reset();
  void      SetMaxSnp(Double_t v) { fMaxSnp = v; }
  Double_t  GetMaxSnp() const { return fMaxSnp; }
  void      SetMaxChi2Cl(Double_t cut)  {fMaxChi2Cl = cut>0 ? cut:9;}
  void      SetMinHits(Int_t n=4)     {fMinHits = n;}
  void      SetMaxNormChi2NDF(Double_t cut=5.) {fMaxNormChi2NDF = cut>0 ? cut:9;}

  Double_t GetMaxChi2Cl()                   const {return fMaxChi2Cl;}
  Double_t GetMaxNormChi2NDF()              const {return fMaxNormChi2NDF;}
  Int_t    GetMinHits()                     const {return fMinHits;}
  void     SetPropagateToOrigin(Bool_t v)        {fPropagateToOrigin = v;}
  Bool_t   GetPropagateToOrigin()           const {return fPropagateToOrigin;}
  

  void PrepareKalmanTrack(Double_t pt, Double_t eta, Double_t mass, int charge, Double_t phi=0,Double_t x=0.,Double_t y=0.,Double_t z=0., Double_t t=0.);
  int TransportKalmanTrackWithMS(R5Probe *probTr, Bool_t applyMatCorr=kTRUE);
  Bool_t PropagateToLayer(R5Probe* trc, const R5Layer* lr, int dir) const;
  Bool_t ExtrapolateToR(R5Probe* probe, Double_t r) const;
  Bool_t UpdateTrack(R5Probe* trc, R5Layer* lr, R5Cluster* cl) const;
  Double_t  Chi2ToCluster(const R5Layer* lr, const R5Probe* trc, const R5Cluster* cl) const;
  
  /*  
  Bool_t SolveSingleTrackViaKalman(Double_t mass, Double_t pt, Double_t eta);
  Bool_t SolveSingleTrackViaKalmanMC(int offset=6);
  Bool_t SolveSingleTrack(Double_t mass, Double_t pt, Double_t eta, TObjArray* sumArr=0, int nMC=10000,int offset=6);
  R5Probe* KalmanSmooth(int actLr, int actMin,int actMax) const;
  R5Probe* KalmanSmoothFull(int actLr, int actMin,int actMax) const; //TBD
  void   EliminateUnrelated();
  //

 
  //
  Bool_t   GetUseBackground()               const {return fUseBackground;}
  void     SetUseBackground(Bool_t v=kTRUE)       {fUseBackground = v;}
  void     CheckTrackProlongations(R5Probe *probe, R5Layer* lr, R5Layer* lrP);
  void     ResetSearchLimits() {fBgYMin=fBgZMin=1e6; fBgYMax=fBgZMax=-1e6; fNBgLimits=0;}
  void     UpdateSearchLimits(R5Probe* probe, R5Layer* lr);
  Int_t    GenBgClusters(R5Layer* lr);
  Bool_t   NeedToKill(R5Probe* probe) const;
  Double_t PropagateBack(R5Probe* trc);
  //
  // definition of reconstructable track
  void     RequirePattern(UInt_t *patt, int groups);
  //
  //
  Double_t GetUpdCalls()                       const {return fUpdCalls;}
  TH2F*    GetHMCLrResidRPhi()                 const {return fHMCLrResidRPhi;}
  TH2F*    GetHMCLrResidZ()                    const {return fHMCLrResidZ;}
  TH2F*    GetHMCLrChi2()                      const {return fHMCLrChi2;}
  //
  void     PrintITS(Option_t* opt="") const {for (int i=0;i<=fLastActiveITSLayer;i++) if (!GetLayer(i)->IsDead()) GetLayer(i)->Print(opt);}
  static void SetVtxConstraint(Double_t d=-1, Double_t z=-1) {fgVtxConstraint[0]=d; fgVtxConstraint[1]=z;}
  //
  void CalcHardSearchLimits(Double_t dzv);
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
  Double_t fBFieldG;             // Magnetic Field in Gauss (set in Tesla)
  Double_t fIntegrationTime;     // electronics integration time

  Int_t fdNdEtaCent;       // Multiplicity
  Double_t fDensFactorEta;                             // density scaling for non-0 eta
  //
  // reconstruction settings
  Double_t fMaxChi2Cl;   // max cluster-track chi2 
  Double_t fMaxNormChi2NDF;// max chi2/NDF to accept
  Int_t    fMinHits;  // min ITS hits in track to accept
  Double_t fMaxSnp;      // max allowe snp
  Bool_t   fPropagateToOrigin; // propagate all tracks to DCA to origin
  Double_t fTOFResolutionPS; // TOF resolution in ps
  //
  R5Probe fProbeInMC0; // initially provided probe
  R5Probe fProbeOutMC; // probe propagated to outer radius with material effects
  R5Probe fProbeInward; // probe after innermost update
  AliESDtrack fESDtrack; // place-holder to transfer reconstructed track as ESDtrack
  //
  static Double_t fgVtxConstraint[2];  // if both positive, the vertex is used as constraint (accounted in chi2 but not in update)
  ClassDef(R5Detector,1);
};


#endif
