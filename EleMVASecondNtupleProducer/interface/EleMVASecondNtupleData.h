#ifndef EleMVASecondNtupleData_h
#define EleMVASecondNtupleData_h

#include "NtuTool/Common/interface/TreeWrapper.h"

class EleMVASecondNtupleData: public virtual TreeWrapper {

 public:

  EleMVASecondNtupleData() {
  }
  virtual ~EleMVASecondNtupleData() {
  }

  void initTree() {
    treeName = "EleMVAsecondTree";
    
    // General event variables
    setBranch("evtNumber", &evtNumber, "evtNumber/I", &b_evtNumber);
    setBranch("evtWeight", &evtWeight, "evtWeight/I", &b_evtWeight);
    setBranch("tightEvent", &tightEvent, "tightEvent/O", &b_tightEvent);
    
    setBranch("nGenB", &nGenB, "nGenB/I", &b_nGenB);

    setBranch("JPsiMuHltBit", &JPsiMuHltBit, "JPsiMuHltBit/O", &b_JPsiMuHltBit);
    setBranch("JPsiTrkTrkHltBit", &JPsiTrkTrkHltBit, "JPsiTrkTrkHltBit/O", &b_JPsiTrkTrkHltBit);
    setBranch("JPsiTrkHltBit", &JPsiTrkHltBit, "JPsiTrkHltBit/O", &b_JPsiTrkHltBit);

    setBranch("iPV", &iPV, "iPV/I", &b_iPV);
    
    // Signal-side variables
    setBranch("BPt", &BPt, "BPt/F", &b_BPt);
    setBranch("BEta", &BEta, "BEta/F", &b_BEta);
    setBranch("BPhi", &BPhi, "BPhi/F", &b_BPhi);
    setBranch("BMass", &BMass, "BMass/F", &b_BMass);
    
    setBranch("BLxy", &BLxy, "BLxy/F", &b_BLxy);
    setBranch("BCt2DBS", &BCt2DBS, "BCt2DBS/F", &b_BCt2DBS);
    
    setBranch("BCt2DPV", &BCt2DPV, "BCt2DPV/F", &b_BCt2DPV);
    setBranch("BCt2DPVErr", &BCt2DPVErr, "BCt2DPVErr/F", &b_BCt2DPVErr);
    setBranch("BCt2DPVSigmaUnit", &BCt2DPVSigmaUnit, "BCt2DPVSigmaUnit/F", &b_BCt2DPVSigmaUnit);
    
    setBranch("BCt3DPV", &BCt3DPV, "BCt3DPV/F", &b_BCt3DPV);
    setBranch("BCt3DPVErr", &BCt3DPVErr, "BCt3DPVErr/F", &b_BCt3DPVErr);
    setBranch("BCt3DPVSigmaUnit", &BCt3DPVSigmaUnit, "BCt3DPVSigmaUnit/F", &b_BCt3DPVSigmaUnit);
    
    setBranch("BiSV", &BiSV, "BiSV/I", &b_BiSV);
    
    setBranch("BidGen", &BidGen, "BidGen/I", &b_BidGen);
    
    // Opposite-side variables
    setBranch("elePt", &elePt, "elePt/F", &b_elePt);
    setBranch("eleEta", &eleEta, "eleEta/F", &b_eleEta);
    setBranch("elePhi", &elePhi, "elePhi/F", &b_elePhi);
    
  }
  
  // Second ntuple data
  int evtNumber, evtWeight; bool tightEvent;
  int nGenB;
  bool JPsiMuHltBit, JPsiTrkTrkHltBit, JPsiTrkHltBit;
  int iPV;
  
  float BPt, BEta, BPhi, BMass;
  float BLxy, BCt2DBS;
  float BCt2DPV, BCt2DPVErr, BCt2DPVSigmaUnit;
  float BCt3DPV, BCt3DPVErr, BCt3DPVSigmaUnit;
  int BiSV;
  int BidGen;
  
  float elePt, eleEta, elePhi;

  // Second ntuple branches
  TBranch *b_evtNumber, *b_evtWeight, *b_tightEvent;
  TBranch *b_nGenB;
  TBranch *b_JPsiMuHltBit, *b_JPsiTrkTrkHltBit, *b_JPsiTrkHltBit;
  TBranch *b_iPV;
  
  TBranch *b_BPt, *b_BEta, *b_BPhi, *b_BMass;
  TBranch *b_BLxy, *b_BCt2DBS;
  TBranch *b_BCt2DPV, *b_BCt2DPVErr, *b_BCt2DPVSigmaUnit;
  TBranch *b_BCt3DPV, *b_BCt3DPVErr, *b_BCt3DPVSigmaUnit;
  TBranch *b_BiSV;
  TBranch *b_BidGen;
  
  TBranch *b_elePt, *b_eleEta, *b_elePhi;
  
 private:

  EleMVASecondNtupleData           ( const EleMVASecondNtupleData& a );
  EleMVASecondNtupleData& operator=( const EleMVASecondNtupleData& a );

};

#endif

