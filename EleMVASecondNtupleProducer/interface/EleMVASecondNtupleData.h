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
    
    // Signal-side variables
    setBranch("tightEvent", &tightEvent, "tightEvent/b", &b_tightEvent);
    
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
    
    // Opposite-side variables
    setBranch("elePt", &elePt, "elePt/F", &b_elePt);
    setBranch("eleEta", &eleEta, "eleEta/F", &b_eleEta);
    setBranch("elePhi", &elePhi, "elePhi/F", &b_elePhi);
    
  }
  
  // Second ntuple data
  unsigned short tightEvent;
  float BPt, BEta, BPhi, BMass;
  float BLxy, BCt2DBS;
  float BCt2DPV, BCt2DPVErr, BCt2DPVSigmaUnit;
  float BCt3DPV, BCt3DPVErr, BCt3DPVSigmaUnit;
  
  float elePt, eleEta, elePhi;

  // Second ntuple branches
  TBranch *b_tightEvent;
  TBranch *b_BPt, *b_BEta, *b_BPhi, *b_BMass;
  TBranch *b_BLxy, *b_BCt2DBS;
  TBranch *b_BCt2DPV, *b_BCt2DPVErr, *b_BCt2DPVSigmaUnit;
  TBranch *b_BCt3DPV, *b_BCt3DPVErr, *b_BCt3DPVSigmaUnit;
  
  TBranch *b_elePt, *b_eleEta, *b_elePhi;
  
 private:

  EleMVASecondNtupleData           ( const EleMVASecondNtupleData& a );
  EleMVASecondNtupleData& operator=( const EleMVASecondNtupleData& a );

};

#endif

