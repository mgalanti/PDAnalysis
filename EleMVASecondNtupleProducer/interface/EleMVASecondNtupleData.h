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
    
    setBranch("eleIdGen", &eleIdGen, "eleIdGen/I", &b_eleIdGen);
    setBranch("eleBMot", &eleBMot, "eleBMot/I", &b_eleBMot);
    
    setBranch("eleIDNIV2Val", &eleIDNIV2Val, "eleIDNIV2Val/F", &b_eleIDNIV2Val);
    setBranch("eleIDIV2Val", &eleIDIV2Val, "eleIDIV2Val/F", &b_eleIDIV2Val);
    setBranch("eleIDHZZV1Val", &eleIDHZZV1Val, "eleIDHZZV1Val/F", &b_eleIDHZZV1Val);
    setBranch("eleIDNIV2Cat", &eleIDNIV2Cat, "eleIDNIV2Cat/I", &b_eleIDNIV2Cat);
    setBranch("eleIDIV2Cat", &eleIDIV2Cat, "eleIDIV2Cat/I", &b_eleIDIV2Cat);
    setBranch("eleIDHZZV1Cat", &eleIDHZZV1Cat, "eleIDHZZV1Cat/I", &b_eleIDHZZV1Cat);
    
    setBranch("eleDxy", &eleDxy, "eleDxy/F", &b_eleDxy);
    setBranch("eleDz", &eleDz, "eleDz/F", &b_eleDz);
    setBranch("eleExy", &eleExy, "eleExy/F", &b_eleExy);
    setBranch("eleEz", &eleEz, "eleEz/F", &b_eleEz);
    
    setBranch("eleDRB", &eleDRB, "eleDRB/F", &b_eleDRB);
    
    setBranch("eleConePt", &eleConePt, "eleConePt/F", &b_eleConePt);
    setBranch("eleConePtRel", &eleConePtRel, "eleConePtRel/F", &b_eleConePtRel);
    setBranch("eleConeDR", &eleConeDR, "eleConeDR/F", &b_eleConeDR);
    setBranch("eleConeEnergyRatio", &eleConeEnergyRatio, "eleConeEnergyRatio/F", &b_eleConeEnergyRatio);
    setBranch("eleConeQ", &eleConeQ, "eleConeQ/F", &b_eleConeQ);
    setBranch("eleConeSize", &eleConeSize, "eleConeSize/F", &b_eleConeSize);
    setBranch("eleConeNF", &eleConeNF, "eleConeNF/F", &b_eleConeNF);
    setBranch("eleConeCF", &eleConeCF, "eleConeCF/F", &b_eleConeCF);
    setBranch("eleConeNCH", &eleConeNCH, "eleConeNCH/F", &b_eleConeNCH);
        
    setBranch("eleConeCleanPt",          &eleConeCleanPt,          "eleConeCleanPt/F",          &b_eleConeCleanPt);
    setBranch("eleConeCleanPtRel",       &eleConeCleanPtRel,       "eleConeCleanPtRel/F",       &b_eleConeCleanPtRel);
    setBranch("eleConeCleanDR",          &eleConeCleanDR,          "eleConeCleanDR/F",          &b_eleConeCleanDR);
    setBranch("eleConeCleanEnergyRatio", &eleConeCleanEnergyRatio, "eleConeCleanEnergyRatio/F", &b_eleConeCleanEnergyRatio);
    setBranch("eleConeCleanQ",           &eleConeCleanQ,           "eleConeCleanQ/F",           &b_eleConeCleanQ);
    setBranch("eleConeCleanSize",        &eleConeCleanSize,        "eleConeCleanSize/F",        &b_eleConeCleanSize);
    setBranch("eleConeCleanNF",          &eleConeCleanNF,          "eleConeCleanNF/F",          &b_eleConeCleanNF);
    setBranch("eleConeCleanCF",          &eleConeCleanCF,          "eleConeCleanCF/F",          &b_eleConeCleanCF);
    setBranch("eleConeCleanNCH",         &eleConeCleanNCH,         "eleConeCleanNCH/F",         &b_eleConeCleanNCH);
        
    // Tagging truth
    setBranch("tagTruth", &tagTruth, "tagTruth/I", &b_tagTruth);
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
  int eleIdGen, eleBMot;
  float eleIDNIV2Val, eleIDIV2Val, eleIDHZZV1Val;
  int eleIDNIV2Cat, eleIDIV2Cat, eleIDHZZV1Cat;
  float eleDxy, eleDz, eleExy, eleEz;
  float eleDRB;
  float eleConePt, eleConePtRel, eleConeDR, eleConeEnergyRatio, eleConeQ, eleConeSize, eleConeNF, eleConeCF, eleConeNCH;
  float eleConeCleanPt, eleConeCleanPtRel, eleConeCleanDR, eleConeCleanEnergyRatio, eleConeCleanQ, eleConeCleanSize, eleConeCleanNF, eleConeCleanCF, eleConeCleanNCH;
  
  int tagTruth;

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
  TBranch *b_eleIdGen, *b_eleBMot;
  TBranch *b_eleIDNIV2Val, *b_eleIDIV2Val, *b_eleIDHZZV1Val;
  TBranch *b_eleIDNIV2Cat, *b_eleIDIV2Cat, *b_eleIDHZZV1Cat;
  TBranch *b_eleDxy, *b_eleDz, *b_eleExy, *b_eleEz;
  TBranch *b_eleDRB;
  TBranch *b_eleConePt, *b_eleConePtRel, *b_eleConeDR, *b_eleConeEnergyRatio, *b_eleConeQ, *b_eleConeSize, *b_eleConeNF, *b_eleConeCF, *b_eleConeNCH;
  TBranch *b_eleConeCleanPt, *b_eleConeCleanPtRel, *b_eleConeCleanDR, *b_eleConeCleanEnergyRatio, *b_eleConeCleanQ, *b_eleConeCleanSize, *b_eleConeCleanNF, *b_eleConeCleanCF, *b_eleConeCleanNCH;
  
  TBranch *b_tagTruth;
  
 private:

  EleMVASecondNtupleData           ( const EleMVASecondNtupleData& a );
  EleMVASecondNtupleData& operator=( const EleMVASecondNtupleData& a );

};

#endif

