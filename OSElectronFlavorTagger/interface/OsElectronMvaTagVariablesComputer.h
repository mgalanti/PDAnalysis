#ifndef OsElectronMvaTagVariablesComputer_h
#define OsElectronMvaTagVariablesComputer_h



#include "TLorentzVector.h"

// #include "PDAnalysis/Ntu/interface/MGBaseTools.h"
#include "PDAnalysis/Ntu/interface/MGRecoTools.h"
#include "OsElectronMvaTagVariables.h"




class OsElectronMvaTagVariablesComputer :
//       public virtual MGBaseTools,
      public virtual MGRecoTools
{
public:
  OsElectronMvaTagVariablesComputer();
  ~OsElectronMvaTagVariablesComputer(){}
  
  void initialize();
  
  inline void setObjectIndexes(const int iEle, const int iB, const int iPV) {iEle_ = iEle; iB_ = iB; iPV_ = iPV; isInitialized_ = true; evtInitialized_ = currentEvt;}
  
  void computeOsElectronMvaTagVariables();
  
  const OsElectronMvaTagVariables getEleTagVars();
  
  int iEle_;
  int iB_;
  int iPV_;  
  
private:

  bool isInitialized_;
  int evtInitialized_;

  OsElectronMvaTagVariables eleTagVars;
};

#endif
