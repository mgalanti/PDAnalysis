#ifndef OSElectronMvaTag_H
#define OSElectronMvaTag_H

#include <vector>
#include <set>
#include <map>
#include <string>

#include "OsElectronMvaTagVariables.h"
#include "OsElectronMvaTagVariablesComputer.h"
#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"

#include "TString.h"
#include "TGraph.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/PyMethodBase.h"

class OSElectronMvaTag : 
      public virtual MGBaseAnalyzer, 
      public virtual OsElectronMvaTagVariablesComputer
{

public:
    OSElectronMvaTag();
    ~OSElectronMvaTag();

    void    inizializeOsElectronTagVars();

    bool    makeOsElectronTagging();
//     int     selectOSElectron(); // defined in MGSelector class 

    int     getOsElectron(){ return osElectronIndex_; }
    int     getOsElectronTag(){ return osElectronTagDecision_; }
    float   getOsElectronTagMvaValue(){ return osElectronTagMvaValue_; }
    float   getOsElectronTagMistagProbRaw(){ return osElectronTagMistagProbRaw_; }
    float   getOsElectronTagMistagProbCalProcess(){ return osElectronTagMistagProbCalProcess_; }
    float   getOsElectronTagMistagProbCalProcessBuBs(){ return osElectronTagMistagProbCalProcessBuBs_; }

    void    setVtxOsElectronTag(int iB, int iPV) { ssIndex_ = iB; pvIndex_ = iPV;}
    void    setOsElectron(int iEle) {osElectronIndex_ = iEle;}
    void    setOsElectronMvaCut(float wp);
    void    setOsElectronDzCut(float dzCut);
    void    inizializeOSElectronMvaReader(TString, TString);
    bool    inizializeOSElectronCalibration(TString process, TString processBuMC, TString processBsMC, TString methodPath);

    int     getNosElectrons(){ return nElectronsSel_; }

private:
    
    TString methodNameFromWeightName();
    void    computeOsElectronTagVariables();

    TMVA::Reader osElectronTagReader_;
    TString weightsFile_;
    TString methodName_;
    TString methodPath_ = "mvaWeights/OsElectronTag/";

    int ssIndex_;
    int pvIndex_;
    int osElectronIndex_;
    int osElectronTrackIndex_;
    int osElectronTagDecision_;

    float osElectronTagMvaValue_;
    float osElectronTagMistagProbRaw_;
    float osElectronTagMistagProbCalProcess_;
    float osElectronTagMistagProbCalProcessBuBs_;

    float wp_;
    float dzCut_;
    float PFIsoCut_;

    int nElectronsSel_;

    //MVA Variables
    float elePt_;
    float eleEta_;
    float eleDxy_;
    float eleExy_;
    float eleDz_;
    float eleEz_;
    float eleIDNIV2Val_;
    int eleIDNIV2Cat_;
    float eleDRB_;
    float elePFIso_;
    float eleConeCleanPt_;
    float eleConeCleanPtRel_;
    float eleConeCleanDr_;
    float eleConeCleanEnergyRatio_;
    float eleConeCleanQ_;
    float eleCharge_;

    float DUMMY_;

    //MISTAG VARIABLES
    TF1 *wCalProcess_;
    TF1 *wCalBuMC_;
    TF1 *wCalBsMC_;
    TF1 *wCalBuBs_;

};

#endif
