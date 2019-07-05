#define UTIL_USE FULL

#include "PDAnalysis/OSElectronFlavorTagger/interface/OsElectronMvaTag.h"

#include "TF1.h"

using namespace std;

OSElectronMvaTag::OSElectronMvaTag():
                osElectronTagReader_("!Color:Silent")
,               ssIndex_(-1)
,               pvIndex_(-1)
,               osElectronIndex_(-1)
,               osElectronTrackIndex_(-1)
,               osElectronTagDecision_(0)
,               osElectronTagMvaValue_(-1.)
,               osElectronTagMistagProbRaw_(-1.)
,               osElectronTagMistagProbCalProcess_(-1.)
,               osElectronTagMistagProbCalProcessBuBs_(-1.)
,               wp_(0.21)
,               dzCut_(1.)
,               nElectronsSel_(0)
{}

OSElectronMvaTag::~OSElectronMvaTag() {}

// =====================================================================================
void OSElectronMvaTag::inizializeOsElectronTagVars()
{
    ssIndex_ = -1;
    pvIndex_ = -1;
    osElectronIndex_ = -1;
    osElectronTrackIndex_ = -1;
    osElectronTagDecision_ = 0;
    osElectronTagMvaValue_ = -1;
    osElectronTagMistagProbRaw_ = -1;
    osElectronTagMistagProbCalProcess_ = -1;
    osElectronTagMistagProbCalProcessBuBs_ = -1;
    nElectronsSel_ = 0;
}

void OSElectronMvaTag::setOsElectronMvaCut(float wp = 0.21)
{
    wp_ = wp;
}

void OSElectronMvaTag::setOsElectronDzCut(float dzCut = 1.)
{
    dzCut_ = dzCut;
}

void OSElectronMvaTag::inizializeOSElectronMvaReader(
    TString methodName = "DNNOsElectronHLTJpsiTrkTrk"
,   TString methodPath = ""
    )
{

    if(methodPath == "") methodPath = methodPath_;
    TMVA::PyMethodBase::PyInitialize();
    methodName_ = methodName;

    osElectronTagReader_.AddVariable("elePt", &elePt_);
    osElectronTagReader_.AddVariable("eleEta", &eleEta_);
    osElectronTagReader_.AddVariable("eleDxy", &eleDxy_);
    osElectronTagReader_.AddVariable("eleExy", &eleExy_);
    osElectronTagReader_.AddVariable("eleDz", &eleDz_);
    osElectronTagReader_.AddVariable("eleEz", &eleEz_);
    osElectronTagReader_.AddVariable("eleIDNIV2Val", &eleIDNIV2Val_);
    osElectronTagReader_.AddVariable("eleIDNIV2Cat", &eleIDNIV2Cat_);
    osElectronTagReader_.AddVariable("eleDRB", &eleDRB_);
    osElectronTagReader_.AddVariable("elePFIsoScaled", &elePFIso_);
    osElectronTagReader_.AddVariable("eleConeCleanPt", &eleConeCleanPt_);
    osElectronTagReader_.AddVariable("eleConeCleanPtRel", &eleConeCleanPtRel_);
    osElectronTagReader_.AddVariable("eleConeCleanDr", &eleConeCleanDr_);
    osElectronTagReader_.AddVariable("eleConeCleanEnergyRatio", &eleConeCleanEnergyRatio_);
    osElectronTagReader_.AddVariable("eleConeCleanQ", &eleConeCleanQ_);
    osElectronTagReader_.BookMVA( methodName_, methodPath + "TMVAClassification_" + methodName_ + ".weights.xml" );
}

bool OSElectronMvaTag::inizializeOSElectronCalibration( 
    TString process = "BuJPsiKData2018"
,   TString processBuMC = "BuJPsiKMC2018"
,   TString processBsMC = "BsJPsiPhiMC2018"
,   TString methodPath = ""  
)
{
    if(methodPath == "") methodPath = methodPath_;
    auto *f   = new TFile(methodPath + "OSElectronTaggerCalibration" + process + ".root");
    auto *fBu = new TFile(methodPath + "OSElectronTaggerCalibration" + processBuMC + ".root");
    auto *fBs = new TFile(methodPath + "OSElectronTaggerCalibration" + processBsMC + ".root");

    if(f->IsZombie()){ cout<<"f IsZombie"<<endl;return false; }
    if(fBu->IsZombie()){ cout<<"fBu IsZombie"<<endl;return false; }
    if(fBs->IsZombie()){ cout<<"fBs IsZombie"<<endl;return false; }

    wCalProcess_ = (TF1*)f->Get("osElectronCal");
    wCalBuMC_    = (TF1*)fBu->Get("osElectronCal");
    wCalBsMC_    = (TF1*)fBs->Get("osElectronCal");

    wCalBuBs_ = new TF1("osElectronCalBuBs","[0]-[1]*[2]/[3]+[2]/[3]*x",0.,1.);
    float qs = wCalBsMC_->GetParameter(0);
    float ms = wCalBsMC_->GetParameter(1);
    float qu = wCalBuMC_->GetParameter(0);
    float mu = wCalBuMC_->GetParameter(1);
    wCalBuBs_->SetParameters(qs, qu, ms, mu);

    delete f;
    delete fBu;
    delete fBs;
    return true;  
}

bool OSElectronMvaTag::makeOsElectronTagging(){
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -999; }
//     selectOsElectron();
    if(osElectronIndex_ < 0){ osElectronTagDecision_ = 0; return 1;}
    else osElectronTagDecision_ = -1*trkCharge->at(osElectronTrackIndex_); 

    computeOsElectronTagVariables();
    osElectronTagMvaValue_ = osElectronTagReader_.EvaluateMVA(methodName_);
    osElectronTagMistagProbRaw_ = 1 - osElectronTagMvaValue_;
    osElectronTagMistagProbCalProcess_ = wCalProcess_->Eval(osElectronTagMistagProbRaw_);
    osElectronTagMistagProbCalProcessBuBs_ = wCalBuBs_->Eval(osElectronTagMistagProbCalProcess_);

    return 1;
}

// int OSMuonMvaTag::selectOsMuon(){
//     int iB = ssIndex_;
//     int iPV = pvIndex_;
// 
//     vector <int> tkSsB = tracksFromSV(iB);
//     TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);
//     if(pvIndex_ < 0) pvIndex_ = GetBestPV(iB, tB);
// 
//     int bestMuIndex = -1;
//     float bestMuPt = 2.;
//     int bestMuTrack = -1;
//     nMuonsSel_ = 0;
// 
//     for(int iMuon = 0; iMuon < nMuons; ++iMuon ){
// 
//         int itkmu = muonTrack( iMuon, PDEnumString::muInner );
//         if(itkmu<0) continue;
// 
//         if(std::find(tkSsB.begin(), tkSsB.end(), itkmu) != tkSsB.end()) continue;
// 
//         if(muoPt->at( iMuon ) < 2.) continue;
//         if(fabs(muoEta->at( iMuon )) > 2.4) continue;
//         if(!IsMvaMuon(iMuon, wp_)) continue;
//         if(fabs(dZ(itkmu, iPV)) > dzCut_) continue;
//         if(deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon)) < 0.4) continue;
//         //if(GetMuoPFiso(iMuon) > PFIsoCut_)  continue;
// 
//         nMuonsSel_++;
// 
//         if(muoPt->at( iMuon ) > bestMuPt){
//             bestMuPt = muoPt->at( iMuon );
//             bestMuIndex = iMuon;
//             bestMuTrack = itkmu;
//         }
//     }
// 
//     osMuonIndex_ = bestMuIndex;
//     osMuonTrackIndex_ = bestMuTrack;
// 
//     return bestMuIndex;
// }

void OSElectronMvaTag::computeOsElectronTagVariables()
{
//   // FIXME: to be rewritten
//     int iB = ssIndex_;
//     int iElectron = osElectronIndex_;
//     int iPV = pvIndex_;
// 
//     int itkmu = muonTrack( iMuon, PDEnumString::muInner );
//     TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);
//     vector <int> tkSsB = tracksFromSV(iB);
// 
//     TLorentzVector tConeClean(0.,0.,0.,0.), tMu;
//     tMu.SetPtEtaPhiM(muoPt->at(iMuon),muoEta->at(iMuon),muoPhi->at(iMuon),MassMu);
// 
//     float kappa = 1;
//     float drCone = 0.4;
// 
//     float qConeClean = 0., ptConeClean = 0.;
//     float muoConeCleanNF = 0., muoConeCleanCF = 0;
//     int muoConeCleanNCH = 0;
// 
//     for(int ipf=0; ipf<nPF; ++ipf){
//         float ptpfc = pfcPt->at(ipf);
//         float etapfc = pfcEta->at(ipf);
//         if( deltaR(etapfc, pfcPhi->at(ipf), muoEta->at(iMuon), muoPhi->at(iMuon)) > drCone) continue;
//         if(ptpfc < 0.5) continue;
//         if(fabs(etapfc) > 3.0) continue;
//         if(pfcCharge->at(ipf) == 0) continue;
//         if(std::find(tkSsB.begin(), tkSsB.end(), pfcTrk->at(ipf)) != tkSsB.end()) continue;
//         if(pfcTrk->at(ipf)<0) continue;
//         if(fabs(dZ(pfcTrk->at(ipf), iPV))>=1.0) continue;
//   
//         TLorentzVector a;
//         a.SetPxPyPzE(pfcPx->at(ipf), pfcPy->at(ipf), pfcPz->at(ipf), pfcE->at(ipf));
//         tConeClean += a;
// 
//         qConeClean += pfcCharge->at(ipf) * pow(ptpfc, kappa);
//         ptConeClean += pow(ptpfc, kappa);
// 
//         if(pfcCharge->at(ipf)==0) muoConeCleanNF += pfcE->at(ipf);
//         if(abs(pfcCharge->at(ipf))==1){
//             muoConeCleanNCH++;
//             muoConeCleanCF += pfcE->at(ipf);
//         }
//     }
// 
//     if(ptConeClean != 0) qConeClean /= ptConeClean;
//     else qConeClean = 1;
//     qConeClean *= trkCharge->at(itkmu);
//     if(tConeClean.E()!=0){
//         muoConeCleanCF /= tConeClean.E();
//         muoConeCleanNF /= tConeClean.E();
//     }
// 
//     muoConeCleanQ_ = qConeClean;
//     muoConeCleanPt_ = tConeClean.Pt();
//     muoConeCleanDr_ = deltaR(tConeClean.Eta(), tConeClean.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
//     if(tConeClean.E()!=0) muoConeCleanEnergyRatio_ = muoE->at(iMuon) / tConeClean.E();
//     else muoConeCleanEnergyRatio_ = 1;
// 
//     tConeClean -= tMu;
//     muoConeCleanPtRel_ = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tConeClean.Vect().Unit());
//     tConeClean += tMu; // for IP sign
// 
//     muoPt_ = muoPt->at( iMuon );
//     muoEta_ = muoEta->at( iMuon );
//     muoCharge_ = trkCharge->at(itkmu);
//     muoDxy_ = dSign(itkmu, tConeClean.Px(), tConeClean.Py())*abs(trkDxy->at(itkmu));
//     muoExy_ = trkExy->at(itkmu);
//     muoDz_ = dZ(itkmu, iPV);
//     muoEz_ = trkEz->at(itkmu);
//     muoSoftMvaValue_ = computeMuonMva(iMuon);
//     muoDrB_ = deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon));
//     muoPFIso_ = GetMuoPFiso(iMuon);
    
}
