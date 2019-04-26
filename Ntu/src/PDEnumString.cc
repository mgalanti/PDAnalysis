#include "PDAnalysis/Ntu/interface/PDEnumString.h"

#include <iostream>
#include <algorithm>

//using namespace std;

std::string PDEnumString::defaultString;
std::map<int,std::string> PDEnumString::    trigPathMapIS;
std::map<std::string,int> PDEnumString::    trigPathMapSI;
std::map<int,std::string> PDEnumString::  trigObjectMapIS;
std::map<std::string,int> PDEnumString::  trigObjectMapSI;
std::map<int,std::string> PDEnumString::  recoObjectMapIS;
std::map<std::string,int> PDEnumString::  recoObjectMapSI;
std::map<int,std::string> PDEnumString::  recoUITypeMapIS;
std::map<std::string,int> PDEnumString::  recoUITypeMapSI;
std::map<int,std::string> PDEnumString::  vertexTypeMapIS;
std::map<std::string,int> PDEnumString::  vertexTypeMapSI;
std::map<int,std::string> PDEnumString::   trackTypeMapIS;
std::map<std::string,int> PDEnumString::   trackTypeMapSI;
std::map<int,std::string> PDEnumString::    muonTypeMapIS;
std::map<std::string,int> PDEnumString::    muonTypeMapSI;
std::map<int,std::string> PDEnumString::     tagTypeMapIS;
std::map<std::string,int> PDEnumString::     tagTypeMapSI;
std::map<int,std::string> PDEnumString::   eleIDTypeMapIS;
std::map<std::string,int> PDEnumString::   eleIDTypeMapSI;
std::map<int,std::string> PDEnumString::  eleMVATypeMapIS;
std::map<std::string,int> PDEnumString::  eleMVATypeMapSI;
std::map<int,std::string> PDEnumString::eleUserFloatMapIS;
std::map<std::string,int> PDEnumString::eleUserFloatMapSI;
std::map<int,std::string> PDEnumString::  eleUserIntMapIS;
std::map<std::string,int> PDEnumString::  eleUserIntMapSI;
std::map<int,std::string> PDEnumString:: elePropertyMapIS;
std::map<std::string,int> PDEnumString:: elePropertyMapSI;

static PDEnumString pdEnumString;

#define mP( NAME ) std::make_pair( #NAME, NAME )

PDEnumString::PDEnumString() {

  defaultString = "";

  trigPathMapSI.clear();
//  trigPathMapSI.insert(
//      mP( dummy ) );

  trigPathMapSI.insert(
      mP( HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi3p5_Muon2_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi_L1_NoOS_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi_Muon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi_NoVertexing_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_LowMass_L1_0er1p5R_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_LowMass_L1_0er1p5_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_LowMass_L1_4R_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_LowMass_L1_4_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_LowMass_L1_TM530_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_LowMass_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Phi_Barrel_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_L1_4p5NoOS_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_L1_4p5_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_L1_4p5er2p0_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_L1_5M_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_L1_5_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_Muon_L1_TM0_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_Muon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0_Upsilon_NoVertexing_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon0er16_Jpsi_NoVertexing_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon10_Jpsi_Barrel_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon10_PsiPrime_Barrel_Seagulls_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon10_Upsilon_Barrel_Seagulls_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon12_Upsilon_eta1p5_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon13_PsiPrime_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon13_Upsilon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon14_Phi_Barrel_Seagulls_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon16_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon18_PsiPrime_noCorrL1_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon18_PsiPrime_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon20_Jpsi_Barrel_Seagulls_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon20_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon24_Phi_noCorrL1_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon24_Upsilon_noCorrL1_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon25_Jpsi_noCorrL1_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon25_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon6_Jpsi_NoVertexing_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon8_PsiPrime_Barrel_v ) );
  trigPathMapSI.insert(
      mP( HLT_Dimuon8_Upsilon_Barrel_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle33_CaloIdL_MW_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle33_CaloIdL_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu18NoFiltersNoVtx_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu20_7_Mass0to30_L1_DM4_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu20_7_Mass0to30_Photon23_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu23NoFiltersNoVtxDisplaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu28NoFiltersNoVtxDisplaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu33NoFiltersNoVtx_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu38NoFiltersNoVtx_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu3_TkMu_DsTau3Mu_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu3_Trk_Tau3mu_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_3_Bs_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_3_Jpsi_Displaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_JpsiTrkTrk_Displaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_JpsiTrk_Displaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_Jpsi_Displaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_Jpsi_NoVertexing_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu4_PsiPrimeTrk_Displaced_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu8_Mass8_PFHT250_v ) );
  trigPathMapSI.insert(
      mP( HLT_DoubleMu8_Mass8_PFHT300_v ) );
  trigPathMapSI.insert(
      mP( HLT_ECALHT800_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele17_CaloIdL_GsfTrkIdVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele20_WPLoose_Gsf_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele27_HighEta_Ele20_Mass55_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele30WP60_Ele8_Mass55_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele30WP60_SC4_Mass55_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v ) );
  trigPathMapSI.insert(
      mP( HLT_L2DoubleMu23_NoVertex_v ) );
  trigPathMapSI.insert(
      mP( HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v ) );
  trigPathMapSI.insert(
      mP( HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu10_CentralPFJet30_BTagCSV_p13_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu16_TkMu0_dEta18_Onia_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu16_TkMu0_dEta18_Phi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_Mu8_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_Mu8_SameSign_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_Mu8_SameSign_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_Mu8_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_TkMu8_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_TrkIsoVVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu17_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu20_Mu10_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu20_Mu10_SameSign_DZ_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu20_Mu10_SameSign_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu20_Mu10_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu20_TkMu0_Phi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu25_TkMu0_Onia_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu25_TkMu0_Phi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu25_TkMu0_dEta18_Onia_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu27_TkMu8_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu30_TkMu0_Onia_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu30_TkMu11_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu3_PFJet40_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu40_TkMu11_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_L2Mu2_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_L2Mu2_Upsilon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_Track2_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_Track2_Upsilon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_Track3p5_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_Track3p5_Upsilon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_Track7_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu7p5_Track7_Upsilon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu8_TrkIsoVVL_v ) );
  trigPathMapSI.insert(
      mP( HLT_Mu8_v ) );
  trigPathMapSI.insert(
      mP( HLT_QuadMuon0_Dimuon0_Jpsi_v ) );
  trigPathMapSI.insert(
      mP( HLT_QuadMuon0_Dimuon0_Upsilon_v ) );
  trigPathMapSI.insert(
      mP( HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v ) );
  trigPathMapSI.insert(
      mP( HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_v ) );
  trigPathMapSI.insert(
      mP( HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_v ) );
  trigPathMapSI.insert(
      mP( HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_v ) );
  trigPathMapSI.insert(
      mP( HLT_Trimuon5_3p5_2_Upsilon_Muon_v ) );
  trigPathMapSI.insert(
      mP( HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v ) );
  trigPathMapSI.insert(
      mP( HLT_TripleMu_12_10_5_v ) );
  trigPathMapSI.insert(
      mP( HLT_TripleMu_5_3_3_v ) );
  trigPathMapSI.insert(
      mP( HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_v  ) );
  trigPathMapSI.insert(
      mP( MC_Ele5_WPTight_Gsf_v  ) );
  revertMap( trigPathMapSI,
             trigPathMapIS );

  trigObjectMapSI.clear();
  trigObjectMapSI.insert( mP( hltJet      ) );
  trigObjectMapSI.insert( mP( hltMuon     ) );
  trigObjectMapSI.insert( mP( hltElectron ) );
  trigObjectMapSI.insert( mP( hltTau      ) );
  trigObjectMapSI.insert( mP( hltTrack    ) );
  revertMap( trigObjectMapSI,
             trigObjectMapIS );

  recoObjectMapSI.clear();
  recoObjectMapSI.insert( mP( recJet      ) );
  recoObjectMapSI.insert( mP( recMuon     ) );
  recoObjectMapSI.insert( mP( recElectron ) );
  recoObjectMapSI.insert( mP( recTau      ) );
  revertMap( recoObjectMapSI,
             recoObjectMapIS );

  recoUITypeMapSI.clear();
  recoUITypeMapSI.insert( mP( puBeta                                                     ) );  // 2018
  recoUITypeMapSI.insert( mP( puBetaStar                                                 ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_loose                     ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_medium                    ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_tight                     ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_veto                      ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_loose                     ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_medium                    ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_tight                     ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_veto                      ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_loose                   ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_medium                  ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_tight                   ) );  // 2018
  recoUITypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_veto                    ) );  // 2018
  recoUITypeMapSI.insert( mP( heepElectronID_HEEPV70                                     ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV1RawValues               ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV1Values                  ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV2RawValues               ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV2Values                  ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV1RawValues             ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV1Values                ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV2RawValues             ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV2Values                ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16GeneralPurposeV1RawValues  ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values     ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16HZZV1RawValues             ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16HZZV1Values                ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalEnergyErrPostCorr                                      ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalEnergyErrPreCorr                                       ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalEnergyPostCorr                                         ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalEnergyPreCorr                                          ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalTrkEnergyErrPostCorr                                   ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalTrkEnergyErrPreCorr                                    ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalTrkEnergyPostCorr                                      ) );  // 2018
  recoUITypeMapSI.insert( mP( ecalTrkEnergyPreCorr                                       ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleDown                                            ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleGainDown                                        ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleGainUp                                          ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleStatDown                                        ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleStatUp                                          ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleSystDown                                        ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleSystUp                                          ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleUp                                              ) );  // 2018
  recoUITypeMapSI.insert( mP( energyScaleValue                                           ) );  // 2018
  recoUITypeMapSI.insert( mP( energySigmaDown                                            ) );  // 2018
  recoUITypeMapSI.insert( mP( energySigmaPhiDown                                         ) );  // 2018
  recoUITypeMapSI.insert( mP( energySigmaPhiUp                                           ) );  // 2018
  recoUITypeMapSI.insert( mP( energySigmaRhoDown                                         ) );  // 2018
  recoUITypeMapSI.insert( mP( energySigmaRhoUp                                           ) );  // 2018
  recoUITypeMapSI.insert( mP( energySigmaUp                                              ) );  // 2018
  recoUITypeMapSI.insert( mP( energySigmaValue                                           ) );  // 2018
  recoUITypeMapSI.insert( mP( energySmearNrSigma                                         ) );  // 2018
  recoUITypeMapSI.insert( mP( heepTrkPtIso                                               ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV1Categories              ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV2Categories              ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV1Categories            ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV2Categories            ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories ) );  // 2018
  recoUITypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16HZZV1Categories            ) );  // 2018
  replaceCharInMap( recoUITypeMapSI, '_', '-'); // Needed to obtain the correct strings
  revertMap( recoUITypeMapSI,
             recoUITypeMapIS );

  vertexTypeMapSI.clear();
  vertexTypeMapSI.insert( mP( svtTagInfo   ) );
  vertexTypeMapSI.insert( mP( svtFitPair   ) );
  vertexTypeMapSI.insert( mP( svtJPsi      ) );
  vertexTypeMapSI.insert( mP( svtKx0       ) );
  vertexTypeMapSI.insert( mP( svtPhi       ) );
  vertexTypeMapSI.insert( mP( svtBuJPsiK   ) );
  vertexTypeMapSI.insert( mP( svtBdJPsiKx  ) );
  vertexTypeMapSI.insert( mP( svtBsJPsiPhi ) );
  vertexTypeMapSI.insert( mP( svtMuMu      ) );
  vertexTypeMapSI.insert( mP( svtBdKxMuMu  ) );
  vertexTypeMapSI.insert( mP( svtK0short   ) );
  vertexTypeMapSI.insert( mP( svtLambda0   ) );
  revertMap( vertexTypeMapSI,
             vertexTypeMapIS );

  trackTypeMapSI.clear();
  trackTypeMapSI.insert( mP( general      ) );
  trackTypeMapSI.insert( mP( packed       ) );
  trackTypeMapSI.insert( mP( pflink       ) );
  trackTypeMapSI.insert( mP( rclink       ) );
  trackTypeMapSI.insert( mP( gplink       ) );
  trackTypeMapSI.insert( mP( muInner      ) );
  trackTypeMapSI.insert( mP( muStandalone ) );
  trackTypeMapSI.insert( mP( muGlobal     ) );
  trackTypeMapSI.insert( mP( muBest       ) );
  trackTypeMapSI.insert( mP( muPF         ) );
  trackTypeMapSI.insert( mP( muReco       ) );
  revertMap( trackTypeMapSI,
             trackTypeMapIS );

  muonTypeMapSI.clear();
  muonTypeMapSI.insert( mP( tmOneStation ) );
  muonTypeMapSI.insert( mP( gmPrompt     ) );
  muonTypeMapSI.insert( mP( pf           ) );
  muonTypeMapSI.insert( mP( loose        ) );
  muonTypeMapSI.insert( mP( medium       ) );
  muonTypeMapSI.insert( mP( tight        ) );
  muonTypeMapSI.insert( mP( soft         ) );
  muonTypeMapSI.insert( mP( highPt       ) );
  muonTypeMapSI.insert( mP( inner        ) );
  muonTypeMapSI.insert( mP( standalone   ) );
  muonTypeMapSI.insert( mP( global       ) );
  revertMap( muonTypeMapSI,
             muonTypeMapIS );

  tagTypeMapSI.insert( mP( pfDeepCSVJetTags_probudsg ) );
  tagTypeMapSI.insert( mP( pfDeepCSVJetTags_probc    ) );
  tagTypeMapSI.insert( mP( pfDeepCSVJetTags_probcc   ) );
  tagTypeMapSI.insert( mP( pfDeepCSVJetTags_probb    ) );
  tagTypeMapSI.insert( mP( pfDeepCSVJetTags_probbb   ) );
  revertMap( tagTypeMapSI,
             tagTypeMapIS );
  
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_loose    ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_medium   ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_tight    ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_veto     ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_loose    ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_medium   ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_tight    ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_veto     ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_loose  ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_medium ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_tight  ) );  // 2018
  eleIDTypeMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_veto   ) );  // 2018
  eleIDTypeMapSI.insert( mP( heepElectronID_HEEPV70                    ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_iso_V1_wp80               ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_iso_V1_wp90               ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_iso_V1_wpLoose            ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_iso_V2_wp80               ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_iso_V2_wp90               ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_iso_V2_wpHZZ              ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_iso_V2_wpLoose            ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_noIso_V1_wp80             ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_noIso_V1_wp90             ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_noIso_V1_wpLoose          ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_noIso_V2_wp80             ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_noIso_V2_wp90             ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Fall17_noIso_V2_wpLoose          ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Spring16_GeneralPurpose_V1_wp80  ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Spring16_GeneralPurpose_V1_wp90  ) );  // 2018
  eleIDTypeMapSI.insert( mP( mvaEleID_Spring16_HZZ_V1_wpLoose          ) );  // 2018
  replaceCharInMap( eleIDTypeMapSI, '_', '-'); // Needed to obtain the correct strings
  revertMap( eleIDTypeMapSI,
              eleIDTypeMapIS );
    
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV1RawValues              ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV1Values                 ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV2RawValues              ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV2Values                 ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV1RawValues            ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV1Values               ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV2RawValues            ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV2Values               ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16GeneralPurposeV1RawValues ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values    ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16HZZV1RawValues            ) );
  eleMVATypeMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16HZZV1Values               ) );
  revertMap( eleMVATypeMapSI,
             eleMVATypeMapIS );

  eleUserFloatMapSI.insert( mP( ecalEnergyErrPostCorr    ) );
  eleUserFloatMapSI.insert( mP( ecalEnergyErrPreCorr     ) );
  eleUserFloatMapSI.insert( mP( ecalEnergyPostCorr       ) );
  eleUserFloatMapSI.insert( mP( ecalEnergyPreCorr        ) );
  eleUserFloatMapSI.insert( mP( ecalTrkEnergyErrPostCorr ) );
  eleUserFloatMapSI.insert( mP( ecalTrkEnergyErrPreCorr  ) );
  eleUserFloatMapSI.insert( mP( ecalTrkEnergyPostCorr    ) );
  eleUserFloatMapSI.insert( mP( ecalTrkEnergyPreCorr     ) );
  eleUserFloatMapSI.insert( mP( energyScaleDown          ) );
  eleUserFloatMapSI.insert( mP( energyScaleGainDown      ) );
  eleUserFloatMapSI.insert( mP( energyScaleGainUp        ) );
  eleUserFloatMapSI.insert( mP( energyScaleStatDown      ) );
  eleUserFloatMapSI.insert( mP( energyScaleStatUp        ) );
  eleUserFloatMapSI.insert( mP( energyScaleSystDown      ) );
  eleUserFloatMapSI.insert( mP( energyScaleSystUp        ) );
  eleUserFloatMapSI.insert( mP( energyScaleUp            ) );
  eleUserFloatMapSI.insert( mP( energyScaleValue         ) );
  eleUserFloatMapSI.insert( mP( energySigmaDown          ) );
  eleUserFloatMapSI.insert( mP( energySigmaPhiDown       ) );
  eleUserFloatMapSI.insert( mP( energySigmaPhiUp         ) );
  eleUserFloatMapSI.insert( mP( energySigmaRhoDown       ) );
  eleUserFloatMapSI.insert( mP( energySigmaRhoUp         ) );
  eleUserFloatMapSI.insert( mP( energySigmaUp            ) );
  eleUserFloatMapSI.insert( mP( energySigmaValue         ) );
  eleUserFloatMapSI.insert( mP( energySmearNrSigma       ) );
  eleUserFloatMapSI.insert( mP( heepTrkPtIso             ) );
  revertMap( eleUserFloatMapSI,
             eleUserFloatMapIS );

  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_loose                     ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_medium                    ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_tight                     ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V1_veto                      ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_loose                     ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_medium                    ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_tight                     ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Fall17_94X_V2_veto                      ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_loose                   ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_medium                  ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_tight                   ) );  // 2018
  eleUserIntMapSI.insert( mP( cutBasedElectronID_Summer16_80X_V1_veto                    ) );  // 2018
  eleUserIntMapSI.insert( mP( heepElectronID_HEEPV70                                     ) );  // 2018
  eleUserIntMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV1Categories              ) );  // 2018
  eleUserIntMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17IsoV2Categories              ) );  // 2018
  eleUserIntMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV1Categories            ) );  // 2018
  eleUserIntMapSI.insert( mP( ElectronMVAEstimatorRun2Fall17NoIsoV2Categories            ) );  // 2018
  eleUserIntMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories ) );  // 2018
  eleUserIntMapSI.insert( mP( ElectronMVAEstimatorRun2Spring16HZZV1Categories            ) );  // 2018
  replaceCharInMap( eleUserIntMapSI, '_', '-'); // Needed to obtain the correct strings
  revertMap( eleUserIntMapSI,
             eleUserIntMapIS );

  elePropertyMapSI.insert( mP( isPF                          ) );
  elePropertyMapSI.insert( mP( passConversionVeto            ) );
  elePropertyMapSI.insert( mP( ambiguous                     ) );
  elePropertyMapSI.insert( mP( isEB                          ) );
  elePropertyMapSI.insert( mP( isEBEEGap                     ) );
  elePropertyMapSI.insert( mP( isEBEtaGap                    ) );
  elePropertyMapSI.insert( mP( isEBGap                       ) );
  elePropertyMapSI.insert( mP( isEBPhiGap                    ) );
  elePropertyMapSI.insert( mP( isEcalEnergyCorrected         ) );
  elePropertyMapSI.insert( mP( isEE                          ) );
  elePropertyMapSI.insert( mP( isEEDeeGap                    ) );
  elePropertyMapSI.insert( mP( isEEGap                       ) );
  elePropertyMapSI.insert( mP( isEERingGap                   ) );
  elePropertyMapSI.insert( mP( isElectron                    ) );
  elePropertyMapSI.insert( mP( isEnergyScaleCorrected        ) );
  elePropertyMapSI.insert( mP( isGap                         ) );
  elePropertyMapSI.insert( mP( isGsfCtfChargeConsistent      ) );
  elePropertyMapSI.insert( mP( isGsfCtfScPixChargeConsistent ) );
  elePropertyMapSI.insert( mP( isGsfScPixChargeConsistent    ) );
  elePropertyMapSI.insert( mP( isSeedSaturated               ) );
  revertMap( elePropertyMapSI,
             elePropertyMapIS );
}


PDEnumString::~PDEnumString() {
}


void
PDEnumString::resetTrigMap() {
  trigPathMapIS.clear();
  trigPathMapSI.clear();
  return;
}


void
PDEnumString::revertTrigMap() {
  revertMap( trigPathMapSI,
             trigPathMapIS );
  return;
}


const std::map<int,std::string>&
PDEnumString::trigMap() {
  return trigPathMapIS;
}


const std::map<int,std::string>&
PDEnumString::eleIDTypeMap() {
  return eleIDTypeMapIS;
}


const std::map<int,std::string>&
PDEnumString::eleMVATypeMap() {
  return eleMVATypeMapIS;
}


const std::map<int,std::string>&
PDEnumString::eleUserFloatMap() {
  return eleUserFloatMapIS;
}


const std::map<int,std::string>&
PDEnumString::eleUserIntMap() {
  return eleUserIntMapIS;
}


const std::map<int,std::string>&
PDEnumString::elePropertyMap() {
  return elePropertyMapIS;
}


std::string PDEnumString::trigBase( const std::string& trigPath ) {
  int pathLength = trigPath.length();
  const char* str = trigPath.c_str() + pathLength - 1;
  while ( ( *str >= '0' ) && ( *str <= '9' ) ) {
    --str;
    --pathLength;
  }
  return trigPath.substr( 0, pathLength );
}


int
PDEnumString::findTrigPath( const std::string& trigPath,
                            bool fixedList, int code ) {
//  std::string trigName = trigBase( trigPath );
//  std::map<std::string,int>::const_iterator
//    iter = trigPathMapSI.find( trigName );
  int path = find( trigPath, trigPathMapSI );
//  std::map<std::string,int>::const_iterator
//    iter = trigPathMapSI.find( trigPath );
//  std::map<std::string,int>::const_iterator
//    iend = trigPathMapSI.end();
//  if ( iter != iend ) return iter->second;
  if ( path >= 0 ) return path;
  if ( fixedList ) return -1;
  int nextSlot = ( code < 0 ? trigPathMapIS.rbegin()->first + 1 : code );
//  std::cout << "add " << nextSlot << " " << trigPath << std::endl;
//  trigPathMapSI.insert( std::make_pair( trigName, nextSlot ) );
//  trigPathMapIS.insert( std::make_pair( nextSlot, trigName ) );
  trigPathMapSI.insert( std::make_pair( trigPath, nextSlot ) );
  trigPathMapIS.insert( std::make_pair( nextSlot, trigPath ) );
  return nextSlot;
}


const std::string&
PDEnumString::findTrigPath( int trigPath ) {
  return find( trigPath, trigPathMapIS );
//  std::map<int,std::string>::const_iterator
//    iter = trigPathMapIS.find( trigPath );
//  std::map<int,std::string>::const_iterator
//    iend = trigPathMapIS.end();
//  if ( iter != iend ) return iter->second;
//  return defaultString;
}


int
PDEnumString::findTrigObject( const std::string& trigObject ) {
  return find( trigObject, trigObjectMapSI );
//  std::map<std::string,int>::const_iterator
//    iter = trigObjectMapSI.find( trigObject );
//  std::map<std::string,int>::const_iterator
//    iend = trigObjectMapSI.end();
//  if ( iter != iend ) return iter->second;
//  return -1;
}


const std::string&
PDEnumString::findTrigObject( int trigObject ) {
  return find( trigObject, trigObjectMapIS );
//  std::map<int,std::string>::const_iterator
//    iter = trigObjectMapIS.find( trigObject );
//  std::map<int,std::string>::const_iterator
//    iend = trigObjectMapIS.end();
//  if ( iter != iend ) return iter->second;
//  return defaultString;
}


int
PDEnumString::findRecoObject( const std::string& recoObject ) {
  return find( recoObject, recoObjectMapSI );
//  std::map<std::string,int>::const_iterator
//    iter = recoObjectMapSI.find( recoObject );
//  std::map<std::string,int>::const_iterator
//    iend = recoObjectMapSI.end();
//  if ( iter != iend ) return iter->second;
//  return -1;
}


const std::string&
PDEnumString::findRecoObject( int recoObject ) {
  return find( recoObject, recoObjectMapIS );
//  std::map<int,std::string>::const_iterator
//    iter = recoObjectMapIS.find( recoObject );
//  std::map<int,std::string>::const_iterator
//    iend = recoObjectMapIS.end();
//  if ( iter != iend ) return iter->second;
//  return defaultString;
}


int
PDEnumString::findRecoUIType( const std::string& recoUIType ) {
  return find( recoUIType, recoUITypeMapSI );
//  std::map<std::string,int>::const_iterator
//    iter = recoUITypeMapSI.find( recoUIType );
//  std::map<std::string,int>::const_iterator
//    iend = recoUITypeMapSI.end();
//  if ( iter != iend ) return iter->second;
//  return -1;
}


const std::string&
PDEnumString::findRecoUIType( int recoUIType ) {
  return find( recoUIType, recoUITypeMapIS );
//  std::map<int,std::string>::const_iterator
//    iter = recoUITypeMapIS.find( recoUIType );
//  std::map<int,std::string>::const_iterator
//    iend = recoUITypeMapIS.end();
//  if ( iter != iend ) return iter->second;
//  return defaultString;
}


int
PDEnumString::findVertexType( const std::string& vtxType ) {
  return find( vtxType, vertexTypeMapSI );
//  std::map<std::string,int>::const_iterator
//    iter = vertexTypeMapSI.find( vtxType );
//  std::map<std::string,int>::const_iterator
//    iend = vertexTypeMapSI.end();
//  if ( iter != iend ) return iter->second;
//  return -1;
}


const std::string&
PDEnumString::findVertexType( int vtxType ) {
  return find( vtxType, vertexTypeMapIS );
//  std::map<int,std::string>::const_iterator
//    iter = vertexTypeMapIS.find( vtxType );
//  std::map<int,std::string>::const_iterator
//    iend = vertexTypeMapIS.end();
//  if ( iter != iend ) return iter->second;
//  return defaultString;
}


int
PDEnumString::findTrackType( const std::string& trkType ) {
  return find( trkType, trackTypeMapSI );
}


const std::string&
PDEnumString::findTrackType( int trkType ) {
  return find( trkType, trackTypeMapIS );
}


int
PDEnumString::findMuonType( const std::string& muoType ) {
  return find( muoType, muonTypeMapSI );
}


const std::string&
PDEnumString::findMuonType( int muoType ) {
  return find( muoType, muonTypeMapIS );
}


int
PDEnumString::findTagType( const std::string& tagType ) {
  return find( tagType, tagTypeMapSI );
}


const std::string&
PDEnumString::findTagType( int tagType ) {
  return find( tagType, tagTypeMapIS );
}


int
PDEnumString::findEleIDType( const std::string& eleIDType ) {
  return find( eleIDType, eleIDTypeMapSI );
}


const std::string&
PDEnumString::findEleIDType( int eleIDType ) {
  return find( eleIDType, eleIDTypeMapIS );
}


int
PDEnumString::findEleMVAType( const std::string& eleMVAType ) {
  return find( eleMVAType, eleMVATypeMapSI );
}


const std::string&
PDEnumString::findEleMVAType( int eleMVAType ) {
  return find( eleMVAType, eleMVATypeMapIS );
}


int
PDEnumString::findEleUserFloat( const std::string& eleUserFloat ) {
  return find( eleUserFloat, eleUserFloatMapSI );
}


const std::string&
PDEnumString::findEleUserFloat( int eleUserFloat ) {
  return find( eleUserFloat, eleUserFloatMapIS );
}


int
PDEnumString::findEleUserInt( const std::string& eleUserInt ) {
  return find( eleUserInt, eleUserIntMapSI );
}


const std::string&
PDEnumString::findEleUserInt( int eleUserInt ) {
  return find( eleUserInt, eleUserIntMapIS );
}


int
PDEnumString::findEleProperty( const std::string& eleProperty ) {
  return find( eleProperty, elePropertyMapSI );
}


const std::string&
PDEnumString::findEleProperty( int eleProperty ) {
  return find( eleProperty, elePropertyMapIS );
}


void PDEnumString::revertMap( const std::map<std::string,int>& mapSI,
                                    std::map<int,std::string>& mapIS ) {
  mapIS.clear();
  std::map<std::string,int>::const_iterator iter = mapSI.begin();
  std::map<std::string,int>::const_iterator iend = mapSI.end();
  while ( iter != iend ) {
    const std::pair<std::string,int>& entry = *iter++;
    mapIS.insert( std::make_pair( entry.second, entry.first ) );
  }
  return;
}


void PDEnumString::replaceCharInMap( std::map<std::string,int>& mapSI, 
                                     const char from, const char to ) {
  std::map<std::string,int> newMapSI;
  std::map<std::string,int>::iterator iter = mapSI.begin();
  std::map<std::string,int>::iterator iend = mapSI.end();
  while ( iter != iend) {
    std::string str = iter->first;
    std::replace( str.begin(), str.end(), from, to );
    newMapSI[str] = iter->second;
    iter++;
  }
  mapSI = newMapSI;
}


int
PDEnumString::find( const std::string& name,
                    const std::map<std::string,int>& eMap ) {
  std::map<std::string,int>::const_iterator iter = eMap.find( name );
  std::map<std::string,int>::const_iterator iend = eMap.end();
  if ( iter != iend ) return iter->second;
  return -1;
}


const std::string&
PDEnumString::find( int code,
                    const std::map<int,std::string>& eMap ) {
  std::map<int,std::string>::const_iterator iter = eMap.find( code );
  std::map<int,std::string>::const_iterator iend = eMap.end();
  if ( iter != iend ) return iter->second;
  return defaultString;
}

