#ifndef PDEnumString_H
#define PDEnumString_H

#include <map>
#include <string>

class PDEnumString {

 public:

  PDEnumString();
  virtual ~PDEnumString();

//  enum trigPath   { dummy = 0 };
  enum trigPath   { HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v = 1,
                    HLT_Dimuon0_Jpsi3p5_Muon2_v,
                    HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v,
                    HLT_Dimuon0_Jpsi_L1_NoOS_v,
                    HLT_Dimuon0_Jpsi_Muon_v,
                    HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v,
                    HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v,
                    HLT_Dimuon0_Jpsi_NoVertexing_v,
                    HLT_Dimuon0_Jpsi_v,
                    HLT_Dimuon0_LowMass_L1_0er1p5R_v,
                    HLT_Dimuon0_LowMass_L1_0er1p5_v,
                    HLT_Dimuon0_LowMass_L1_4R_v,
                    HLT_Dimuon0_LowMass_L1_4_v,
                    HLT_Dimuon0_LowMass_L1_TM530_v,
                    HLT_Dimuon0_LowMass_v,
                    HLT_Dimuon0_Phi_Barrel_v,
                    HLT_Dimuon0_Upsilon_L1_4p5NoOS_v,
                    HLT_Dimuon0_Upsilon_L1_4p5_v,
                    HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v,
                    HLT_Dimuon0_Upsilon_L1_4p5er2p0_v,
                    HLT_Dimuon0_Upsilon_L1_5M_v,
                    HLT_Dimuon0_Upsilon_L1_5_v,
                    HLT_Dimuon0_Upsilon_Muon_L1_TM0_v,
                    HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v,
                    HLT_Dimuon0_Upsilon_Muon_v,
                    HLT_Dimuon0_Upsilon_NoVertexing_v,
                    HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_v,
                    HLT_Dimuon0er16_Jpsi_NoVertexing_v,
                    HLT_Dimuon10_Jpsi_Barrel_v,
                    HLT_Dimuon10_PsiPrime_Barrel_Seagulls_v,
                    HLT_Dimuon10_Upsilon_Barrel_Seagulls_v,
                    HLT_Dimuon12_Upsilon_eta1p5_v,
                    HLT_Dimuon13_PsiPrime_v,
                    HLT_Dimuon13_Upsilon_v,
                    HLT_Dimuon14_Phi_Barrel_Seagulls_v,
                    HLT_Dimuon16_Jpsi_v,
                    HLT_Dimuon18_PsiPrime_noCorrL1_v,
                    HLT_Dimuon18_PsiPrime_v,
                    HLT_Dimuon20_Jpsi_Barrel_Seagulls_v,
                    HLT_Dimuon20_Jpsi_v,
                    HLT_Dimuon24_Phi_noCorrL1_v,
                    HLT_Dimuon24_Upsilon_noCorrL1_v,
                    HLT_Dimuon25_Jpsi_noCorrL1_v,
                    HLT_Dimuon25_Jpsi_v,
                    HLT_Dimuon6_Jpsi_NoVertexing_v,
                    HLT_Dimuon8_PsiPrime_Barrel_v,
                    HLT_Dimuon8_Upsilon_Barrel_v,
                    HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v,
                    HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v,
                    HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v,
                    HLT_DoubleEle33_CaloIdL_MW_v,
                    HLT_DoubleEle33_CaloIdL_v,
                    HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v,
                    HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_v,
                    HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v,
                    HLT_DoubleMu18NoFiltersNoVtx_v,
                    HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v,
                    HLT_DoubleMu20_7_Mass0to30_L1_DM4_v,
                    HLT_DoubleMu20_7_Mass0to30_Photon23_v,
                    HLT_DoubleMu23NoFiltersNoVtxDisplaced_v,
                    HLT_DoubleMu28NoFiltersNoVtxDisplaced_v,
                    HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v,
                    HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi_v,
                    HLT_DoubleMu33NoFiltersNoVtx_v,
                    HLT_DoubleMu38NoFiltersNoVtx_v,
                    HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v,
                    HLT_DoubleMu3_TkMu_DsTau3Mu_v,
                    HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v,
                    HLT_DoubleMu3_Trk_Tau3mu_v,
                    HLT_DoubleMu4_3_Bs_v,
                    HLT_DoubleMu4_3_Jpsi_Displaced_v,
                    HLT_DoubleMu4_JpsiTrkTrk_Displaced_v,
                    HLT_DoubleMu4_JpsiTrk_Displaced_v,
                    HLT_DoubleMu4_Jpsi_Displaced_v,
                    HLT_DoubleMu4_Jpsi_NoVertexing_v,
                    HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v,
                    HLT_DoubleMu4_PsiPrimeTrk_Displaced_v,
                    HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v,
                    HLT_DoubleMu8_Mass8_PFHT250_v,
                    HLT_DoubleMu8_Mass8_PFHT300_v,
                    HLT_ECALHT800_v,
                    HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_v,
                    HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v,
                    HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v,
                    HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v,
                    HLT_Ele17_CaloIdL_GsfTrkIdVL_v,
                    HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v,
                    HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v,
                    HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v,
                    HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Ele20_WPLoose_Gsf_v,
                    HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v,
                    HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v,
                    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v,
                    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v,
                    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Ele27_HighEta_Ele20_Mass55_v,
                    HLT_Ele30WP60_Ele8_Mass55_v,
                    HLT_Ele30WP60_SC4_Mass55_v,
                    HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v,
                    HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v,
                    HLT_L2DoubleMu23_NoVertex_v,
                    HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v,
                    HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_v,
                    HLT_Mu10_CentralPFJet30_BTagCSV_p13_v,
                    HLT_Mu16_TkMu0_dEta18_Onia_v,
                    HLT_Mu16_TkMu0_dEta18_Phi_v,
                    HLT_Mu17_Mu8_DZ_v,
                    HLT_Mu17_Mu8_SameSign_DZ_v,
                    HLT_Mu17_Mu8_SameSign_v,
                    HLT_Mu17_Mu8_v,
                    HLT_Mu17_TkMu8_DZ_v,
                    HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v,
                    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v,
                    HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v,
                    HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v,
                    HLT_Mu17_TrkIsoVVL_v,
                    HLT_Mu17_v,
                    HLT_Mu20_Mu10_DZ_v,
                    HLT_Mu20_Mu10_SameSign_DZ_v,
                    HLT_Mu20_Mu10_SameSign_v,
                    HLT_Mu20_Mu10_v,
                    HLT_Mu20_TkMu0_Phi_v,
                    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Mu25_TkMu0_Onia_v,
                    HLT_Mu25_TkMu0_Phi_v,
                    HLT_Mu25_TkMu0_dEta18_Onia_v,
                    HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL_v,
                    HLT_Mu27_TkMu8_v,
                    HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v,
                    HLT_Mu30_TkMu0_Onia_v,
                    HLT_Mu30_TkMu11_v,
                    HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL_v,
                    HLT_Mu3_PFJet40_v,
                    HLT_Mu40_TkMu11_v,
                    HLT_Mu7p5_L2Mu2_Jpsi_v,
                    HLT_Mu7p5_L2Mu2_Upsilon_v,
                    HLT_Mu7p5_Track2_Jpsi_v,
                    HLT_Mu7p5_Track2_Upsilon_v,
                    HLT_Mu7p5_Track3p5_Jpsi_v,
                    HLT_Mu7p5_Track3p5_Upsilon_v,
                    HLT_Mu7p5_Track7_Jpsi_v,
                    HLT_Mu7p5_Track7_Upsilon_v,
                    HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v,
                    HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_v,
                    HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v,
                    HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v,
                    HLT_Mu8_TrkIsoVVL_v,
                    HLT_Mu8_v,
                    HLT_QuadMuon0_Dimuon0_Jpsi_v,
                    HLT_QuadMuon0_Dimuon0_Upsilon_v,
                    HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v,
                    HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_v,
                    HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_v,
                    HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_v,
                    HLT_Trimuon5_3p5_2_Upsilon_Muon_v,
                    HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v,
                    HLT_TripleMu_12_10_5_v,
                    HLT_TripleMu_5_3_3_v,
                    HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_v,
                    MC_Ele5_WPTight_Gsf_v,
                    DUMMY_LAST_TRIGGER };
  enum trigObject { hltJet = 1,
                    hltMuon,
                    hltElectron,
                    hltTau,
                    hltTrack };
  enum recoObject { recJet = 1,
                    recMuon,
                    recElectron,
                    recTau,
                    recSvt,
                    recPV };
  // MG: this now starts from 10 because otherwise its values overlap
  //     with the ones from the eleIDType map. Values from both maps are
  //     put into the recoUIType map so we need unique int identifiers.
  enum infoType   { puBeta = 10,
                    puBetaStar };
  enum vertexType { svtTagInfo = 1,
                    svtFitPair,
                    svtJPsi,
                    svtKx0,
                    svtPhi,
                    svtBuJPsiK,
                    svtBdJPsiKx,
                    svtBsJPsiPhi,
                    svtMuMu,
                    svtBdKxMuMu,
                    svtK0short,
                    svtLambda0,
                    DUMMY_LAST_DECAY };
  enum trackType  { general      =     1,
                    packed       =     2,
                    pflink       =     4,
                    rclink       =     8,
                    gplink       =    16,
                    muInner      =  1024,
                    muStandalone =  2048,
                    muGlobal     =  4096,
                    muBest       =  8192,
                    muPF         = 16384,
                    muReco       = 32768 };
  enum muonType   { tmOneStation =    1,
                    gmPrompt     =    2,
                    pf           = pflink,
                    loose        =    8,
                    medium       =   16,
                    tight        =   32,
                    soft         =   64,
                    highPt       =  128,
                    inner        = muInner,
                    standalone   = muStandalone,
                    global       = muGlobal };
  enum tagType    { pfDeepCSVJetTags_probudsg = 11,
                    pfDeepCSVJetTags_probc    = 21,
                    pfDeepCSVJetTags_probcc   = 22,
                    pfDeepCSVJetTags_probb    = 31,
                    pfDeepCSVJetTags_probbb   = 32,
                    DUMMY_LAST_TAG };

  // MG: eleIDs found in 2018 MC
  enum eleIDType  { cutBasedElectronID_Fall17_94X_V1_loose    =         1,    // 2018
                    cutBasedElectronID_Fall17_94X_V1_medium   =         2,    // 2018
                    cutBasedElectronID_Fall17_94X_V1_tight    =         4,    // 2018
                    cutBasedElectronID_Fall17_94X_V1_veto     =         8,    // 2018
                    cutBasedElectronID_Fall17_94X_V2_loose    =        16,    // 2018
                    cutBasedElectronID_Fall17_94X_V2_medium   =        32,    // 2018
                    cutBasedElectronID_Fall17_94X_V2_tight    =        64,    // 2018
                    cutBasedElectronID_Fall17_94X_V2_veto     =       128,    // 2018
                    cutBasedElectronID_Summer16_80X_V1_loose  =       256,    // 2018
                    cutBasedElectronID_Summer16_80X_V1_medium =       512,    // 2018
                    cutBasedElectronID_Summer16_80X_V1_tight  =      1024,    // 2018
                    cutBasedElectronID_Summer16_80X_V1_veto   =      2048,    // 2018
                    heepElectronID_HEEPV70                    =      4096,    // 2018
                    mvaEleID_Fall17_iso_V1_wp80               =      8192,    // 2018
                    mvaEleID_Fall17_iso_V1_wp90               =     16384,    // 2018
                    mvaEleID_Fall17_iso_V1_wpLoose            =     32768,    // 2018
                    mvaEleID_Fall17_iso_V2_wp80               =     65536,    // 2018
                    mvaEleID_Fall17_iso_V2_wp90               =    131072,    // 2018
                    mvaEleID_Fall17_iso_V2_wpHZZ              =    262144,    // 2018
                    mvaEleID_Fall17_iso_V2_wpLoose            =    524288,    // 2018
                    mvaEleID_Fall17_noIso_V1_wp80             =   1048576,    // 2018
                    mvaEleID_Fall17_noIso_V1_wp90             =   2097152,    // 2018
                    mvaEleID_Fall17_noIso_V1_wpLoose          =   4194304,    // 2018
                    mvaEleID_Fall17_noIso_V2_wp80             =   8388608,    // 2018
                    mvaEleID_Fall17_noIso_V2_wp90             =  16777216,    // 2018
                    mvaEleID_Fall17_noIso_V2_wpLoose          =  33554432,    // 2018
                    mvaEleID_Spring16_GeneralPurpose_V1_wp80  =  67108864,    // 2018
                    mvaEleID_Spring16_GeneralPurpose_V1_wp90  = 134217728,    // 2018
                    mvaEleID_Spring16_HZZ_V1_wpLoose          = 268435456 };  // 2018

  enum eleMVAType { ElectronMVAEstimatorRun2Fall17IsoV1RawValues = 10001,
                    ElectronMVAEstimatorRun2Fall17IsoV1Values,
                    ElectronMVAEstimatorRun2Fall17IsoV2RawValues,
                    ElectronMVAEstimatorRun2Fall17IsoV2Values,
                    ElectronMVAEstimatorRun2Fall17NoIsoV1RawValues,
                    ElectronMVAEstimatorRun2Fall17NoIsoV1Values,
                    ElectronMVAEstimatorRun2Fall17NoIsoV2RawValues,
                    ElectronMVAEstimatorRun2Fall17NoIsoV2Values,
                    ElectronMVAEstimatorRun2Spring16GeneralPurposeV1RawValues,
                    ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values,
                    ElectronMVAEstimatorRun2Spring16HZZV1RawValues,
                    ElectronMVAEstimatorRun2Spring16HZZV1Values };
  
  enum eleUserFloat { ecalEnergyErrPostCorr = 20001,
                      ecalEnergyErrPreCorr,
                      ecalEnergyPostCorr,
                      ecalEnergyPreCorr,
                      ecalTrkEnergyErrPostCorr,
                      ecalTrkEnergyErrPreCorr,
                      ecalTrkEnergyPostCorr,
                      ecalTrkEnergyPreCorr,
                      energyScaleDown,
                      energyScaleGainDown,
                      energyScaleGainUp,
                      energyScaleStatDown,
                      energyScaleStatUp,
                      energyScaleSystDown,
                      energyScaleSystUp,
                      energyScaleUp,
                      energyScaleValue,
                      energySigmaDown,
                      energySigmaPhiDown,
                      energySigmaPhiUp,
                      energySigmaRhoDown,
                      energySigmaRhoUp,
                      energySigmaUp,
                      energySigmaValue,
                      energySmearNrSigma,
                      heepTrkPtIso };
                      
  enum eleUserInt { ElectronMVAEstimatorRun2Fall17IsoV1Categories = 30001,
                    ElectronMVAEstimatorRun2Fall17IsoV2Categories,
                    ElectronMVAEstimatorRun2Fall17NoIsoV1Categories,
                    ElectronMVAEstimatorRun2Fall17NoIsoV2Categories,
                    ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories,
                    ElectronMVAEstimatorRun2Spring16HZZV1Categories };
  
  enum eleProperty { isPF                               = 1,
                     passConversionVeto                 = 2,
                     ambiguous                          = 4,
                     isEB                               = 8,
                     isEBEEGap                         = 16,
                     isEBEtaGap                        = 32,
                     isEBGap                           = 64,
                     isEBPhiGap                       = 128,
                     isEcalEnergyCorrected            = 256,
                     isEE                             = 512,
                     isEEDeeGap                      = 1024,
                     isEEGap                         = 2048,
                     isEERingGap                     = 4096,
                     isElectron                      = 8192,
                     isEnergyScaleCorrected         = 16384,
                     isGap                          = 32768,
                     isGsfCtfChargeConsistent       = 65536,
                     isGsfCtfScPixChargeConsistent = 131072,
                     isGsfScPixChargeConsistent    = 262144,
                     isSeedSaturated               = 524288 };
  
  static void  resetTrigMap();
  static void revertTrigMap();
  static const std::map<int,std::string>& trigMap();
  static const std::map<int,std::string>& eleIDTypeMap();  // For Paolo: move these declarations below if you think it is more appropriate
  static const std::map<int,std::string>& eleMVATypeMap();
  static const std::map<int,std::string>& eleUserFloatMap();
  static const std::map<int,std::string>& eleUserIntMap();
  static const std::map<int,std::string>& elePropertyMap();
  static std::string trigBase( const std::string& trigPath );

  static int                findTrigPath  ( const std::string& trigPath,
                                            bool fixedList = true,
                                            int code = -1 );
  static const std::string& findTrigPath    ( int                trigPath     );
  static int                findTrigObject  ( const std::string& trigObject   );
  static const std::string& findTrigObject  ( int                trigObject   );
  static int                findRecoObject  ( const std::string& recoObject   );
  static const std::string& findRecoObject  ( int                recoObject   );
  static int                findRecoUIType  ( const std::string& recoUIType   );
  static const std::string& findRecoUIType  ( int                recoUIType   );
  static int                findVertexType  ( const std::string& vtxType      );
  static const std::string& findVertexType  ( int                vtxType      );
  static int                findTrackType   ( const std::string& trkType      );
  static const std::string& findTrackType   ( int                trkType      );
  static int                findMuonType    ( const std::string& muoType      );
  static const std::string& findMuonType    ( int                muoType      );
  static int                findTagType     ( const std::string& tagType      );
  static const std::string& findTagType     ( int                tagType      );
  static int                findEleIDType   ( const std::string& eleIDType    );
  static const std::string& findEleIDType   ( int                eleIDType    );
  static int                findEleMVAType  ( const std::string& eleMVAType   );
  static const std::string& findEleMVAType  ( int                eleMVAType   );
  static int                findEleUserFloat( const std::string& eleUserFloat );
  static const std::string& findEleUserFloat( int                eleUserFloat );  
  static int                findEleUserInt  ( const std::string& eleUserInt   );
  static const std::string& findEleUserInt  ( int                eleUserInt   );
  static int                findEleProperty ( const std::string& eleProperty  );
  static const std::string& findEleProperty ( int                eleProperty  );
  
 private:

  static std::string defaultString;
  static std::map<int,std::string>     trigPathMapIS;
  static std::map<std::string,int>     trigPathMapSI;
  static std::map<int,std::string>   trigObjectMapIS;
  static std::map<std::string,int>   trigObjectMapSI;
  static std::map<int,std::string>   recoObjectMapIS;
  static std::map<std::string,int>   recoObjectMapSI;
  static std::map<int,std::string>   recoUITypeMapIS;
  static std::map<std::string,int>   recoUITypeMapSI;
  static std::map<int,std::string>   vertexTypeMapIS;
  static std::map<std::string,int>   vertexTypeMapSI;
  static std::map<int,std::string>    trackTypeMapIS;
  static std::map<std::string,int>    trackTypeMapSI;
  static std::map<int,std::string>     muonTypeMapIS;
  static std::map<std::string,int>     muonTypeMapSI;
  static std::map<int,std::string>      tagTypeMapIS;
  static std::map<std::string,int>      tagTypeMapSI;
  static std::map<int,std::string>    eleIDTypeMapIS;
  static std::map<std::string,int>    eleIDTypeMapSI;
  static std::map<int,std::string>   eleMVATypeMapIS;
  static std::map<std::string,int>   eleMVATypeMapSI;
  static std::map<int,std::string> eleUserFloatMapIS;
  static std::map<std::string,int> eleUserFloatMapSI;
  static std::map<int,std::string>   eleUserIntMapIS;
  static std::map<std::string,int>   eleUserIntMapSI;
  static std::map<int,std::string>  elePropertyMapIS;
  static std::map<std::string,int>  elePropertyMapSI;

  static void revertMap( const std::map<std::string,int>& mapSI,
                               std::map<int,std::string>& mapIS );
  
  // Added by Mario. Needed to obtain the correct 
  // strings for various electron-related type names
  static void replaceCharInMap( std::map<std::string,int>& mapSI, 
                                const char from, const char to );

  static int                find( const std::string& name,
                                  const std::map<std::string,int>& eMap );
  static const std::string& find( int                code,
                                  const std::map<int,std::string>& eMap );

};


#endif

