#define UTIL_USE FULL

#include "TSystem.h"

#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"

#define insertI( MAP, VAR )  MAP.insert(std::make_pair<std::string, int*>(#VAR,&VAR))

#define insertVF( MAP, VAR ) MAP.insert(std::make_pair<std::string, std::vector<float>**>(#VAR,&VAR))
#define insertVI( MAP, VAR ) MAP.insert(std::make_pair<std::string, std::vector<int>**>(#VAR,&VAR))
#define insertVB( MAP, VAR ) MAP.insert(std::make_pair<std::string, std::vector<bool>**>(#VAR,&VAR))



MGBaseAnalyzer::MGBaseAnalyzer()
{
  std::cout << "MGBaseAnalyzer::MGBaseAnalyzer(): E R R O R ! Calling default constructor without arguments!\n";
  std::cout << "                                  Please make sure to call the constructor that takes the\n";
  std::cout << "                                  class name as argument!\n";
  std::cout << "                                  Exiting...\n";
  exit(1);
}



MGBaseAnalyzer::MGBaseAnalyzer(const std::string name)
{
  className = name;
  baseInitialized = false;
}


void MGBaseAnalyzer::beginJob()
{
  if(baseInitialized)
    return;
  PDAnalyzerUtil::beginJob();

  getUserParameter( "verbose", verbose );
  
  int tries = 0;
  treeListName = getUserParameter("treeListName");
  evtSelection = getUserParameter("evtSelection");
  sampleName = treeListName;

  do
  {
    sampleName = sampleName.substr(sampleName.find("/") + 1);
  } 
  while(sampleName.find("/") != std::string::npos);
  
  sampleName = sampleName.substr(0, sampleName.find("."));
  
  histOutFileName = "his__" + className + "__" + sampleName + "__" + evtSelection + "__" + std::to_string(tries++) + ".root";
  std::cout << "MGBaseAnalyzer::beginJob(): setting variables:\n";
  std::cout << "                            className       = " << className << std::endl;
  std::cout << "                            treeListName    = " << treeListName << std::endl;
  std::cout << "                            sampleName      = " << sampleName << std::endl;
  std::cout << "                            evtSelection    = " << evtSelection << std::endl;
  std::cout << "                            histOutFileName = " << histOutFileName << std::endl;
  bool origNameNotOk = false;
  while(!gSystem->AccessPathName(gSystem->ExpandPathName(histOutFileName.c_str())))
  {
    origNameNotOk = true;
    if(tries == 1)
    {
      std::cout << "W A R N I N G! Output file \"" << histOutFileName << "\" already exists!\n";
      std::cout << "               Trying a different name...\n";
    }
    histOutFileName = "his__" + className + "__" + sampleName + "__" + evtSelection + "__" + std::to_string(tries++) + ".root";
    if(tries > 9999)
      break;
  }
  if(!gSystem->AccessPathName(gSystem->ExpandPathName(histOutFileName.c_str())))
  {
    std::cout << "E R R O R! Output file \"" << histOutFileName << "\" already exists!\n";
    std::cout << "           Exiting...\n";
    exit(1);
  }
  if(origNameNotOk)
    std::cout << "                            New histOutFileName = " << histOutFileName << std::endl;
  
  insertI( mAllNObjects,  nElectrons         );
  
  if(use_electrons_sphe)
  {
    insertVF(mAllEleFloats, elePt              );
    insertVF(mAllEleFloats, eleEta             );
    insertVF(mAllEleFloats, elePhi             );
    insertVF(mAllEleFloats, eleGsfPt           );
    insertVF(mAllEleFloats, eleGsfEta          );
    insertVF(mAllEleFloats, eleGsfPhi          );
    insertVF(mAllEleFloats, eleGsfPtAtVtx      );
    insertVF(mAllEleFloats, eleGsfEtaAtVtx     );
    insertVF(mAllEleFloats, eleGsfPhiAtVtx     );
    insertVF(mAllEleFloats, eleCtfPt           );
    insertVF(mAllEleFloats, eleCtfEta          );
    insertVF(mAllEleFloats, eleCtfPhi          );
  }
  if(use_electrons_cart)
  {
    insertVF(mAllEleFloats, elePx              );
    insertVF(mAllEleFloats, elePy              );
    insertVF(mAllEleFloats, elePz              );
    insertVF(mAllEleFloats, eleGsfPx           );
    insertVF(mAllEleFloats, eleGsfPy           );
    insertVF(mAllEleFloats, eleGsfPz           );
    insertVF(mAllEleFloats, eleGsfPxAtVtx      );
    insertVF(mAllEleFloats, eleGsfPyAtVtx      );
    insertVF(mAllEleFloats, eleGsfPzAtVtx      );
    insertVF(mAllEleFloats, eleCtfPx           );
    insertVF(mAllEleFloats, eleCtfPy           );
    insertVF(mAllEleFloats, eleCtfPz           );
  }
  insertVF(mAllEleFloats, eleGsfChi2         );
  insertVF(mAllEleFloats, eleGsfNdof         );
  insertVI(mAllEleInts,   eleGsfLostHits     );
  insertVI(mAllEleInts,   eleGsfValidHits    );
  insertVI(mAllEleInts,   eleGsfPVtx         );
  insertVI(mAllEleInts,   eleGsfQuality      );
  insertVI(mAllEleInts,   eleGsfHitPattern   );
  insertVI(mAllEleInts,   eleGsfLayPattern   );
  insertVF(mAllEleFloats, eleGsfPtError      );
  insertVF(mAllEleFloats, eleGsfDxy          );
  insertVF(mAllEleFloats, eleGsfDxyPV        );
  insertVF(mAllEleFloats, eleGsfDz           );
  insertVF(mAllEleFloats, eleGsfDzPV         );
  insertVF(mAllEleFloats, eleGsfExy          );
  insertVF(mAllEleFloats, eleGsfEz           );
  insertVI(mAllEleInts,   eleCtfTrk          );
  insertVF(mAllEleFloats, eleCtfChi2         );
  insertVF(mAllEleFloats, eleCtfNdof         );
  insertVI(mAllEleInts,   eleCtfLostHits     );
  insertVI(mAllEleInts,   eleCtfValidHits    );
  insertVI(mAllEleInts,   eleCtfPVtx         );
  insertVI(mAllEleInts,   eleCtfQuality      );
  insertVI(mAllEleInts,   eleCtfHitPattern   );
  insertVI(mAllEleInts,   eleCtfLayPattern   );
  insertVF(mAllEleFloats, eleCtfPtError      );
  insertVF(mAllEleFloats, eleCtfDxy          );
  insertVF(mAllEleFloats, eleCtfDxyPV        );
  insertVF(mAllEleFloats, eleCtfDz           );
  insertVF(mAllEleFloats, eleCtfDzPV         );
  insertVF(mAllEleFloats, eleCtfExy          );
  insertVF(mAllEleFloats, eleCtfEz           );
  insertVF(mAllEleFloats, eleE               );
  insertVI(mAllEleInts,   eleCharge          );
  insertVI(mAllEleInts,   eleTrg             );
  insertVF(mAllEleFloats, eleChaIso          );
  insertVF(mAllEleFloats, eleNeuIso          );
  insertVF(mAllEleFloats, elePhoIso          );
  insertVF(mAllEleFloats, elePCHIso          );
  insertVF(mAllEleFloats, eleSumCHpt         );
  insertVF(mAllEleFloats, eleSumCPpt         );
  insertVF(mAllEleFloats, eleSumNHet         );
  insertVF(mAllEleFloats, eleSumPHet         );
  insertVF(mAllEleFloats, eleSumPUpt         );
  insertVF(mAllEleFloats, eleAbsEta          );
  insertVF(mAllEleFloats, eleAEff            );
  insertVI(mAllEleInts,   eleID              );
  insertVB(mAllEleBools,  eleConsCha         );
  insertVB(mAllEleBools,  eleEBEEGap         );
  insertVF(mAllEleFloats, eleDb              );
  insertVI(mAllEleInts,   eleProps           );
  insertVF(mAllEleFloats, elePartIso         );
  insertVF(mAllEleFloats, eleEcalIso         );
  insertVF(mAllEleFloats, eleHcalIso         );
  insertVF(mAllEleFloats, eleTrackIso        );
  insertVI(mAllEleInts,   eleClass           );
  insertVF(mAllEleFloats, eleCtfGsfOv        );
  insertVF(mAllEleFloats, eleDEtaElClTkAtCalo);
  insertVF(mAllEleFloats, eleDEtaSeSlTkAtCalo);
  insertVF(mAllEleFloats, eleDEtaSeClTkAtVtx );
  insertVF(mAllEleFloats, eleDEtaSuClTkAtVtx );
  insertVF(mAllEleFloats, eleDPhiElClTkAtCalo);
  insertVF(mAllEleFloats, eleDPhiSeClTkAtCalo);
  insertVF(mAllEleFloats, eleDPhiSuClTkAtVtx );
  insertVF(mAllEleFloats, eleEcalE           );
  insertVF(mAllEleFloats, eleEcalEE          );
  insertVF(mAllEleFloats, eleEElClOvPOut     );
  insertVF(mAllEleFloats, eleESeClOvP        );
  insertVF(mAllEleFloats, eleESeClOvPOut     );
  insertVF(mAllEleFloats, eleESuClOvP        );
  insertVF(mAllEleFloats, eleFBrem           );
  insertVF(mAllEleFloats, eleHcalOvEcal      );
  insertVF(mAllEleFloats, eleF5x5E1x5        );
  insertVF(mAllEleFloats, eleF5x5E2x5B       );
  insertVF(mAllEleFloats, eleF5x5E2x5L       );
  insertVF(mAllEleFloats, eleF5x5E2x5M       );
  insertVF(mAllEleFloats, eleF5x5E2x5R       );
  insertVF(mAllEleFloats, eleF5x5E2x5T       );
  insertVF(mAllEleFloats, eleF5x5E5x5        );
  insertVF(mAllEleFloats, eleF5x5EB          );
  insertVF(mAllEleFloats, eleF5x5EL          );
  insertVF(mAllEleFloats, eleF5x5ER          );
  insertVF(mAllEleFloats, eleF5x5ET          );
  insertVF(mAllEleFloats, eleF5x5HcalOvEcal  );
  insertVF(mAllEleFloats, eleF5x5HcalOvEcalV );
  insertVF(mAllEleFloats, eleF5x5R9          );
  insertVF(mAllEleFloats, eleF5x5SEtaEta     );
  insertVF(mAllEleFloats, eleF5x5SIetaIeta   );
  insertVF(mAllEleFloats, eleF5x5SIphiIphi   );
  insertVF(mAllEleFloats, eleSCEtaWidth      );
  insertVF(mAllEleFloats, eleSCPhiWidth      );
  insertVF(mAllEleFloats, eleSCPreshowerE    );
  insertVF(mAllEleFloats, eleSCRawE          );
  insertVF(mAllEleFloats, eleSCEta           );
  insertVF(mAllEleFloats, eleConvDcot        );
  insertVF(mAllEleFloats, eleConvDist        );
  insertVI(mAllEleInts,   eleConvFlags       );
  insertVF(mAllEleFloats, eleConvRadius      );
  insertVF(mAllEleFloats, eleConvVtxProb     );
  
  baseInitialized = true;
}



void MGBaseAnalyzer::save(const std::string& oldHistName)
{
  // oldHistName is unused...
  TreeWrapper::save(histOutFileName);
}



void MGBaseAnalyzer::checkBranches()
{
  has_hltlist = checkBranch("nHLTPaths", "HLT path names", use_hltlist);
  has_hlts = checkBranch("nHLTStatus", "HLT status", use_hlts);
  has_hlto = checkBranch("nHLTObjects", "HLT objects", use_hlto);
  has_hltm = checkBranch("nHLTMatches", "HLT matches", use_hltm);
  has_bspot = checkBranch("bsX", "Beam spot", use_bspot);
  has_met = checkBranch("mEt", "Missing Et", use_met);
  has_muons = checkBranch("nMuons", "Muons", use_muons);
  has_electrons = checkBranch("nElectrons", "Electrons", use_electrons);
  has_taus = checkBranch("nTaus", "Taus", use_taus);
  has_jets = checkBranch("nJets", "Jets", use_jets);
  has_tags = checkBranch("nTags", "Jet tags", use_tags);
  has_info = checkBranch("nUserInfo", "User info", use_info);
  has_pflow = checkBranch("nPF", "Particle flow candidates", use_pflow);
  has_tracks = checkBranch("nTracks", "Tracks", use_tracks);
  has_trkper = checkBranch("nTrkPer", "Track perigee parameters", use_trkper);
  has_pvts = checkBranch("nPVertices", "Primary vertices", use_pvts);
  has_svts = checkBranch("nSVertices", "Secondary vertices", use_svts);
  has_vsub = checkBranch("nCompVts", "Composite particle vertices", use_vsub);
  has_tkips = checkBranch("nTkIPs", "Track impact parameters", use_tkips);
  has_vtxps = checkBranch("nVtxPs", "Track momenta at vertices", use_vtxps);
  has_puwgt = checkBranch("puWeight", "PU weight", use_puwgt);
  has_gen = checkBranch("nGenP", "Gen particles", use_gen);
  has_gpj = checkBranch("nGenJets", "Gen jets", use_gpj);  
}


