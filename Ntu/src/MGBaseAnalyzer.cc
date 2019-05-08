#define UTIL_USE FULL

#include "TSystem.h"

#include "PDAnalysis/Ntu/interface/MGBaseAnalyzer.h"



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
  has_pvts = checkBranch("nPVertices", "Primary vertices", use_pvts);
  has_svts = checkBranch("nSVertices", "Secondary vertices", use_svts);
  has_vsub = checkBranch("nCompVts", "Composite particle vertices", use_vsub);
  has_tkips = checkBranch("nTkIPs", "Track impact parameters", use_tkips);
  has_vtxps = checkBranch("nVtxPs", "Track momenta at vertices", use_vtxps);
  has_puwgt = checkBranch("puWeight", "PU weight", use_puwgt);
  has_gen = checkBranch("nGenP", "Gen particles", use_gen);
  has_gpj = checkBranch("nGenJets", "Gen jets", use_gpj);  
}


