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
    std::cout << "W A R N I N G! Output file \"" << histOutFileName << "\" already exists!\n";
    std::cout << "               Trying a different name...\n";
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
