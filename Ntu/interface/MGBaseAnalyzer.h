#ifndef MGBaseAnalyzer_h
#define MGBaseAnalyzer_h



// #include "PDAnalysis/Ntu/interface/PDAnalyzerUtil.h"
#include "PDAnalysis/Ntu/interface/MGGenTools.h"
#include "PDAnalysis/Ntu/interface/MGRecoTools.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"


class MGBaseAnalyzer: 
    public virtual MGGenTools,
    public virtual MGRecoTools,
    public virtual MGSelector
//     public virtual PDAnalyzerUtil
    
{
  public:
    MGBaseAnalyzer();
    MGBaseAnalyzer(const std::string name);
    virtual ~MGBaseAnalyzer() {};
    
    virtual void beginJob();
    
    // Override here a method of the TreeWrapper in order to set a custom output file name from here
    virtual void save(const std::string& oldHistName);
    
    // Method to check whether a given collection is present in the input tree
    inline bool checkBranch(const std::string branchName, const std::string humanReadableName, const bool useFlag)
    {
      bool branchInTree = currentTree->GetBranchStatus(branchName.c_str());
      if(!branchInTree && useFlag)
      {
        if(useFlag)
        {
          std::cout << "W A R N I N G ! " << humanReadableName << " collection is not present in input tree but its corresponding use flag is set to true!\n";
        }
        else if(verbose)
        {
        std::cout << "I N F O . " << humanReadableName << " collection is not present in input tree.\n";
        }
      }
      if(branchInTree && !useFlag)
      {
        std::cout << "W A R N I N G ! " << humanReadableName << " collection is present in input tree but its corresponding use flag is set to false!\n";
      }
      return branchInTree;
    }
    
    void checkBranches();
    
    // Used in derived classes
    std::string evtSelection;
    
    // Poor man's introspection... Enough for my purposes (dynamically set the output name based on tha analyzer name)
    std::string className;

    std::string sampleName;    
    std::string treeListName;    
    std::string histOutFileName;
    
    // Status of various branches in input tree
    bool has_hltlist;
    bool has_hlts;
    bool has_hlto;
    bool has_hltm;
    bool has_bspot;
    bool has_met;
    bool has_muons;
    bool has_electrons;
    bool has_taus;
    bool has_jets;
    bool has_tags;
    bool has_info;
    bool has_pflow;
    bool has_tracks;
    bool has_trkper;
    bool has_pvts;
    bool has_svts;
    bool has_vsub;
    bool has_tkips;
    bool has_vtxps;
    bool has_puwgt;
    bool has_gen;
    bool has_gpj;
    
    std::map<std::string, int*> mAllNObjects;

    std::map<std::string, std::vector<float>**> mAllEleFloats;
    std::map<std::string, std::vector<int>**> mAllEleInts;
    std::map<std::string, std::vector<bool>**> mAllEleBools;

    

  private:
    bool baseInitialized;
    
};



#endif // MGBaseAnalyzer_h
