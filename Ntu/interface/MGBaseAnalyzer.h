#ifndef MGBaseAnalyzer_h
#define MGBaseAnalyzer_h



#include "PDAnalysis/Ntu/interface/PDAnalyzerUtil.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"


class MGBaseAnalyzer : public virtual PDAnalyzerUtil
{
  public:
    MGBaseAnalyzer();
    MGBaseAnalyzer(const std::string name);
    virtual ~MGBaseAnalyzer() {};
    
    virtual void beginJob();
    
    // Override here a method of the TreeWrapper in order to set a custom output file name from here
    virtual void save(const std::string& oldHistName);

    // Used in derived classes
    std::string evtSelection;
    
  private:
    // Poor man's introspection... Enough for my purposes (dynamically set the output name based on tha analyzer name)
    std::string className;

    std::string sampleName;    
    std::string treeListName;    
    std::string histOutFileName;
};



#endif // MGBaseAnalyzer_h
