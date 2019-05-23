#ifndef MGBaseTools_h
#define MGBaseTools_h


#include "PDAnalysis/Ntu/interface/constants.h"
#include "PDAnalysis/Ntu/interface/PDAnalyzerUtil.h"


class MGBaseTools :
      public virtual PDAnalyzerUtil
{
  public:
    virtual ~MGBaseTools() {};
    
    virtual void beginJob();
    
    bool verbose;
    
  private:
    // Nothing for the moment
    
};


#endif
