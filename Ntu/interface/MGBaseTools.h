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
    
    void setCurrentEvt(const int iEvt) {currentEvt = iEvt;}    
    
    bool verbose;

    // Index of current entry
    int currentEvt;
    
  private:
    // Nothing for the moment
    
};


#endif
