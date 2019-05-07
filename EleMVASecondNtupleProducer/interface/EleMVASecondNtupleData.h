#ifndef EleMVASecondNtupleData_h
#define EleMVASecondNtupleData_h

#include "NtuTool/Common/interface/TreeWrapper.h"

class EleMVASecondNtupleData: public virtual TreeWrapper {

 public:

  EleMVASecondNtupleData() {
  }
  virtual ~EleMVASecondNtupleData() {
  }

  void initTree() {
    treeName = "EleMVAsecondTree";
    
    // Signal-side variables
    setBranch( "tightB", &tightB, "tightB/b", &b_tightB );
    
    // Signal-side variables
    setBranch( "elePt", &elePt, "elePt/F", &b_elePt );
    setBranch( "eleEta", &eleEta, "eleEta/F", &b_eleEta );
    setBranch( "elePhi", &elePhi, "elePhi/F", &b_elePhi );
    
  }

  unsigned short tightB;
  
  float elePt, eleEta, elePhi;

  TBranch *b_tightB;
  TBranch *b_elePt, *b_eleEta, *b_elePhi;
  
 private:

  EleMVASecondNtupleData           ( const EleMVASecondNtupleData& a );
  EleMVASecondNtupleData& operator=( const EleMVASecondNtupleData& a );

};

#endif

