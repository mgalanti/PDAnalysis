#ifndef EleMVASecondNtupleWriter_h
#define EleMVASecondNtupleWriter_h

#include "EleMVASecondNtupleData.h"
#include "NtuTool/Common/interface/TreeWriter.h"
#include "TFile.h"
#include <string>

class EleMVASecondNtupleWriter: public EleMVASecondNtupleData, public TreeWriter {

 public:

  EleMVASecondNtupleWriter() {
  }
  virtual ~EleMVASecondNtupleWriter() {
  }

  void open( const std::string& name, const std::string mode = "CREATE" ) {
    TDirectory* current = gDirectory;
    file = new TFile( name.c_str(), mode.c_str() );
    initTree();
    initWrite( file );
    current->cd();
    return;
  }
 
  void close() {
    TreeWriter::close();
    file->Close();
  }
  
  void Reset() {
    reset();
  }

 private:

  TFile* file;

  EleMVASecondNtupleWriter           ( const EleMVASecondNtupleWriter& a );
  EleMVASecondNtupleWriter& operator=( const EleMVASecondNtupleWriter& a );

};

#endif

