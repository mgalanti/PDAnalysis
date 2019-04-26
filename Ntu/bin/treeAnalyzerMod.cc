#define UTIL_USE FULL
#include "PDAnalyzer.h"
#include "NtuAnalysis/Read/interface/NtuFlexibleAnalyzer.h"
#include "NtuTool/Read/bin/treeAnalyze.cc"

#include "NtuAnalysis/Read/interface/NtuReader.h"
#include "NtuAnalysis/Read/interface/NtuEDMReader.h"

/** \class NtuModifiedAnalyzer
 *
 *  Description: 
 *    Modified class to steer the tree analysis letting 
 *          - treeListName set from cnfig.file
 *          - access it from PDAnalyzer
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "NtuTool/Read/interface/TreeStandardAnalyzer.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "NtuTool/Read/interface/TreeReader.h"
#include "TChain.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
class TreeReader;

//              ---------------------
//              -- Class Interface --
//              ---------------------

template<class T>
class NtuModifiedAnalyzer: public TreeStandardAnalyzer {

 public:

  /** Constructor
   */
  NtuModifiedAnalyzer( const std::string& process,
                       const std::string& producer ):
   ntuProdProcess ( process  ),
   ntuProdProducer( producer ) {
  }

  /** Destructor
   */
  virtual ~NtuModifiedAnalyzer() {}

  /** Operations
   */
  /// run the application
  virtual int run( int argc, char* argv[] ) {

    std::cout << "NtuModifiedAnalyzer::run" << std::endl;

    char** argp = argv;
    char** argl = argp + argc;
    std::string key = "";
    std::string val = "";
    enum ntuType { ntu, edm };
    ntuType type = ntu;
    while ( argp < argl ) {
      std::string args( *argp++ );
      if ( args == "-v" ) {
        key = *argp++;
        val = *argp++;
        if ( key == "process"  ) ntuProdProcess  = val;
        if ( key == "producer" ) ntuProdProducer = val;
        if ( key == "ntuType"  ) {
          if      ( val == "ntu" ) type = ntu;
          else if ( val == "edm" ) type = edm;
          else                     std::cout << "invalid ntuple type: " << val
                                             << " choose \"ntu\" or \"edm\""
                                             << std::endl;
        }
      }
    }
    
    if ( type == ntu ) new NtuReader<T>;
    if ( type == edm ) new NtuEDMReader<T>( ntuProdProcess, ntuProdProducer );

    TreeReader* tr = TreeReader::getInstance();
    std::string treeListName( "treeList" );
    std::string histFileName( "hist.root" );
    bool listMissing = true;
    bool histMissing = true;
    char** argq = argv;
    char** argj = argq++ + argc;
    char** argn = new char*[argc + 2];
    char** argm = argn;
    *argm++ = *argv;
    while ( argq < argj ) {
      char* cptr = *argq++;
      const std::string args( cptr );
      if ( *cptr == '-' ) *argm++ = cptr;
      if ( args == "-n" ) {
	tr->setUserParameter( "-n" , *argm++ = *argq++ );
        continue;
      }
      if ( args == "-s" ) {
	tr->setUserParameter( "-s" , *argm++ = *argq++ );
        continue;
      }
      if ( args == "-a" ) {
	tr->setUserParameter( "-a" , *argm++ = *argq++ );
        continue;
      }
      if ( args == "-c" ) {
        tr->setConfiguration( *argq );
	tr->setUserParameter( "-c" , *argm++ = *argq++ );
        continue;
      }
      if ( args == "-v" ) {
        key = *argm++ = *argq++;
        val = *argm++ = *argq++;
        tr->setUserParameter( key, val );
        continue;
      }
      const std::string& tlnConf = tr->getUserParameter( "treeListName" );
      const std::string& hfnConf = tr->getUserParameter( "histFileName" );
      if ( tlnConf.length() )     treeListName = tlnConf;
      else tr->setUserParameter( "treeListName", treeListName );
      if ( hfnConf.length() )     histFileName = hfnConf;
      else tr->setUserParameter( "histFileName", histFileName );
      if ( listMissing ) {
        tr->setUserParameter( "treeListName", args );
        listMissing = false;
        continue;
      }
      if ( histMissing ) {
        tr->setUserParameter( "histFileName", args );
        histMissing = false;
        continue;
      }
    }
    if ( listMissing ) ++argc;
    if ( histMissing ) ++argc;
    const std::string& tlnConf = tr->getUserParameter( "treeListName" );
//    if ( histMissing && tlnConf.length() ) histFileName = treeListName;
    const std::string& hfnConf = tr->getUserParameter( "histFileName" );
    if ( tlnConf.length() )     treeListName = tlnConf;
//    else tr->setUserParameter( "treeListName", treeListName );
    if ( hfnConf.length() )     histFileName = hfnConf;
//    else tr->setUserParameter( "histFileName", histFileName );

    const char* sptr;
    char* dptr;
    sptr = treeListName.c_str();
    dptr = *argm++ = new char[treeListName.length() + 1];
    while ( ( *dptr++ = *sptr++ ) );
    sptr = histFileName.c_str();
    dptr = *argm++ = new char[histFileName.length() + 1];
    while ( ( *dptr++ = *sptr++ ) );

    return TreeStandardAnalyzer::run( argc, argn );

  }

 private:

  NtuModifiedAnalyzer( const NtuModifiedAnalyzer& t );
  NtuModifiedAnalyzer& operator=( const NtuModifiedAnalyzer& t );

  std::string ntuProdProcess;
  std::string ntuProdProducer;

};



static NtuModifiedAnalyzer<PDAnalyzer> nfa( "pdAnalysis", "pdAnalyzer" );
