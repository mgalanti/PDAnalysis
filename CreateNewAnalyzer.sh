#!/bin/bash

if [[ $1 == '-h' || $1 == '--help' ]]
  then echo Instructions:
  # FIXME: to be written
fi

if [[ $1 == '' || $2 != '' ]]
  then 
    echo Script accepts exactly one argument!
    echo Usage:
    echo   $0 name
    echo Write \"$0 --help\" for more help.
  exit
fi

# Some common words that can be present in the analyzer name, to use in string substitution below
declare -a -r MainStrings=("Analyzer" "Reader" "Writer" "Producer")
# declare -a -r LCStrings=("analyzer" "reader" "writer" "producer")
declare -a -r VerbStrings=("analyze" "read" "write" "produce")
declare -a -r NounStrings=("Analysis" "Reading" "Writing" "Production")

name=$1
LCName=${name,}

declare -i i
i=0
found=false
for MainString in "${MainStrings[@]}"
do
  if [[ $name == *${MainString} ]]
  then
    found=true
    substrLen=${#MainStrings[i]}
    strippedName=${name::${#name}-$substrLen}
    verbName=${VerbStrings[i]}$strippedName
    nounName=$strippedName${NounStrings[i]}
  fi     
  i+=1
done

if [[ "$found" == false ]]
then
  verbName="analyze"$name
  nounName=$name"Analysis"
fi

echo Name of the analyzer to be created: $name
echo Verb name: $verbName
echo Noun name: $nounName
echo lower-case name: $LCName

echo Creating directory $name
if [[ -e $name ]]
  then echo ERROR: directory $name already exists!
  echo Aborting analyzer creation!
  exit
fi
mkdir $name

echo Creating file $name/BuildFile.xml
if [[ -e $name/BuildFile.xml ]]
  then echo ERROR: analyzer build file $name/BuildFile.xml already exists!
  echo Aborting analyzer creation!
  exit
fi
cp templates/Ntu_BuildFile_xml.template $name/BuildFile.xml

echo Creating directory $name/bin
if [[ -e $name/bin ]]
  then echo ERROR: directory $name/bin already exists!
  echo Aborting analyzer creation!
  exit
fi
mkdir $name/bin

echo Creating file $name/bin/$name.cc
if [[ -e $name/bin/$name.cc ]]
  then echo ERROR: analyzer source file $name/bin/$name.cc already exists!
  echo Aborting analyzer creation!
  exit
fi
sed "s#ANALYZERNAME#$name#g" templates/Ntu_bin_Analyzer_cc.template > $name/bin/$name.cc

echo Creating file $name/bin/$name.h
if [[ -e $name/bin/$name.h ]]
  then echo ERROR: analyzer header file $name/bin/$name.h already exists!
  echo Aborting analyzer creation!
  exit
fi
sed "s#ANALYZERNAME#$name#g" templates/Ntu_bin_Analyzer_h.template > $name/bin/$name.h

echo Creating file $name/bin/BuildFile.xml
if [[ -e $name/bin/BuildFile.xml ]]
  then echo ERROR: analyzer build file $name/bin/BuildFile.xml already exists!
  echo Aborting analyzer creation!
  exit
fi
sed -e "s#ANALYZERNAME#$name#g" -e "s#VERBNAME#$verbName#g" templates/Ntu_bin_BuildFile_xml.template > $name/bin/BuildFile.xml

echo Creating file $name/bin/$verbName.cc
if [[ -e $name/bin/$verbName.cc ]]
  then echo ERROR: analyzer source file $name/bin/$verbName.cc already exists!
  echo Aborting analyzer creation!
  exit
fi
sed -e "s#ANALYZERNAME#$name#g" -e "s#ANALYSISNAME#$nounName#g" -e "s#ANALYZERLCNAME#$LCName#g" templates/Ntu_bin_analyze_cc.template > $name/bin/$verbName.cc

echo Creating file $name/bin/$verbName.cfg
if [[ -e $name/bin/$verbName.cfg ]]
  then echo ERROR: analyzer configuration file $name/bin/$verbName.cfg already exists!
  echo Aborting analyzer creation!
  exit
fi
cp templates/Ntu_bin_analyze_cfg.template $name/bin/$verbName.cfg

echo Creating file $name/bin/${name}FW.cc
if [[ -e $name/bin/${name}FW.cc ]]
  then echo ERROR: analyzer source file $name/bin/${name}FW.cc already exists!
  echo Aborting analyzer creation!
  exit
fi
sed "s#ANALYZERNAME#$name#g" templates/Ntu_bin_AnalyzerFW_cc.template > $name/bin/${name}FW.cc

echo Creating file $name/bin/compile.sh
if [[ -e $name/bin/compile.sh ]]
  then echo ERROR: analyzer compilation script file $name/bin/compile.sh already exists!
  echo Aborting analyzer creation!
  exit
fi
sed "s#VERBNAME#$verbName#g" templates/Ntu_bin_compile_sh.template > $name/bin/compile.sh

echo Creating directory $name/src
if [[ -e $name/src ]]
  then echo ERROR: directory $name/src already exists!
  echo Aborting analyzer creation!
  exit
fi
mkdir $name/src

echo Creating file $name/src/${name}FW.cc
if [[ -e $name/src/${name}FW.cc ]]
  then echo ERROR: analyzer source file $name/src/${name}FW.cc already exists!
  echo Aborting analyzer creation!
  exit
fi
sed "s#ANALYZERNAME#$name#g" templates/Ntu_src_AnalyzerFW_cc.template > $name/src/${name}FW.cc
