#!/bin/bash

eventsPerJob=$1
nJobs=$2

timeNow=`date +%Y%m%d_%H%M%S`

libDir=${CMSSW_BASE}/lib/${SCRAM_ARCH}
srcDir=${CMSSW_BASE}/src/PDAnalysis
exeDir=${CMSSW_BASE}/bin/${SCRAM_ARCH}
jobDir=${srcDir}/ElectronEfficiencyAnalyzer/bin/analyzeElectronEfficiency_job_${timeNow}

libFile=lib_${timeNow}.tar.gz
cfgFile=analyzeElectronEfficiency_${timeNow}.cfg
scriptFile=analyzeElectronEfficiency_${timeNow}.sh
jobConfigFile=analyzeElectronEfficiency_${timeNow}.jobconfig

# Create directory tree to store job input/output files
mkdir ${jobDir}
mkdir ${jobDir}/his
mkdir ${jobDir}/out

# Compress and save custom libs to a tar.gz file
# Needed only when not using a filesystem shared between master and slaves
echo Creating tar file ${libFile} with custom libs...
tar -czvf ${jobDir}/${libFile} -C ${libDir}/ .

# Read analyzer configuration file and modify it to be used with condor
echo Creating analyzer configuration file ${jobDir}/${cfgFile}  ...
touch ${jobDir}/${cfgFile}

hasSecondNtuple=false

while IFS= read -r inputLine
do
  if [[ ${inputLine} == treeListName* ]]
  then
    firstItem=`echo ${inputLine} | awk '{print $1}'`
    secondItem=`echo ${inputLine} | awk '{print $2}'`
    treeListFileWithDir=${secondItem}
    treeListFile=`echo ${secondItem} | awk -F "/" '{print $NF}'`
    outputLine=`echo ${firstItem} ${treeListFile}`
  elif [[ ${inputLine} == evtSelection* ]]
  then
    evtSelection=`echo ${inputLine} | awk '{print $2}'`
    outputLine=${inputLine}
  elif [[ ${inputLine} == secondNtupleBaseName* ]]
  then
    hasSecondNtuple=true
    secondNtupleBaseName=`echo ${inputLine} | awk '{print $2}'`
    outputLine=${inputLine}
  else
    outputLine=${inputLine}
  fi
  echo $outputLine >> ${jobDir}/${cfgFile}
done < analyzeElectronEfficiency.cfg

if [[ ${hasSecondNtuple} == "true" ]]
then
  mkdir ${jobDir}/ntu
  outFiles="his,ntu"
else
  outFiles="his"
fi

sampleName="${treeListFile%.*}"
echo treeListFileWithDir is ${treeListFileWithDir}
echo treeListFile is ${treeListFile}
echo sampleName is ${sampleName}
echo evtSelection is ${evtSelection}

# Create jobconfig from template inside the job directory
echo Creating file ${jobDir}/${jobConfigFile}
sed -e "s#JOBDIR#${jobDir}#g" -e "s#SCRIPTFILE#${scriptFile}#g" -e "s#EVENTSPERJOB#${eventsPerJob}#g" -e "s#LIBFILE#${libFile}#g" -e "s#EXEDIR#${exeDir}#g" -e "s#CFGFILE#${cfgFile}#g" -e "s#TREELISTFILE#${treeListFile}#g" -e "s#OUTFILES#${outFiles}#g" -e "s#NJOBS#${nJobs}#g" analyzeElectronEfficiency_template.jobconfig > ${jobDir}/${jobConfigFile}

# Create shell script executable from template inside the job directory
echo Creating file ${jobDir}/${scriptFile}
sed -e "s#LIBFILE#${libFile}#g" -e "s#CFGFILE#${cfgFile}#g" -e "s#SAMPLENAME#${sampleName}#g" -e "s#EVTSELECTION#${evtSelection}#g" -e "s#SECONDNTUPLEBASENAME#${secondNtupleBaseName}#g" -e "s#TIMENOW#${timeNow}#g" analyzeElectronEfficiency_template.sh > ${jobDir}/${scriptFile}

# Copy other needed files to job directory
echo Copying tree list file ${treeListFileWithDir} to job directory
cp ${treeListFileWithDir} ${jobDir}/${treeListFile}

echo Copying executable ${exeDir}/analyzeElectronEfficiency to job directory
cp ${exeDir}/analyzeElectronEfficiency ${jobDir}/analyzeElectronEfficiency
