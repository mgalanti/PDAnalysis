#!/bin/bash

eventsPerJob=$1
nJobs=$2

timeNow=`date +%Y%m%d_%H%M%S`

libDir=${CMSSW_BASE}/lib/${SCRAM_ARCH}
srcDir=${CMSSW_BASE}/src/PDAnalysis
exeDir=${CMSSW_BASE}/bin/${SCRAM_ARCH}
jobDir=${srcDir}/EleMVASecondNtupleProducer/bin/produceEleMVASecondNtuple_job_${timeNow}

libFile=lib_${timeNow}.tar.gz
cfgFile=produceEleMVASecondNtuple_${timeNow}.cfg
scriptFile=produceEleMVASecondNtuple_${timeNow}.sh
jobConfigFile=produceEleMVASecondNtuple_${timeNow}.jobconfig

# Create directory to store job input/output files
mkdir ${jobDir}

# Compress and save custom libs to a tar.gz file
# Needed only when not using a filesystem shared between master and slaves
echo Creating tar file ${libFile} with custom libs...
tar -czvf ${jobDir}/${libFile} -C ${libDir}/ .

# Read analyzer configuration file and modify it to be used with condor
echo Creating analyzer configuration file ${jobDir}/${cfgFile}  ...
touch ${jobDir}/${cfgFile}
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
    secondNtupleBaseName=`echo ${inputLine} | awk '{print $2}'`
    outputLine=${inputLine}
  else
    outputLine=${inputLine}
  fi
  echo $outputLine >> ${jobDir}/${cfgFile}
done < produceEleMVASecondNtuple.cfg

sampleName="${treeListFile%.*}"
echo treeListFileWithDir is ${treeListFileWithDir}
echo treeListFile is ${treeListFile}
echo sampleName is ${sampleName}
echo evtSelection is ${evtSelection}

# Create jobconfig from template inside the job directory
echo Creating file ${jobDir}/${jobConfigFile}
sed -e "s#JOBDIR#${jobDir}#g" -e "s#SCRIPTFILE#${scriptFile}#g" -e "s#EVENTSPERJOB#${eventsPerJob}#g" -e "s#LIBFILE#${libFile}#g" -e "s#EXEDIR#${exeDir}#g" -e "s#CFGFILE#${cfgFile}#g" -e "s#TREELISTFILE#${treeListFile}#g" -e "s#NJOBS#${nJobs}#g" produceEleMVASecondNtuple_template.jobconfig > ${jobDir}/${jobConfigFile}

# Create shell script executable from template inside the job directory
echo Creating file ${jobDir}/${scriptFile}
sed -e "s#LIBFILE#${libFile}#g" -e "s#CFGFILE#${cfgFile}#g" -e "s#SAMPLENAME#${sampleName}#g" -e "s#EVTSELECTION#${evtSelection}#g" -e "s#SECONDNTUPLEBASENAME#${secondNtupleBaseName}#g" -e "s#TIMENOW#${timeNow}#g" produceEleMVASecondNtuple_template.sh > ${jobDir}/${scriptFile}

# Copy other needed files to job directory
echo Copying tree list file ${treeListFileWithDir} to job directory
cp ${treeListFileWithDir} ${jobDir}/${treeListFile}
