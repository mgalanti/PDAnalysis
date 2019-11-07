#!/bin/bash

timeNow=`date +%Y%m%d_%H%M%S`

logFile=produceEleMVACalibration_${timeNow}.log

# This automatically copies script's output to a log file
exec > $logFile
exec 2>&1

inputCfgFile=produceEleMVACalibration.cfg
nEvents=-1
verbose=1

while getopts ":c:n:qh" opt; do
  case $opt in
    h)
      echo "Usage: produceEleMVACalibration.sh [-c inputCfgFile] [-n nEvents] [-q] [-h]"
      echo ""
      echo "Meaning of optional arguments:"
      echo "-c: uses a non-standard name for the input cfg file."
      echo "    [default: produceEleMVACalibration.cfg]"
      echo "-n: allows to run on a given number of events."
      echo "    [default: runs on all events]"
      echo "-q: quiet execution. Switches off all verbose messages."
      echi "    [default: verbose messages are on.]"
      echo "-h: shows this help."
      exit 0
    ;;
    c)
#       echo "c was triggered! parameter $OPTARG"
      inputCfgFile=$OPTARG
    ;;
    n)
#       echo "n was triggered! parameter $OPTARG"
      nEvents=$OPTARG
    ;;
    q)
      verbose=0
    ;;
    \?) 
      echo "Error! Invalid option -$OPTARG" >&2
      exit 1
    ;;
    :)
      echo "Error! Option -$OPTARG requires an argument." >&2
      exit 1  
    ;;
  esac
done

cfgFileBaseName=`echo ${inputCfgFile} | rev | cut -d '.' -f 2- | rev`

# echo inputCfgFile = ${inputCfgFile}
# echo cfgFileBaseName = ${cfgFileBaseName}
# echo nEvents = ${nEvents}

if [[ ${verbose} == 1 ]]
then
  echo Modifying cfg file name...
fi
  
cfgFile=${cfgFileBaseName}_${timeNow}.cfg

if [[ ${verbose} == 1 ]]
then
  echo "     cfgFile = $cfgFile"
fi 

cp "$inputCfgFile" "$cfgFile"

if [[ ${verbose} == 1 ]]
then
  echo ""
  echo Parsing cfg file...
fi

treeListFileWithDir=""
sampleName=""
mvaMethod=""
mvaInputPath=""

while IFS= read -r inputLine
do
  if [[ ${inputLine} == treeListName* ]]
  then
    firstItem=`echo ${inputLine} | awk '{print $1}'`
    secondItem=`echo ${inputLine} | awk '{print $2}'`
    treeListFileWithDir=${secondItem}
    treeListFile=`echo ${secondItem} | awk -F "/" '{print $NF}'`
    sampleName=`echo ${treeListFile} | cut -d "_" -f 2- | rev | cut -d "." -f 2- | rev | cut -c 2-`
  elif [[ ${inputLine} == mvaMethod* ]]
  then
    firstItem=`echo ${inputLine} | awk '{print $1}'`
    secondItem=`echo ${inputLine} | awk '{print $2}'`
    mvaMethod=${secondItem}
  elif [[ ${inputLine} == mvaInputPath* ]]
  then
    firstItem=`echo ${inputLine} | awk '{print $1}'`
    secondItem=`echo ${inputLine} | awk '{print $2}'`
    mvaInputPath=${secondItem}
  fi 
done < $cfgFile

if [[ ${verbose} == 1 ]]
then
  echo "     cfg file parsed." 
fi
  
if [[ $treeListFileWithDir == "" ]]
then
  echo Error! Could not parse treeListName entry in configuration!
  echo Please fix the cfg file...
  exit 1
fi

if [[ $sampleName == "" ]]
then
  echo Error! Could not retrieve sampleName from configuration!
  echo Please fix the cfg file...
  exit 1
fi

if [[ $mvaMethod == "" ]]
then
  echo Error! Could not parse mvaMethod entry in configuration!
  echo Please fix the cfg file...
  exit 1
fi

if [[ $mvaInputPath == "" ]]
then
  echo Error! Could not parse mvaInputPath entry in configuration!
  echo Please fix the cfg file...
  exit 1
fi

if [[ ${verbose} == 1 ]]
then
  echo All needed variables retrieved:
  echo "     "treeListFileWithDir = $treeListFileWithDir
  echo "     "treeListFile = $treeListFile
  echo "     "sampleName = $sampleName
  echo "     "mvaMethod = $mvaMethod
  echo "     "mvaInputPath = $mvaInputPath
  echo ""
fi

# Create directory to host the job output
outDir=`echo jobs/job_${timeNow}`
if [[ ${verbose} == 1 ]]
then
  echo Creating job output directory \"${outDir}\"...
fi
mkdir $outDir

# The job is launched and run from the bin/ directory, however its inputs are copied to
# the job directory for future reference
if [[ ${verbose} == 1 ]]
then
  echo Copy input files to job output directory...
  echo "     "cp "$cfgFile" "$outDir"
fi
cp "$cfgFile" "$outDir"

# Now let's make sure that the correct weeights are used
if [[ ${verbose} == 1 ]]
then
  echo Sym-linking MVA weights directory...
fi

if [[ -L dataset/weights ]]
then
  if [[ ${verbose} == 1 ]]
  then
    echo "     "Found existing directory \"dataset/weights\". It is a symlink to:
    echo "          "\"`readlink -f dataset/weights`\".
    echo "     "Unlinking it...
  fi
  unlink dataset/weights
fi

if [[ -f dataset/weights ]]
then
  echo Error! Directory \"dataset/weights\" already exists and is not a symlink!
  echo Please remove it before launching the calibration...
  exit 1
fi

if [[ ${verbose} == 1 ]]
then
  echo Checking existence of mvaInputPath directory
fi

if [[ ! -e $mvaInputPath ]]
then
  echo Error! MVA input path does not exist!
  exit 1
fi

if [[ ${verbose} == 1 ]]
then
  echo Creating symbolic link to mva input path...
  echo "     "ln -s "`realpath ${mvaInputPath}`" dataset/weights
fi
ln -s "`realpath ${mvaInputPath}`" dataset/weights

clOptions="${treeListFileWithDir} his_produceEleMVACalibration__${sampleName}_${timeNow}.root -c ${cfgFile}" 

if [[ $nEvents > 0 ]]
then
  clOptions+=" -n ${nEvents}"
fi

# # produceEleMVACalibration ${treeListFileWithDir} his_produceEleMVACalibration__${sampleName}.root -c ${cfgFile} -n ${nEvents} > produceEleMVACalibration.log

echo ""
echo Launching calibration job...
echo ""
if [[ ${verbose} == 1 ]]
then
  echo "produceEleMVACalibration ${clOptions}"
fi
produceEleMVACalibration ${clOptions}

echo ""
echo Calibration job finished.

echo Copying output files to output directory...
# Let's first retrieve the process name used in this job...
# The logic below is the same as that used inside the executable
if [[ ${verbose} == 1 ]]
then
  echo Retrieving process name from job configuration...
fi

processName=""
if [[ "$treeListFileWithDir" == *"Bs"* ]]
then
  process="BsJPsiPhi"
  if [[ "$treeListFileWithDir" == *"DG0"* ]]
  then
    process+="DG0"
  fi
elif [[ "$treeListFileWithDir" == *"Bu"* ]]
then
  process="BuJPsiK"
fi

if [[ "$treeListFileWithDir" == *"MC"* ]]
then
  process+="MC"
elif [[ "$treeListFileWithDir" == *"Data"* ]]
then
  process+="Data"
fi

if [[ "$treeListFileWithDir" == *"2017"* ]]
then
  process+="2017"
elif [[ "$treeListFileWithDir" == *"2018"* ]]
then
  process+="2018"
fi

if [[ ${verbose} == 1 ]]
then
  echo "     "Process name retrieved: \"${process}\"
fi

calPlotFileName="calibration_${process}.pdf"
dnnPlotFileName="dnnDistribution_${process}.pdf"
calRootFileName="OSElectronTaggerCalibration${process}.root"
hisRootFileName="his_produceEleMVACalibration__${sampleName}_${timeNow}.root"

if [[ ${verbose} == 1 ]]
then
  echo Moving calibration plot file to output directory...
  echo "     "mv \"${calPlotFileName}\" \"${outDir}\"
fi
mv "${calPlotFileName}" "${outDir}"

if [[ ${verbose} == 1 ]]
then
  echo Moving dnn distribution plot file to output directory...
  echo "     "mv \"${dnnPlotFileName}\" \"${outDir}\"
fi
mv "${dnnPlotFileName}" "${outDir}"

if [[ ${verbose} == 1 ]]
then
  echo Moving calibration ROOT file to output directory...
  echo "     "mv \"${calRootFileName}\" \"${outDir}\"
fi
mv "${calRootFileName}" "${outDir}"

if [[ ${verbose} == 1 ]]
then
  echo Moving histogram ROOT file to output directory...
  echo "     "mv \"${hisRootFileName}\" \"${outDir}\"
fi
mv "${hisRootFileName}" "${outDir}"

if [[ ${verbose} == 1 ]]
then
  echo Removing input configuration file from execution directory...
  echo "     "rm \"${cfgFile}\"
fi
if [[ -e ${cfgFile} ]]
then
  if [[ ${cfgFile} != ${inputCfgFile} ]]
  then
    rm "${cfgFile}"
  else
    echo "Warning! Cfg file name is the same as input cfg file name! Not removing it."
  fi
else
  echo "Warning! Cfg file does not exist anymore!"
fi

if [[ ${verbose} == 1 ]]
then
  echo Moving log file to output directory...
  if [[ -e ${logFile} ]]
  then 
    echo "     "cp \"${logFile}\" \"${outDir}\"
    echo "     "rm \"${logFile}\"
  fi
fi

if [[ -e ${logFile} ]]
then 
  cp "${logFile}" "${outDir}"
  sleep 2
  rm "${logFile}"
fi

exit 0
