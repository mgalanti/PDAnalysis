#!/bin/bash

echo Extracting custom library files...
filesInLib=`tar tf LIBFILE`
tar xvfz LIBFILE

echo " "

echo Setting environment: export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH

echo " "

echo Setting job parameters...
eventsPerJob=$1
skipEvents=$(( $2 * $1 ))
lastEvent=$(( ${skipEvents} + ${eventsPerJob} ))
echo     eventsPerJob is $eventsPerJob
echo     skipEvents is $skipEvents
echo     lastEvent is $lastEvent

echo " "

echo Running analyzer: ./produceEleMVASecondNtuple -c CFGFILE -n $eventsPerJob -s $skipEvents
./produceEleMVASecondNtuple -c CFGFILE -n $eventsPerJob -s $skipEvents

echo " "

echo Cleaning up files...
for filename in `echo $filesInLib`
do 
  if [[ $filename != ./ ]]
  then 
    rm -v $filename
  fi
done

echo " "

echo Moving histogram output file: mv his__EleMVASecondNtupleProducer__SAMPLENAME__EVTSELECTION__0.root his/his__EleMVASecondNtupleProducer__SAMPLENAME__EVTSELECTION__TIMENOW__${skipEvents}__${lastEvent}.root
mkdir his
mv his__EleMVASecondNtupleProducer__SAMPLENAME__EVTSELECTION__0.root his/his__EleMVASecondNtupleProducer__SAMPLENAME__EVTSELECTION__TIMENOW__${skipEvents}__${lastEvent}.root

echo " "

if [[ -e "SECONDNTUPLEBASENAME__SAMPLENAME__EVTSELECTION__0.root" ]]
then
  echo Moving second ntuple output file: mv SECONDNTUPLEBASENAME__SAMPLENAME__EVTSELECTION__0.root ntu/SECONDNTUPLEBASENAME__SAMPLENAME__EVTSELECTION__TIMENOW__${skipEvents}__${lastEvent}.root
  mkdir ntu
  mv SECONDNTUPLEBASENAME__SAMPLENAME__EVTSELECTION__0.root ntu/SECONDNTUPLEBASENAME__SAMPLENAME__EVTSELECTION__TIMENOW__${skipEvents}__${lastEvent}.root
fi
