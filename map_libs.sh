#!/usr/bin/env bash

USE_LIBS="
 SUSYTools GoodRunsLists\
 ApplyJetResolutionSmearing ApplyJetCalibration CalibrationDataInterface\
 egammaAnalysisUtils\
 egammaEvent JetResolution JetUncertainties\
 MissingETUtility MuonEfficiencyCorrections MuonMomentumCorrections\
 PileupReweighting ReweightUtils\
 TauCorrUncert\
 TrigMuonEfficiency ElectronEfficiencyCorrection PATCore TileTripReader"

USE_HEADER="SUSYTools MissingETUtility egammaAnalysisUtils GoodRunsLists\
 CalibrationDataInterface TauCorrUncert"


if [[ $1 == '-l' ]]
then
    for l in $USE_LIBS
    do
	echo -n -L$2/${l}/StandAlone -Wl,-rpath,$2/${l}/StandAlone -l${l}\ 
    done
elif [[ $1 == '-i' ]]
then
    for i in $USE_HEADER
    do 
	echo -n $2/$i\ 
    done
else 
    echo "ERROR: $0 needs flag -i or -l" >&2 
    false 
fi

