#! /bin/bash

cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}

export RECREL=CMSSW_10_6_8
export SCRAM_ARCH=slc7_amd64_gcc700

export InFile=$2
export OutFile="DoubleEG_16H_${1}.root"

echo $InFile

echo $OutFile

echo "--> Environment"
date
hostname
uname -a
#df -kl
#limit coredumpsize 0

source /cvmfs/cms.cern.ch/cmsset_default.sh
echo "-> which edg-gridftp-ls"
which edg-gridftp-ls
echo "-> which globus-url-copy"
which globus-url-copy
#echo "-> which srmcp"
#which srmcp

pwd
echo "--> End of env testing"


# ----------------------------------------------------------------------
# -- Setup CMSSW
# ----------------------------------------------------------------------
echo "--> Setup CMSSW"
pwd
date
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 project CMSSW CMSSW_10_6_8`
cd CMSSW_10_6_8/src
#scp /afs/cern.ch/work/p/prsaha/public/service_task/egamma/CMSSW_10_6_8/src/Setup.tar.gz .

tar -xvf ../../Setup.tar.gz

eval `scramv1 runtime -sh`

cp XGBoostCMSSW/XGBoostInterface/toolbox/*xml $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/

scram setup rabit

scram setup xgboost

scramv1 b

cd EgammaIdCorrection/DifferentialEgammaIdCorrector/test/

cp ../../../../../DifferentialEgammaIdCorrector_cfg.py .

cmsRun DifferentialEgammaIdCorrector_cfg.py

ls -rtl
xrdcp DoubleEG_16H_${1}.root root://eosuser.cern.ch///eos/user/p/prsaha/ServiceTask/egamma_output/Data_16/



echo "run: This is the end, my friend"
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_6_8
