#!/bin/sh

training_name=11 #el
method1=MLP
#method2=Cuts

input_path=/atlas/data1/userdata/snezana/TprimeAnalysis/2011_7TeV_full/TRCR-11-00-00-05/VLQTree_my
lepton=elmu
analysis_channel=combined
btconf=_1btin

################## Running TMVClassification

cd weights
mkdir Training_${lepton}_${training_name}
cd Training_${lepton}_${training_name}
mkdir plots
cd .. #Training_${lepton}_${training_name}
cd .. #weights

./TMVAClassification ${method1} -p ${input_path} ${lepton} ${analysis_channel} ${btconf} Training_${lepton}_${training_name} > training.log
echo './TMVAClassification ${method1} -p ${input_path} ${lepton} ${analysis_channel} ${btconf} > training.log'

#tail -n 25 training.log > training_summary.log
cd weights
mv TMVAClassification* Training_${lepton}_${training_name}
mv ../training.log Training_${lepton}_${training_name}
cp ../TMVA.root Training_${lepton}_${training_name}
cd .. #weights

