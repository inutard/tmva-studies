#!/bin/sh

training_name=00
training_type=b
method1=MLP

input_path=/mnt/xrootdb/alister/MVA_studies/samples # add your input path here
lepton=elmu
analysis_channel=combined
btconf=_1btin

################## Running TMVClassification

cd weights
mkdir -p Training_${lepton}_${training_name}${training_type}
cd Training_${lepton}_${training_name}${training_type}
mkdir -p plots
cd .. #Training_${lepton}_${training_name}
cd .. #weights

./TMVAClassification ${method1} -p ${input_path} ${lepton} ${analysis_channel} ${btconf} ${training_type} Training_${lepton}_${training_name}${training_type} > training.log
echo './TMVAClassification ${method1} -p ${input_path} ${lepton} ${analysis_channel} ${btconf} > training.log'

#tail -n 25 training.log > training_summary.log
cd weights
cp ../TMVAClassification* Training_${lepton}_${training_name}${training_type}
mv ../training.log Training_${lepton}_${training_name}${training_type}
cp ../TMVA.root Training_${lepton}_${training_name}${training_type}
cd .. #weights

