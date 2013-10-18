#!/bin/sh

home=${PWD}
training_name=07 #el
#training_name=165 #el
#training_name=00 #mu
method1=MLP
#method2=Cuts

input_path=/atlas/data1/userdata/snezana/TprimeAnalysis/2011_7TeV_full/TRCR-11-00-00-05/VLQTree
lepton=elmu
analysis_channel=combined
btconf=_1btin

################## Running TMVAClassificationApplication
weights_file=weights/Training_${lepton}_${training_name}/TMVAClassification_MLP.weights.xml

#nevts=-1 
###nevts=10000

echo '------------------------------'
echo 'Basic setup:'
echo 'input_path = ' ${input_path}
echo 'output_path = ' ${output_path}
echo 'weights_file = ' ${weights_file}
echo 'running over n events = ' ${nevts}
#echo '------------------------------'


# define the input and the output path
input_path=/atlas/data1/userdata/snezana/TprimeAnalysis/2011_7TeV_full/TRCR-11-00-00-05/VLQTree
output_path=/atlas/data1/userdata/snezana/TprimeAnalysis/2011_7TeV_full/TRCR-11-00-00-05/MVAOutput650_training_${lepton}_${training_name}

# making MVAOutput folders
cd /atlas/data1/userdata/snezana/TprimeAnalysis/2011_7TeV_full/TRCR-11-00-00-05
mkdir MVAOutput650_training_${lepton}_${training_name}
cd MVAOutput650_training_${lepton}_${training_name}
##### Make MVA_output folder
for channel in el mu
  do
  mkdir ${channel}
  cd ${channel}
  for jbin in combined
    do
    mkdir ${jbin}
    cd ${jbin}
    for btbin in 1btin
      do
      mkdir ${btbin}
    done
    cd .. #jbin
  done
  cd .. #channel
done

########### make subfolders
for channel in el
  do
  cd ${channel}
  for jbin in combined
    do
    cd ${jbin}
    for btbin in 1btin
      do
      cd ${btbin}
      for syst in jer jesd jesu jmr jmsd jmsu jvfsfd jvfsfu loose nominal AF2 qcd btagsfu btagsfd ctagsfu ctagsfd mistagsfu mistagsfd eltrigsfu eltrigsfd elrecosfu elrecosfd elidsfu elidsfd
	do
	mkdir ${syst}
      done
      cd .. #btbin
    done
    cd .. #jbin
  done
  cd .. #channel
done

for channel in mu
  do
  cd ${channel}
  for jbin in combined
    do
    cd ${jbin}
    for btbin in 1btin
      do
      cd ${btbin}
      for syst in jer jesd jesu jmr jmsd jmsu jvfsfd jvfsfu nominal AF2 qcd btagsfu btagsfd ctagsfu ctagsfd mistagsfu mistagsfd mutrigsfu mutrigsfd murecosfu murecosfd muidsfu muidsfd
	do
	mkdir ${syst}
      done
      cd .. #btbin
    done
    cd .. #jbin
  done
  cd .. #channel
done

#######################################################
# Applying NN

cd ${home}
for channel in el mu
  do
  
  for syst in nominal 
    do
    for file in data.root qcd.root tprime_500_lessPS.root tprime_500_morePS.root
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      #./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      #mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
  
  for syst in nominal jer jesd jesu jmr jmsd jmsu jvfsfd jvfsfu btagsfu btagsfd ctagsfu ctagsfd mistagsfu mistagsfd #qcd #nominal
    do
    for file in tprime_700.root tprime_750.root diboson.root singletop.root tprime_350.root tprime_400.root tprime_450.root tprime_500.root tprime_500_1M.root tprime_550.root tprime_600.root tprime_650.root ttbar.root wjets.root zjets.root
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      #./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      #mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
  
  for syst in AF2 
    do
    for file in `ls -1 ${input_path}/${channel}/${syst}`
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      #./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      #mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
done 

for channel in el 
  do 
  for syst in eltrigsfu eltrigsfd elidsfu elidsfd # elrecosfu elrecosfd
    do 
    for file in tprime_700.root tprime_750.root diboson.root singletop.root tprime_350.root tprime_400.root tprime_450.root tprime_500.root tprime_500_1M.root tprime_550.root tprime_600.root tprime_650.root ttbar.root wjets.root zjets.root #`ls ${input_path}/${channel}/${syst}` 
      do 
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      ./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done 
  done
done

for channel in mu
  do
  for syst in muidsfu muidsfd #murecosfu murecosfd mutrigsfu mutrigsfd 
    do
    for file in tprime_700.root tprime_750.root diboson.root singletop.root tprime_350.root tprime_400.root tprime_450.root tprime_500.root tprime_500_1M.root tprime_550.root tprime_600.root tprime_650.root ttbar.root wjets.root zjets.root #`ls ${input_path}/${channel}/${syst}`
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      #./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      #mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
done
