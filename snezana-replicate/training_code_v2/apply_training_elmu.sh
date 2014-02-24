#!/bin/sh

home=${PWD}
training_name=00
method1=MLP

lepton=elmu
analysis_channel=combined
btconf=_1btin

################## Running TMVAClassificationApplication
weights_file=weights/Training_${lepton}_${training_name} #/TMVAClassification_MLP.weights.xml

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
input_path=/mnt/xrootdb/alister/MVA_studies/samples # same as in run_training_elmu.sh
output_path=${home}/outputs # where you want your outputs to be written

# making MVAOutput folders
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
      for syst in jer jesd jesu jmr jmsd jmsu jvfsfd jvfsfu loose nominal AF2 qcd btagsfu btagsfd ctagsfu ctagsfd mistagsfu mistagsfd eltrigsfu eltrigsfd elrecosfu elrecosfd elidsfu elidsfd zjxsd zjxsu wjxsd wjxsu wjbbcccu wjbbcccd wjbbccu wjbbccd wjbb4u wjbb4d wjc4u wjc4d wjcau wjcad 
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
      for syst in jer jesd jesu jmr jmsd jmsu jvfsfd jvfsfu nominal AF2 qcd btagsfu btagsfd ctagsfu ctagsfd mistagsfu mistagsfd mutrigsfu mutrigsfd murecosfu murecosfd muidsfu muidsfd zjxsd zjxsu wjxsd wjxsu wjbbcccu wjbbcccd wjbbccu wjbbccd wjbb4u wjbb4d wjc4u wjc4d wjcau wjcad 
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
    for file in tprime_500_1M.root data.root qcd.root
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      ./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
  
  for syst in nominal jer jesd jesu jmr jmsd jmsu jvfsfd jvfsfu btagsfu btagsfd ctagsfu ctagsfd mistagsfu mistagsfd
    do
    for file in ttbar.root wjets.root zjets.root diboson.root singletop.root tprime_500.root tprime_550.root tprime_600.root tprime_650.root tprime_700.root tprime_750.root 
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      ./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
  
  for syst in AF2 zjxsu zjxsd wjxsu wjxsd wjbbcccu wjbbcccd wjbbccu wjbbccd wjbb4u wjbb4d wjc4u wjc4d wjcau wjcad 
    do
    for file in `ls -1 ${input_path}/${channel}/${syst}`
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      ./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
done 

for channel in el 
  do 
  for syst in eltrigsfu eltrigsfd elidsfu elidsfd elrecosfu elrecosfd
    do 
    for file in diboson.root singletop.root tprime_500.root tprime_550.root tprime_600.root tprime_650.root tprime_700.root tprime_750.root ttbar.root wjets.root zjets.root  
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
  for syst in muidsfu muidsfd murecosfu murecosfd mutrigsfu mutrigsfd 
    do
    for file in diboson.root singletop.root tprime_500.root tprime_550.root tprime_600.root tprime_650.root tprime_700.root tprime_750.root ttbar.root wjets.root zjets.root 
      do
      echo 'running over ' ${input_path}/${channel}/${syst}/${file}
      ./TMVAClassificationApplication ${input_path}/${channel}/${syst}/${file} ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file} ${weights_file} ${nevts}
      mv TMVApp.root ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/$${syst}/{file}
      echo 'moving output to ' ${output_path}/${channel}/${analysis_channel}/${btconf#*_}/${syst}/${file}
    done
  done
done
