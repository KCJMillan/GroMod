#uncomment if in clusters

#module load gcc

mkdir ../Executable

mkdir ../Executable/ConfigFiles

cp /CheckEng/boxConf.box ../

gcc ../CheckEng/*.cpp ../CheckEng/timer-master/timer.cpp std=c++11

