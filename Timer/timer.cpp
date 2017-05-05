#include "timer.h"

Timer::Timer(){
state = false;
}

Timer::~Timer(){}

//Time point of the start
void Timer::start(){
if(state == false){
	state = true;
	startPoint = system_clock::now();
}

else{
	//cerr<<"Error: Timer not stopped"<<endl;
}
}

//Time point of the end of execution
void Timer::stop(){
if(state==true){
	endPoint = system_clock::now();
	state = false;
}

else{
	//cerr<<"Error: Timer never started"<<endl;
}
}

//Return elapsed time in given unit
int Timer::getElapsedTime(timeUnit unit){
if(state == false){
	if(unit == SEC)
		return (int)duration_cast<seconds>(endPoint-startPoint).count();

	else if(unit == MILLISEC)
		return (int)duration_cast<milliseconds>(endPoint-startPoint).count();

	else if(unit == MICROSEC)
		return (int)duration_cast<microseconds>(endPoint-startPoint).count();
}

else{
	//cerr<<"Error: Timer still ticking";
	if(unit == SEC)
		return (int)duration_cast<seconds>(system_clock::now()-startPoint).count();
	else if(unit == MILLISEC)
		return (int)duration_cast<milliseconds>(system_clock::now()-startPoint).count();
	else if(unit == MICROSEC)
		return (int)duration_cast<microseconds>(system_clock::now()-startPoint).count();
}

}

//Print the returned finished date and time using ctime(&end_time)
time_t Timer::getFinishedTime(){
if(state == false){
	time_t end_time = system_clock::to_time_t(endPoint);
	return end_time;
}
else
	return (time_t)NULL;
}
