/******************************************************
*Class for calculating time with microsecond precision*
******************************************************/

#include <iostream>
#include <chrono>
#include <ctime>

using namespace std;
//New namespace set in C++11
using namespace chrono;

//Units of time
enum timeUnit{
SEC,
MILLISEC,
MICROSEC
};

class Timer{
public:
Timer();

~Timer();

void start();
void stop();

int getElapsedTime(timeUnit unit);
time_t getFinishedTime();

private:
//Time point private variables and a flag for timer state
bool state;
time_point<system_clock> startPoint, endPoint;

};
