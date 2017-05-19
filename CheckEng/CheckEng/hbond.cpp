#include "hbond.h"
#include <iostream>
#include <string>
#include <fstream>
#include "const.h"
#include "timer-master/timer.h"

using namespace std;

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
////                                                /////
////            MAIN FUNCTION                       /////
////                                                /////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


void hbond(){

//Timer
Timer tim;
float finaltime=0,mean=0;
float chronot=0;
int temp5;
int hours;
int minutes;

//LOAD CONFIG FROM CONFIG FILE

string opath,ipath,oname,boxpath; int nbAtomsToLook;
LoadConfig(&opath,&ipath, &nbAtomsToLook,&oname,&boxpath);

float CritRadius=0, CritAngle=0;string oxy_name,h1_name,h2_name;
LoadHbondConfig(&opath, &ipath, &CritRadius, &oname, &CritAngle,&oxy_name,&h1_name,&h2_name);

int atomMax=0;
LoadBuffers(&atomMax);


/*//DEBUG
cout<<"DEF OUT PATh :"<<opath<<endl<<"DEF IN PATH :"<<ipath<<endl<<"DEF CRIT RAD :"<<CritRadius<<endl<<"  DEF output name"<<oname<<endl<<"DEF NB ATOMS TO LOOK"<<nbAtomsToLook<<endl;
*/


//USER INPUTS

string NeighborPath, temp;
cout <<endl<<endl<< "       *** HBOND RESEARCH ***"<<endl;
cout <<"Input the path of the neighbor data file : "<<endl<<endl;
cout <<ipath;
getline(cin,temp);
NeighborPath = ipath + temp ;

if(temp.length()>4){temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );
oname=temp+oname;}

else {oname=temp;}

//VAR DECLARATION

    //constant vars for one frame
int nbAtoms=0, nbFrames;
float dtFrames;

    //non constant vars
float ToNeighborLenght[nbAtomsToLook];
string ToNeighborName[nbAtomsToLook];
int indexAtom; int ToNeighborIndex[nbAtomsToLook];
string atomName;
float frameToLoad=1000;
int indexAtomToLookAt=1;

//////////////
//////////////
//Loading constants
//////////////
//////////////


cout<<"Loading parameters"<<endl;
LoadNeighborDataFileConstants(NeighborPath, &nbAtoms,&nbFrames,&dtFrames);


        if(nbAtoms!=0){



cout<<"nbAtoms : "<<nbAtoms<<"  nbFrames : "<<nbFrames<<"  dtFrames: "<<dtFrames<<endl;
cout<<" Total time available : "<<nbFrames*dtFrames<<" ps"<<endl;

float frames[nbFrames];
LoadNeighborFrameFromDataFile(NeighborPath, nbFrames, frames);

    // cout<<"first frame ; "<<frames[0]<<" second : "<<frames[1]<<" Last ; "<<frames[nbFrames-1];

cout<<endl<<endl;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////                                                                              /////////////////////////////////////////
////////                                 BODY                                         ///////////////////////////////////////////
////////                        RESEARCH THE NUMBER OF HBONDS                         //////////////////////////////////////
////////                                                                              //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


string hbondPath = opath + oname ;
cout<<hbondPath<<endl;
string tru="no";

ofstream hFile(hbondPath.c_str());


if(hFile){


    cout<<"---Starting Hydrogen Bond Research---"<<endl;

                hFile<<"                                  # HBond Data File #"<<endl;
                hFile<<nbAtoms<<" Atoms   "<<nbFrames<<"  Frames   "<<dtFrames<<"  FrameStep "<<nbFrames*dtFrames<<"  ps"<<endl;








}






else{cout<<"Error while creating hbond data file"<<endl;}



}

else{cout<<"Could not open Neighbor Data File"<<endl;}

}
