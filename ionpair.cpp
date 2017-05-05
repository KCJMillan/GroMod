#include "ionpair.h"
#include "const.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "timer-master/timer.h"
#include <vector>

using namespace std;

void test3(){

    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
    //                                                                                                                                              //
    //                                                      ANALYSIS OF THE neighbor data FILE                                                      //
    //                                                                 Ion pairing                                                                  //
    //                                                                                                                                              //
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//

cout<<endl<<" ------- Ion Pairing Analysis --------"<<endl;

int tempfll;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////                                                                              /////////////////////////////////////////
////////                                 HEAD                                         ///////////////////////////////////////////
////////                  VAR DEFINITIONS + USER INPUTS                               //////////////////////////////////////
////////                                                                              //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Timer

Timer tim;
float finaltime=0,mean=0;
float chronot=0;
int temp5;
int hours;
int minutes;

//LOAD CONFIGS

string pairOutputName, opath, temppath, temppath212,ipath,oname,boxpath; double CritRadius; int nbAtomsToLook;

LoadConfig(&temppath,&temppath212, &nbAtomsToLook,&oname,&boxpath);
LoadPairConfig(&opath, &ipath, &CritRadius,&pairOutputName);

cout<<"Critical Radius set to : "<<CritRadius<<" nm"<<endl;
    /*cout<<" opath : "<<opath<<endl<<" CritRadius : "<<CritRadius<<endl;*/
        cout<<endl;
//Vars




//USER INPUTS

string bondpath, temp, bondname;
cout <<"Input the name of the Neighbor .eng file : "<<endl<<endl;
cout <<ipath;
getline(cin,temp);
bondpath = ipath + temp ;

if(temp.length()>4){temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );
bondname=temp;}

else { bondname="default";}

/*debug*/
/*cout<<".gro path : " << gropath<<endl;*/
cout<<endl;




    //constant vars
int nbAtoms, nbFrames=0;
float dtFrames;
string DataFilePath=bondpath;

    //non constant vars
float frameToLoad=1000;
int indexAtomToLookAt=1;


        //LOAD BUFFER

        int atomMax=0;
        LoadBuffers(&atomMax);

////////
////////
//Asking user for ions names
////////
////////


int nbDifferentIons=2;

string userAtomName[nbDifferentIons]; string userTemp;


for ( int p=0; p < nbDifferentIons ; p++ ) { cout<<"Name of the specie "<< p+1 <<" of ions : "; cin>>userTemp ; userAtomName[p]=userTemp; }

cout<<endl;

    /*cout << " First specie of ions : "<< userAtomName[0] <<endl;*/

//////////////
//////////////
//Loading constants
//////////////
//////////////

cout<<"Loading parameters"<<endl;
LoadNeighborDataFileConstants(DataFilePath, &nbAtoms,&nbFrames,&dtFrames);
cout<<"nbAtoms : "<<nbAtoms<<"  nbFrames : "<<nbFrames<<"  dtFrames: "<<dtFrames<<endl;
cout<<" Total time available : "<<nbFrames*dtFrames<<" ps"<<endl;

    if(nbFrames!=0){

float frames[nbFrames];
cout<<"Saving Frames in memory "<<endl;
LoadNeighborFrameFromDataFile(DataFilePath, nbFrames, frames);

  /*cout<<"first frame ; "<<frames[0]<<" second : "<<frames[1]<<" Last ; "<<frames[nbFrames-1];*/

string atomName[atomMax];
int indexAtom[atomMax];

//loading nb atoms each frame

int nbAtomsAtFrame[nbFrames];
cout<<"Saving Atoms per frame in memory "<<endl;
LoadAtomPerFrames(DataFilePath, nbAtomsAtFrame, nbFrames);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////                                                                              /////////////////////////////////////////
////////                                 BODY                                         ///////////////////////////////////////////
////////                 RESEARCH THE NUMBER OF ION PAIRING                           //////////////////////////////////////
////////                                                                              //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//CREATE IONPAIR OUTPUT FILE + loop vars

string pairPath = opath + bondname + pairOutputName;
string tru="no";
ofstream pairFile(pairPath.c_str());


int l=0, k,a;
int op=0;
float nbGlobalPairs=0;
float nbGlobalBondLenght=0;
float pairPerFrame[nbFrames];
int nbPairs=0;



float ToNeighborLenght[nbAtomsToLook*atomMax];
string ToNeighborName[nbAtomsToLook*atomMax];
int ToNeighborIndex[nbAtomsToLook*atomMax];

vector<float> pairingWith;
vector<string> pairingName;


if(pairFile ){    cout<<"---Starting Ion Pair Research---"<<endl;


                pairFile<<"                                  # Ion Pair Data File #"<<endl;
                pairFile<<nbAtoms<<" Atoms   "<<nbFrames<<"  Frames   "<<dtFrames<<"  FrameStep "<<nbDifferentIons<<"  Different Ions"<<endl;


                fstream DataFilep( DataFilePath.c_str() );

                //Loading file for each frames

    for(l=0;l<nbFrames;l++){

                pairFile<<"Ion      " <<"Index    "<<"Pairing ?          "<<" How much          ";
                for(int n=0;n<nbAtomsToLook;n++){pairFile<<" With index"<<"       Name"<<"     ";}
                pairFile<<endl;
                pairFile<<"Time= "<< frames[l] << endl;

                tim.start();

                nbAtoms=nbAtomsAtFrame[l];

            LoadNeighborDataFile(DataFilep, DataFilePath, nbAtoms , atomName, indexAtom, ToNeighborIndex, ToNeighborLenght,ToNeighborName,nbAtomsToLook,dtFrames,frames[l] );

        //LOOP ; FOR EACH ATOM WE LOOK THE NEIGHBOORHOOD AND SEE IF THERE IS PAIRS


    for(k=0;k<nbAtoms;k++){

   // LoadNeighborDataFile(DataFile3, DataFilePath, k , nbAtoms ,&atomName,&indexAtom, ToNeighborIndex, ToNeighborLenght,ToNeighborName,nbAtomsToLook,dtFrames,frames[l] );


                //CAN DO A LOOK ON SPECIE A AND LOOK ALL SPECIE B TO NOT HAVE A FIXED SIZE ON IONS TO LOOK>MAYBYE ILL DO THAT AFTER

                //IF SPECIE 0 DO BOND WITH SPECIE 1

                if(atomName[k]==userAtomName[0]){


                        for(a=(0+k*nbAtomsToLook);(a<nbAtomsToLook+nbAtomsToLook*k);a++){


                            if(ToNeighborName[a]==userAtomName[1] && ToNeighborLenght[a]<=CritRadius )

                                                        {
                                                nbPairs++; nbGlobalPairs++; nbGlobalBondLenght+=ToNeighborLenght[a]; pairingWith.push_back(ToNeighborIndex[a]); pairingName.push_back(ToNeighborName[a]);
                                                        }

                                                    }

                            if(nbPairs!=0){tru="yes" ;}

                                            }

                //IF SPECIE 1 DO BOND WITH SPECIE 0
                else if(atomName[k]==userAtomName[1]){


                        for(a=(0+k*nbAtomsToLook);(a<nbAtomsToLook+nbAtomsToLook*k);a++){

                            if(ToNeighborName[a]==userAtomName[0] && ToNeighborLenght[a]<=CritRadius )

                                                        {
                                                nbPairs++;nbGlobalPairs++; nbGlobalBondLenght+=ToNeighborLenght[a]; pairingWith.push_back(ToNeighborIndex[a]); pairingName.push_back(ToNeighborName[a]);
                                                        }

                                                    }

                        if(nbPairs!=0){tru="yes" ; }

                                                    }



                //WRITTING IN FILE WITH LAYOUT
                if(atomName[k]==userAtomName[0]){op=1;}
                if(atomName[k]==userAtomName[1]){op=1;}

                if(tru=="yes" ){


                       if(indexAtom[k]<10){

                                    pairFile<<atomName[k]<<"         "<<indexAtom[k]<<"        "<<tru<<"                       "<<nbPairs<<"                     ";
                                    for(a=0;a<nbPairs;a++){pairFile<<pairingWith[a]<<"        "<<pairingName[a]<<"          "; }

                                        }

                       else             {
                                    pairFile<<atomName[k]<<"         "<<indexAtom[k]<<"       "<<tru<<"                       "<<nbPairs<<"                     ";
                                    for(a=0;a<nbPairs;a++){pairFile<<pairingWith[a]<<"        "<<pairingName[a]<<"          "; }

                                        }

                             }

                if(tru=="no" ){

                    pairFile<<atomName[k]<<"         "<<indexAtom[k]<<"         "<<tru<<"                       "<<nbPairs<<"         ";

                            }


                //RESET THE VAR USED

                nbPairs=0;
                pairingName.clear();
                pairingWith.clear();
                tru="no";
                op=0;


                pairFile<<endl;


                        } //end loop on each atom

              //GET TIME USED OF ONE FRAME LOOP

        finaltime= tim.getElapsedTime(MILLISEC)/1000.0;

        tim.stop();

        //We calculate an average of the time to do one loop and then multiply by the loops remaining
        mean  += finaltime  ;


        //Refresh the chrono every x time for long times and saving efficiency

         if( (l) % 10 == 0) {

         chronot= (mean/(float)l)*(float)(nbFrames-l);


            //Layout for longer times
         if(chronot>60 && chronot<3600){

                        minutes=chronot/60; hours = minutes / 60; temp5=(int)chronot%60;

                        cout<<"\rReading Time: "<< frames[l] <<" ------ Time Remaning :"<<(int)minutes<<" minutes and "<<temp5<<" secs "<<flush;

                                       }

          else  {cout<<"\rReading Time: "<< frames[l] <<" ------ Time Remaning :"<<chronot<<" secs "<<flush;}



                           }



                   pairPerFrame[l]=nbGlobalPairs/(float)l;


                    }//End loop over frames

        cout<<"Completed"<<endl;
                nbGlobalBondLenght=nbGlobalBondLenght/nbGlobalPairs;
                nbGlobalPairs=nbGlobalPairs/nbFrames;
                cout<<endl;

        cout<<endl<<"Writting Analysis file and xvg file"<<endl;


}


else {cout<<"Error creating the output file"<<endl ; }




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////                                                                              /////////////////////////////////////////
////////                             ANALYSIS                                         ///////////////////////////////////////////
////////                 RESEARCH THE NUMBER OF ION PAIRING                           //////////////////////////////////////
////////                                                                              //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


string ApairPath = opath + bondname + "IonAnalysis.eng";
string ApairVStime = opath + bondname + "IonPairTime.xvg";
ofstream ApairFile(ApairPath.c_str());
ofstream ApairVStimeFile(ApairVStime.c_str());

//VARS
float averageNbPairs;


if(ApairFile){

        ApairFile<<"Analysis of the ion pairs in the system"<<endl<<endl;

        ApairFile<<"Average Number of Pairs over "<<nbFrames*dtFrames<<" ps : "<<nbGlobalPairs<<" pairs"<<endl;
        ApairFile<<"Average Bond Pair Lenght over "<<nbFrames*dtFrames<<" ps : "<<nbGlobalBondLenght<<" nm" <<" / "<<nbGlobalBondLenght*10<<" A"<<endl;
        ApairFile<<"Average Number of double pairs over "<<nbFrames*dtFrames<<" Frames : "<<" double pairs" <<endl;

}
else{cout<<"Error creating Analysis ion pair output file"<<endl;}

if(ApairVStimeFile)
{
    ApairVStimeFile<<"nb of Pairs"<<"    "<<"  Frame"<<endl;
    for(int l=0;l<nbFrames;l++){ApairVStimeFile<<frames[l]<<"    "<<pairPerFrame[l]<<endl;}

}

else{cout<<"Error creatin Analysios ion pair as function of time file"<<endl;}







        cout<<"Completed"<<endl<<endl;

        }


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void LoadAtomPerFrames(std::string DataFilePath,int nbAtomsAtFrame[],int nbFrames){

    ifstream file(DataFilePath.c_str());
    char temp='j';
    int tempatom=0;
    float t=0;
    string temp2;

    if(file){

    for(int j=0;j<nbFrames;j++){      if(j%15==0)cout<<"\r frame -> "<< j <<"      "<<flush;

        while(temp != ':'){ file.get(temp);  }

        file>>tempatom;
        nbAtomsAtFrame[j]=tempatom;

        file.get(temp);

                            }    cout<<"\r"<<flush;
        }

else{cout<<"ERROR "<<endl;}


}

