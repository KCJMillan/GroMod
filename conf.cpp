#include "conf.h"
#include "const.h"
#include "sort.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "timer-master/timer.h"
#include <vector>

using namespace std;


//FUNCTIONS

void  bondlengthAnalysis()

{


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

/*DEBUG*//*
cout<<"DEF OUT PATh :"<<opath<<endl<<"DEF IN PATH :"<<ipath<<endl<<"DEF CRIT RAD :"<<CritRadius<<endl<<"  DEF output name"<<oname<<endl<<"DEF NB ATOMS TO LOOK"<<nbAtomsToLook<<endl;
*/

//USER INPUTS

string gropath, temp, groname;
cout <<endl<<endl<< "       ** Bond Analysis in .gro file **"<<endl;
cout <<"Input the name of the .gro file (ex : conf.gro): "<<endl<<endl;
cout <<ipath;
getline(cin,temp);
gropath = ipath + temp ;

if(temp.length()>4){temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );
groname=temp;}
else {groname=temp;}

/*debug*/
/*cout<<".gro path : " << gropath<<endl;*/
cout<<endl;

//\\\\\\\\\\\\\\\\\\\\\\\\\\//
//   LOAD GRO FILE AND BOX  //
//\\\\\\\\\\\\\\\\\\\\\\\\\\//


//Count number of atoms in the system
int nbAtoms=0;
cout<<"Counting Atoms"<<endl;
LoadnbAtoms(gropath,&nbAtoms);
cout<<"There is  "<<nbAtoms<<" Atoms"<<endl;

//Make sure the .gro file is ok before starting
        if(nbAtoms!=0){


//Count number of frames in the system
cout<<"Counting Frames";cout<<endl;

int nbFrames=0;
LoadnbFrames(gropath, &nbFrames,nbAtoms);
float frames[nbFrames] , dtframe = 0 ;
cout<<"There is  "<<nbFrames<<" frames"<<endl;
cout<<"Saving Frames in memory"<<endl;
Loadframe(gropath,frames,nbAtoms,&dtframe);
cout<<"Total time available : "<<nbFrames*dtframe<<" ps"<<endl;

//Load box
//The size of the box represented by a,b,c and angles by alpha, beta, teta
float a,b,c;
float alpha, beta, gamma;
cout<<endl<<"Initialising Box"<<endl;
LoadBoxSizeAngles(gropath,boxpath,&a,&b,&c,&alpha,&beta,&gamma,nbAtoms);

/*DEBUG*//*
cout<<"a : "<<a<<"  b : "<<b<<"  c: "<<c<<endl<<"alpha : "<<alpha<<"  beta : "<<beta<<"  teta : "<<teta<<endl;
*/

////
//  Coordinates transformation init
////

//Calculate cosines and sines for the coordinate transform
float v=0; //v is an intermediate var to simplify the calculus
float m11=0,m12=0,m13=0,m22=0,m23=0,m33=0; //transform inverse matrix elements
float bm11=0,bm12=0,bm13=0,bm22=0,bm23=0,bm33=0; //transform matrix elements
float cos_alpha=0, sin_alpha=0, cos_beta=0,sin_beta=0, cos_gamma=0, sin_gamma=0;
float xtemp,ytemp,ztemp; //temporary var to replace the list of x,y,z



//calculate sin and cos
CalcTrigo(alpha,beta,gamma,&cos_alpha,&sin_alpha,&cos_beta,&sin_beta,&cos_gamma, &sin_gamma, &v);


//Calculare matrix elements
CalcMatrixElements(a,b,c,v,cos_alpha,sin_alpha,cos_beta,sin_beta,cos_gamma,sin_gamma,&bm11,&bm12,&bm13,&bm22,&bm23,&bm33);
//calculate inverse matrix elements
CalcInverseMatrixElements(a,b,c,v,cos_alpha,sin_alpha,cos_beta,sin_beta,cos_gamma,sin_gamma,&m11,&m12,&m13,&m22,&m23,&m33);

//ready to transform coordinates



//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////
//
//                                                                      BODY
//
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//\\



        //LOAD BUFFER

        int atomMax=0;
        LoadBuffers(&atomMax);

//Var needing nbAtoms (fixed size is more efficient than vector)
float x[atomMax],y[atomMax],z[atomMax],vx[atomMax],vy[atomMax],vz[atomMax]; string names[atomMax]; int indexAtom[atomMax];

//Creating the file and layout of the output file (neighbors data file)
string outputname=oname;
string path = opath+groname+outputname;

//Var needed for neighbor research
float bondLenght, closestLenght; int closestIndex, closestIndexTab[nbAtomsToLook] ;

//Vector since we want to modify items in this list
vector<float> lenghtsTab(atomMax,0);
vector<float> closestLenghtTab(nbAtomsToLook, 1500);

//VECTOR OF NEIGHBOR Rij=(xj-xi,yj-yi,zj-zi)

float Rij[3];

//CREATE OUTPUT file
ofstream file(path.c_str());

cout<<endl<<"Creating output : Bond data file"<<endl<<endl;

if(file)

{

    file<<"                                  # Bond analysis file #"<<endl;
    file<<nbAtoms<<" Atoms   "<<nbFrames<<"  Frames   "<<dtframe<<"  FrameStep"<<endl;

    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
    //                                                                                                                                              //
    //                                                      RESEARCH OF THE CLOSEST NEIGHBOR                                                        //
    //  2 Loops. For one atom you look at the distance of all the other atoms. Depending on the conf.conf file, we will seek n number of neighbors  //
    //                                                                                                                                              //
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//


        cout<<"--- Starting the Neighbor research ---"<<endl<<"Number of Neighbor research set to : "<<nbAtomsToLook<<endl<<endl;

    fstream filegro(gropath.c_str());

    cout<<"Loading number of atoms for each frame"<<endl<<endl;

            int nbAtomsAtFrame[nbFrames];
    cout<<"Start"<<endl;

          loadNbAtomEachFrame(gropath, nbAtomsAtFrame, nbFrames);

        //FRAME LOOP : BIG LOOP FOR EACH FRAME
    for( int l=0 ; l<nbFrames ; l++ )  {

    nbAtoms=nbAtomsAtFrame[l];

                //LAYOUT

    file<<"Atom number"<<"   Name"<<"   Nearest neighbor"<<"   Name"<<"    Bond Lenght(nm)";

                    for(int a=0;a<nbAtomsToLook-1;a++){ file<<"  "<<" Neighbor"<<a+2<<"      "<<"Name"<<"      "<<"Bond Lenght(nm)"; } file<<endl;

            //timer start
        tim.start();
            //load the coordinates for one frame
        LoadGroFile(filegro, gropath, nbAtoms,x,y,z,vx,vy,vz, names, indexAtom, frames[l] , dtframe );


            /*                                                              */
            //transform all thoses coordinates
      //  TransformAllCoordinates(nbAtoms,x,y,z,m11,m12,m13,m22,m23,m33);



            //writting the fame time we are looking at
        file<<"NbAtoms: "<< nbAtoms << " -- ";
        file<<"Time= "<< frames[l] << endl;


       //First LOOP : Fixed atom. The if and elses are for layout.

    for(int i=0;i<nbAtoms;i++){

        if(i<9){file<<"     "<<indexAtom[i]<<"         "<<names[i];}
        else if(i>8 && i<99){file<<"    "<<indexAtom[i]<<"         "<<names[i];}
        else if(i>98 && i<999){file<<"   "<<indexAtom[i]<<"         "<<names[i];}
        else if(i>998 && i<9999){file<<"  "<<indexAtom[i]<<"         "<<names[i];}

            //We have to put lenghts high to find always the first nearest atoms. If not, the value of a previous atom could be. We refresh Lenghts.
        closestLenght=999.0;
        bondLenght=9929.0;

        for(int k=0;k<nbAtomsToLook;k++){ closestLenghtTab[k]=1500; }

            //FIRST TIME OF SECOND LOOP : Just the research of the closest neighbor. (saved in the output file)

        for(int j=0;j<nbAtoms;j++)   {

            if(j!=i){

                        Rij[0]=x[j]-x[i];//cout<<"Before Rij[0]: "<<Rij[0]<<endl;
                        Rij[1]=y[j]-y[i];
                        Rij[2]=z[j]-z[i];

                        TransformCoordinates(&Rij[0],&Rij[1],&Rij[2],m11,m12,m13,m22,m23,m33);
                   //   cout<<"before pbc correction Rij[0]: "<<Rij[0]<<endl;


                      if(Rij[0]>0.5){ Rij[0] = Rij[0]-1.0 ; }
                      if(Rij[1]>0.5){ Rij[1] = Rij[1]-1.0 ; }
                      if(Rij[2]>0.5){ Rij[2] = Rij[2]-1.0 ; }

                      if(Rij[0]<-0.5){ Rij[0] = Rij[0]+1.0 ; }
                      if(Rij[1]<-0.5){ Rij[1] = Rij[1]+1.0 ; }
                      if(Rij[2]<-0.5){ Rij[2] = Rij[2]+1.0 ; }
              //  cout<<"After pbc correction Rij[0]: "<<Rij[0]<<endl;



                      TransformCoordinates(&Rij[0],&Rij[1],&Rij[2],bm11,bm12,bm13,bm22,bm23,bm33);

              //        cout<<"After Rij[0]: "<<Rij[0]<<endl; cin>>temp;

                      bondLenght = sqrt( Rij[0]*Rij[0]+Rij[1]*Rij[1]+Rij[2]*Rij[2] ) ;

                      lenghtsTab[j] = bondLenght ;


                    }


            if(j==i){ lenghtsTab[j]=999; }


                //The closest neighbor, its also in lenghtab[0], but for clarty, lets define the nearest like one var

            if( bondLenght < closestLenght ){ closestLenght = bondLenght; closestIndex=j+1;}


                                     }

             //LOOK WHICH *nbAtomsToLook* ARE AROUND AND CLOSEST

            for(int a=0;a<nbAtomsToLook;a++){

                for(int k=0;k<nbAtoms;k++){

                                        if( lenghtsTab[k] < closestLenghtTab[a] ) { closestLenghtTab[a] = lenghtsTab[k] ; closestIndexTab[a]=k+1; }

                                            //we remove this closest atom from the list, then the next loop will not count him again

                                        lenghtsTab[closestIndexTab[a]-1]=999;

                                            }

                                            }


            //WRITING THE OUTPUT FINAL with layout (if and elses)

        if(closestIndex<10){file<<"               "<<closestIndex;}
        else if(closestIndex>9 && closestIndex<99){file<<"              "<<closestIndex;}
        else if(closestIndex>98 && closestIndex<999){file<<"              "<<closestIndex;}
        else if(closestIndex>998 && closestIndex<9999){file<<"             "<<closestIndex;}

        file<<"         "<<names[closestIndex-1]<<"           "<<closestLenght;

                //Frame adapted depending the number of atoms to look
        for(int a=1;a<nbAtomsToLook;a++){


            if(a==1){
        if(closestIndexTab[a]<9){file<<"           "<<closestIndexTab[a];}
        else if(closestIndexTab[a]>8 && closestIndexTab[a]<99){file<<"          "<<closestIndexTab[a];}
        else if(closestIndexTab[a]>98 && closestIndexTab[a]<999){file<<"          "<<closestIndexTab[a];}
        else if(closestIndexTab[a]>998 && closestIndexTab[a]<9999){file<<"         "<<closestIndexTab[a];}

        file<<"       "<<names[closestIndexTab[a]-1]<<"           "<<closestLenghtTab[a];}

                    else{

        if(closestIndexTab[a]<9){file<<"           "<<closestIndexTab[a];}
        else if(closestIndexTab[a]>8 && closestIndexTab[a]<99){file<<"          "<<closestIndexTab[a];}
        else if(closestIndexTab[a]>98 && closestIndexTab[a]<999){file<<"          "<<closestIndexTab[a];}
        else if(closestIndexTab[a]>998 && closestIndexTab[a]<9999){file<<"         "<<closestIndexTab[a];}

        file<<"            "<<names[closestIndexTab[a]-1]<<"           "<<closestLenghtTab[a];

                    }


                                        }
        file<<endl;

                                }


        //GET TIME USED OF ONE FRAME LOOP

        finaltime=tim.getElapsedTime(MILLISEC)/1000.0;

        tim.stop();

        //We calculate an average of the time to do one loop and then multiply by the loops remaining
        mean  += finaltime  ;

            chronot= (mean/(float)l)*(nbFrames-l);


        //Refresh the chrono every x time for long times and saving efficiency



         if( (l) % 10 == 0) {

            //Layout for longer times
         if(chronot>60 && chronot<3600){

                        minutes=chronot/60; hours = minutes / 60; temp5=(int)chronot%60;

                        cout<<"\rReading Time: "<< frames[l] <<" ------ Time Remaning :"<<(int)minutes<<" minutes and "<<temp5<<" secs                "<<flush;

                                       }

        else if(chronot>3600 ){

                        minutes=chronot/60; hours = minutes / 60; temp5=(int)chronot%60;

                        cout<<"\rReading Time: "<< frames[l] <<" ps ------ Time Remaning :"<<(int)hours<<" hours and "<<temp5<<" minutes                         "<<flush<<flush;

                                       }

        else  {cout<<"\rReading Time: "<< frames[l] <<" ps ------ Time Remaning :"<<chronot<<" secs "<<flush;}




                           }


}

cout<<endl<<endl<<"Output created"<<endl<<endl;
cout<<"Neighbor research completed"<<endl;

cout<<">>";

}




else{cout<<"Error while creating the processed file"<<endl;}

}

else{cout<<endl<<"Neighbor Research aborted.."<<endl<<endl<<"************";}


}






void ionPairAnalysis(){

    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
    //                                                                                                                                              //
    //                                                      ANALYSIS OF THE neighbor data FILE                                                      //
    //                                                                 Ion pairing                                                                  //
    //                                                                                                                                              //
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//

cout<<endl<<" ------- Ion Pairing Analysis --------"<<endl;




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
LoadNeighborFrameFromDataFile(DataFilePath, nbFrames, frames);

  /*cout<<"first frame ; "<<frames[0]<<" second : "<<frames[1]<<" Last ; "<<frames[nbFrames-1];*/

cout<<endl<<endl;

string atomName[atomMax];
int indexAtom[atomMax];

//loading nb atoms each frame

      //   int nbAtomsAtFrame[nbFrames];

        //  loadNbAtomEachFrame(DataFilePath, nbAtomsAtFrame, nbFrames);


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
             //   nbAtoms=nbAtomsAtFrame[l];

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

        finaltime=tim.getElapsedTime(MILLISEC)/1000.0;

        tim.stop();

        //We calculate an average of the time to do one loop and then multiply by the loops remaining
        mean  += finaltime  ;


        //Refresh the chrono every x time for long times and saving efficiency

         if( (l) % 10 == 0) {

         chronot= (mean/(float)l)*(nbFrames-l);


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
















    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
    //                                                                                                                                              //
    //                                                           SORT   .GRO FILE                                                                   //
    //                                                                 Ion                                                                          //
    //                                                                                                                                              //
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//


void sortConf(){


cout<<endl<<" ---------- Sort .GRO File -----------"<<endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////                                                                              /////////////////////////////////////////
////////                                 HEAD                                         ///////////////////////////////////////////
////////                  VAR DEFINITIONS + USER INPUTS                               //////////////////////////////////////
////////                                                                              //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//LOAD CONFIGS

string gropath,opath,ipath,oname,boxpath; int nbAtomsToLook;
cout<<"      --------Type help for help----------"<<endl<<endl;
cout<<"";
LoadConfig(&opath,&ipath, &nbAtomsToLook,&oname,&boxpath);


//USER INPUTS
cout <<endl<<"Path of Configuration File .gro to modify : "<<endl<<endl;
cout<<ipath;
string temp, groname;
getline(cin,temp);

gropath=ipath+temp;

if(temp.length()>4){temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );
groname=temp;}
else{groname="default";}
cout<<endl;

    //Loading constant vars
int nbAtoms=0;
cout<<"Counting Atoms"<<endl;
LoadnbAtoms(gropath,&nbAtoms);
cout<<"There is  "<<nbAtoms<<" Atoms"<<endl;

//Count number of frames in the system
cout<<"Counting Frames"<<endl;
int nbFrames=0;
LoadnbFrames(gropath, &nbFrames,nbAtoms);
float frames[nbFrames] , dtframe = 0 ;
cout<<"There is  "<<nbFrames<<" Frames"<<endl;
Loadframe(gropath,frames,nbAtoms,&dtframe);
cout<<"Time available "<<nbFrames*dtframe<<" ps"<<endl;





//////////////////
/////////////
////////        MENU
//////////////
///////////////////

    if(nbAtoms!=0){


//COMMANDS

string command;
cout <<endl<<"You want to : (remove, replace, add, open,exit ... ?)"<<endl<<endl;
bool maon=true;

//init conf for programs
string prg1,prg2;

LoadProgramConfig(&prg1,&prg2);

prg1=prg1+" ";
prg2=prg2+" ";

while (maon)

{

    getline(cin,command);

    if(command=="remove"){sortRemove(gropath,groname,opath,nbAtoms,nbFrames,frames,dtframe);}

    else if(command=="change gro path"){cout<<"**Be careful, change only gro files with same amount of frames than the first one loaded. If you want to load a new gro file with differents number of frames just type: exit and sort"<<endl;
                                        changeGroFilePath(&gropath, &groname,ipath); loadNewGroConstants(&nbAtoms,&nbFrames,&dtframe,frames,gropath);}

    else if (command== "open original"){string program = prg1 + gropath + " &"; system(program.c_str()); }
    else if (command == "open sorted"){string program = prg1 + opath + groname + "Processed.gro" + " &"; system(program.c_str()); }

    else if (command== "exit"){maon=false;}


}



}



else{cout<<"Error loading .gro File "<<endl;}








}














