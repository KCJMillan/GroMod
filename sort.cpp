
#include <string>
#include "const.h"
#include <iostream>
#include <vector>
#include "sort.h"
#include <fstream>
#include "timer-master/timer.h"




using namespace std;

/////////////////////////
////            MENU FUNCTIONS
/////////////////////////


void sortRemove(string gropath,string groname, string opath, int nbAtoms,int nbFrames,float frames[], float dtframe)
{

    string temp; int nbToRemove=1;
    string NameToRemove;
    bool cont=true;



                    cout<<"Which one ? (name): "<<endl;
                    cin>>NameToRemove;


    cout<<endl<<"Criteria of removal ? (all, z< , z>, ...) : "<<endl<<"type help to see all criterias"<<endl<<endl;

    while(cont){

    cin>>temp;
    if(temp == "all"){removeAll(gropath,opath, groname, nbAtoms,nbFrames,frames,dtframe,nbToRemove,NameToRemove);  cont=false;   }
    else if (temp=="z<"){removeZinf(gropath,opath, groname, nbAtoms,nbFrames,frames,dtframe,nbToRemove,NameToRemove); cont=false; }

        }



}






////////////////////////////////////////
//////////////////////////
////                           SORT FUNCTIONS
/////////////////////////
//////////////////////////////////////

//REMOVE ALL
void removeAll(string gropath,string opath,string groname, int nbAtoms,int nbFrames,float frames[], float dtframe, int nbToRemove, string NameToRemove)
{

        string ProcessedOpath=opath+groname+"_Processed.gro";
        string lopath=opath+"Temp.gro";
        string temp;


        //LOAD BUFFER

        int atomMax=0;
        LoadBuffers(&atomMax);


        //CREATE First file : We remove all the atoms selected
            //Create and then close to use fstream
        ofstream createTemp(lopath.c_str());
        createTemp.close();

        ofstream createTemp2(ProcessedOpath.c_str());
        createTemp2.close();


        fstream file(lopath.c_str());
        fstream file2(ProcessedOpath.c_str());

        int cnt=0;
        string lineIn[nbAtoms+1];
        float x[nbAtoms],y[nbAtoms],z[nbAtoms],vx[nbAtoms],vy[nbAtoms],vz[nbAtoms];        string resNumbName[nbAtoms]; string names[nbAtoms]; int indexAtom[nbAtoms];

        if(file){

            fstream GroFile(gropath.c_str());


            ///////LOADING SECTION

          string firstline=getFirstLineFromGro(gropath);

          int bufferSize=(atomMax+1)*nbFrames;

          cout<<"Loading number of atoms for each frame"<<endl;

            int nbAtomsAtFrame[nbFrames];

          loadNbAtomEachFrame(gropath, nbAtomsAtFrame, nbFrames);

          cout<<"Loading all lines in gro file"<<endl;

            vector<string> lineIn(bufferSize);
          LoadGroLines(gropath, lineIn, nbAtoms,nbFrames);


          cout<<"Loading coordinates"<<endl;

           vector<string> resNumbName(bufferSize); vector<string> names(bufferSize); vector<int> indexAtom(bufferSize);

           vector<float> z(bufferSize);

           LoadGro(gropath,resNumbName, nbAtoms,z, names, indexAtom, nbFrames );

           int it=0;int whereItWas=0;  int tempp;


           for(int l=0;l<nbFrames;l++){

                        nbAtoms=nbAtomsAtFrame[l];

                file<<firstline; file<<"   "<<frames[l]<<endl<<" "<<nbAtoms<<endl;


                for(int p= whereItWas ; p<whereItWas+nbAtoms ;p++){

                if(names[p]!=NameToRemove){ file<<lineIn[p]<<endl; }

                else if(names[p]==NameToRemove && (l == 0) ){cout<<names[p]<<" "<<indexAtom[p]<<" removed"<<endl;cnt++;}

                                  tempp=p;     }//end loop on one frame

                    whereItWas=tempp+2;

                file<<lineIn[whereItWas-1]<<endl;

                if(l==0){cout<<endl<<cnt<<" "<<NameToRemove<<" removed"<<endl;
                         cout<<nbAtoms<<" Entities remaining"<<endl<<endl;

                        }

                if(l==1)cout<<"Processing for each Frame -- Please Wait"<<endl;                   }//end loop on each frame



            nbAtoms=nbAtoms-cnt;



                        cout<<" "<<flush;
                        cout<<"\rSaving File - please wait"<<endl<<endl;

             file.seekp(0, ios::beg);


            for(int l=0;l<nbFrames;l++){

                    getline(file,temp);
                    file2<<temp<<endl;
                    file2<<"  "<<nbAtoms<<endl;
                    getline(file,temp);

                for(int m=0;m<nbAtoms+1;m++){

                    getline(file,temp);
                    file2<<temp<<endl;

                                          }

                                }

            cout<<endl<<"Completed"<<endl;

                if( remove( lopath.c_str() ) != 0 )
                cout<< "Error deleting temporary file"<<endl;
                    else{}


                }

        else{cout<<"Could not create the output file in function removeAll()."<<endl<<" We tried to create a file at this path : "<<lopath<<endl;}


}


//REMOVE Z<


void removeZinf(string gropath,string opath,string groname, int nbAtoms,int nbFrames,float frames[], float dtframe, int nbToRemove, string NameToRemove)

{


//Timer
Timer tim;
float finaltime=0,mean=0;
float chronot=0;
int temp5;
int hours;
int minutes;


        string ProcessedOpath=opath+"_Processed.gro";
        string lopath=opath+"Temp.gro";
        string temp;
        float Criteria=0;

        //LOAD BUFFER

        int atomMax=0;
        LoadBuffers(&atomMax);

        cout<<"z< to ? (nm) "<<endl;
        cin>>Criteria;

        //CREATE First file : We remove all the atoms selected
            //Create and then close to use fstream
        ofstream createTemp(lopath.c_str());
        createTemp.close();

        ofstream createTemp2(ProcessedOpath.c_str());
        createTemp2.close();

        ////////////////////
        //Open output files
        fstream file(lopath.c_str());
        fstream file2(ProcessedOpath.c_str());

        int cnt=0,cnt1=0;int nbAtomsEachFrame[nbFrames];

        string pk;



        cout<<"Creating Buffer Array"<<endl;


        if(file){

        int nbAtomsInit=nbAtoms;


          ///////LOADING SECTION
          string firstline=getFirstLineFromGro(gropath);

          int bufferSize=(atomMax+1)*nbFrames;

          cout<<"Loading number of atoms for each frame"<<endl;

            int nbAtomsAtFrame[nbFrames];

          loadNbAtomEachFrame(gropath, nbAtomsAtFrame, nbFrames);

          cout<<"Loading all lines in gro file"<<endl;

            vector<string> lineIn(bufferSize);
          LoadGroLines(gropath, lineIn, nbAtoms,nbFrames);


          cout<<"Loading coordinates"<<endl;

           vector<string> resNumbName(bufferSize);
           vector<string> names(bufferSize); vector<int> indexAtom(bufferSize);

           vector<float> z(bufferSize);

           LoadGro(gropath,resNumbName, nbAtoms,z, names, indexAtom, nbFrames );


        int it=0;int whereItWas=0;  int tempp;

        for(int l=0;l<nbFrames;l++){

                //timer start
                tim.start();

                nbAtoms=nbAtomsAtFrame[l];

                file<<firstline; file<<"   "<<frames[l]<<endl<<" "<<nbAtoms<<endl;

             for(int p= whereItWas ; p<whereItWas+nbAtoms ;p++){

                if(names[p] != NameToRemove ){ file<<lineIn[p]<<endl; }

                if(names[p]==NameToRemove && z[p]<Criteria){ if(l==0&& nbAtoms<150){cout<<names[p]<<" "<<indexAtom[p]<<" removed"<<endl;}cnt++;}

                if(names[p]==NameToRemove && z[p]>=Criteria){ file<<lineIn[p]<<endl;}

                                       tempp=p;}//end loop on one frame


                            whereItWas=tempp+2;
                  file<<lineIn[whereItWas-1]<<endl;
                  it++;

                  nbAtomsEachFrame[l]=nbAtoms-cnt;

                  if(l==0){cnt1=cnt;}

                  cnt=0;

                    if(l==0){ cout<<endl<<cnt1<<" "<<NameToRemove<<" removed"<<endl;
                    cout<<nbAtomsEachFrame[0]<<" Entities remaining"<<endl<<endl;}


                    if(l==0){cout<<"Processing for each frame -- Please Wait"<<endl<<endl; }


                           //GET TIME USED OF ONE FRAME LOOP

        finaltime=(tim.getElapsedTime(MILLISEC)/1000.0);

        tim.stop();

        //We calculate an average of the time to do one loop and then multiply by the loops remaining
        mean  += finaltime  ;

            chronot= (mean/(float)l)*(nbFrames-l);


        //Refresh the chrono every x time for long times and saving efficiency



         if( (l) % 10 == 0) {

            //Layout for longer times
         if(chronot>60 && chronot<3600){

                        minutes=chronot/60; hours = minutes / 60; temp5=(int)chronot%60;

                        cout<<"\rReading Time: "<< frames[l] <<" ------ Time Remaning :"<<(int)minutes<<" minutes and "<<temp5<<" secs "<<flush;

                                       }

         else if(chronot>3600 ){

                        minutes=chronot/60; hours = minutes / 60; temp5=(int)minutes%60;

                        cout<<"\rReading Time: "<< frames[l] <<" ------ Time Remaning :"<<(int)hours<<" hours and "<<temp5<<" minutes                          "<<flush;

                                       }

          else  {cout<<"\rReading Time: "<< frames[l] <<" ------ Time Remaning :"<<chronot<<" secs "<<flush;}




                           }

                           }  //end loop on each frame



                        cout<<" "<<flush;
                        cout<<"\rSaving File - please wait"<<endl<<endl;


             file.seekp(0, ios::beg);

             for(int l=0;l<nbFrames;l++){

                    getline(file,temp);
                    file2<<temp<<endl;
                    file2<<"  "<<nbAtomsEachFrame[l]<<endl;
                    getline(file,temp);

                for(int m=0;m<nbAtomsEachFrame[l]+1;m++){

                    getline(file,temp);
                    file2<<temp<<endl;

                                                        }

                                        }

               /* if( remove( lopath.c_str() ) != 0 )
                cout<< "Error deleting temporary file"<<endl;
                    else{}*/
cout<<"Completed"<<endl;

                }

        else{cout<<"Could not create the output file in function removeAll()."<<endl<<" We tried to create a file at this path : "<<lopath<<endl;}


}





//////////////
///////////     GET FUNCTIONS
///////////
///////////////



void loadNbAtomEachFrame(string gropath, int nbAtomAtFrame[],int nbFrames){

    ifstream file(gropath.c_str());
    char temp='j';
    int tempatom=0;
    float t=0;
    string temp2;

    if(file){

    for(int j=0;j<nbFrames;j++){      if(j%15==0)cout<<"\r frame -> "<< j <<"      "<<flush;

        while(temp != '='){ file.get(temp);  }
        file>>t;
        file>>tempatom;

        nbAtomAtFrame[j]=tempatom;
        getline(file,temp2);
        file.get(temp);

                            }    cout<<"\r"<<flush;
}

else{cout<<"ERROR "<<endl;}




}

string getFirstLineFromGro(string gropath){

    ifstream file(gropath.c_str());
    string temp;
    getline(file,temp);
    while(temp.back()!='=' ){ temp.erase(temp.size() - 1 ); }
    return temp;
}


void LoadGroLines( string gropath, vector<string> &lineIn, int nbAtoms, int nbFrames){

    ifstream file(gropath.c_str());

    char temp;
    bool continuee=true;
    int tempatom=0;
    float t=0;
    string temp2;
    int whereitwas=-1;
    int y, tempy;

    if(file){

    for(int k=0;k<nbFrames;k++){

        if(k%15==0)cout<<"\r frame -> "<< k <<"      "<<flush;

        getline(file,temp2);

        file>>tempatom;

        getline(file,temp2);

        for (y = whereitwas+1 ; y < whereitwas+tempatom+2 ; y++){ getline(file,temp2); lineIn[y]=temp2; tempy=y; }whereitwas=tempy;

        }
        cout<<"\r"<<flush;



                    }

        else{cout<<"error loading lines"<<endl;}

}



int getNbAtomAtFrame(std::string gropath, float frameToLook, float dt){

    ifstream file(gropath.c_str());
    char temp;
    bool continuee=true;
    int tempatom=0;
    float t=0;
    string temp2;


if(file){

    while(continuee)

        {

    while(temp != '='){ file.get(temp); }

    file>>t;

    if(t != frameToLook){

                        t=t+dt;

                        file.get(temp);

                        }

    else {

    file>>tempatom;

    continuee=false;
                            }



        }

        return tempatom;




        }

    else{cout<<"fr"<<endl;}

}



////////////////////////////////////////////////////////////////


void LoadGro(string path,vector<string> &resNumbName, int nbAtoms, vector<float> &z, vector<string> &names, std::vector<int>  &indexAtom, float nbFrames)

{


ifstream file(path.c_str());

if(file)
{


   string temp;
   string temp2; int i=0; float t=0;
   char tempc;
   bool continuee=true;
   int tempatom;
    int it=0,tempi=0;
    int whereitwas=-2;

    for(int k=0;k<nbFrames;k++){

        if(k%15==0)cout<<"\r frame -> "<< k <<"      "<<flush;

        getline(file,temp2);
        file>>tempatom;
        getline(file,temp2);

     for (int i = whereitwas+2 ; i < whereitwas+tempatom+2 ; i++){

    file>>resNumbName[i];
    file>>names[i]; //ALB
    file>>indexAtom[i]; //numero atome
    file>>temp;//x
    file>>temp; //y
    file>>z[i]; //z
    file>>temp; //vx
    file>>temp; //vy
    file>>temp; //vz
        tempi=i;
                        }whereitwas=tempi;

        getline(file,temp2); getline(file,temp2);



                                }      cout<<"\r"<<flush;

}

else{cout<<"Error loading .gro file"<<endl;}



}














//// OTHER STUFF

void changeGroFilePath(string *gropath, string *groname,string ipath){

    cout <<endl<<"Path of Configuration File .gro to modify : "<<endl<<endl;
    cout<<ipath;
    string temp;

    getline(cin,temp);

    *gropath=ipath+temp;
    temp=*gropath;

    if (temp.find(".gro") != string::npos) {
            cout << "found!" << '\n';

    temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );temp.erase(temp.size() - 1 );

    *groname=temp;

    cout<<endl;

}

    else{ cout<<"No such file found"<<endl;}

}


void loadNewGroConstants(int *nbAtoms,int *nbFrames,float *dtframe,float frames[], std::string gropath)
{

    string temp=gropath;
    bool cond=false;
    int tempNbAtoms=*nbAtoms;

    fstream test(gropath);
    if(test){cond=true;}
        test.close();


if(cond){
int temp;
cout<<"Counting Atoms"<<endl;
LoadnbAtoms(gropath,&temp);
*nbAtoms=temp;
cout<<"There is  "<<*nbAtoms<<" Atoms"<<endl;

//Count number of frames in the system
cout<<"Counting Frames"<<endl;
LoadnbFrames(gropath, &temp,tempNbAtoms);
*nbFrames=temp;
cout<<"There is  "<<*nbFrames<<" Time/frames"<<endl;

float temp2;
float temp3[*nbFrames];

Loadframe(gropath,temp3,*nbAtoms,&temp2);

*dtframe=temp2;
for(int l=0;l<*nbFrames;l++){frames[l]=temp3[l];}

}

}








