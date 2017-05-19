#include "const.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>


using namespace std;






///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////
/////  LOAD AND EDIT CONFIG AND BOX FUNCTIONS
/////
//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void LoadConfig(string *opath, string *ipath, int *nbAtomsTolook, string *oname, string *boxpath){



//ifstream file("/usr/local/GroMod/include/config.conf");

ifstream file("./ConfigFiles/config.conf");


if(file)

{
  cout<<"Config File Successfully loaded"<<endl;


  string temp;
  int N=100;


  for (int i=0;i<N;i++){

        getline(file,temp);

        if(temp=="#Default output Path"){getline(file,temp);while(temp==""){getline(file,temp);} *opath=temp;}

        if(temp=="#Default input Path"){getline(file,temp);while(temp==""){getline(file,temp);}  *ipath=temp;}

        if(temp=="#Default output name"){getline(file,temp);while(temp==""){getline(file,temp);}  *oname=temp;}

        if(temp=="#Default Box path"){getline(file,temp);while(temp==""){getline(file,temp);}  *boxpath=temp;}

        if(temp=="#Number of atoms to look around during the neighbord shearch"){getline(file,temp);while(temp==""){getline(file,temp);} *nbAtomsTolook = atof(temp.c_str()) ; }

  }

  file.close();

}

else{

cout<<"Error loading config file"<<endl;

}

}


//LOAD ION PAIR RESEARCH CONFIG FILE

void LoadPairConfig(string *opath, string *ipath, double *CritRadius, string *oname){



//ifstream file("/usr/local/GroMod/include/ionPairConfig.conf");

ifstream file("./ConfigFiles/ionPairConfig.conf");


if(file)

{
  cout<<"Ion Config File Successfully loaded"<<endl;


  string temp;
  int N=100;


  for (int i=0;i<N;i++){

        getline(file,temp);

        if(temp=="#Default output Path"){getline(file,temp);while(temp==""){getline(file,temp);} *opath=temp;}

        if(temp=="#Default input Path"){getline(file,temp);while(temp==""){getline(file,temp);}  *ipath=temp;}

        if(temp=="#Default output Name"){getline(file,temp);while(temp==""){getline(file,temp);}  *oname=temp;}

        if(temp=="#Critical Bond Lenght for ion pairing (nm)"){getline(file,temp);while(temp==""){getline(file,temp);} *CritRadius = atof(temp.c_str()) ; }

  }

  file.close();

}

else{

cout<<"Error loading config file"<<endl;

}



}






//LOAD EXTERNAL PROGRAMS CONFIG FILE

void LoadProgramConfig(std::string *prog1, std::string *prog2){


//ifstream file("/usr/local/GroMod/include/prog.conf");

ifstream file("./ConfigFiles/prog.conf");


if(file)

{
  string temp;
  int N=100;


  for (int i=0;i<N;i++){

        getline(file,temp);

        if(temp=="#To read file"){getline(file,temp);while(temp==""){getline(file,temp);} *prog1=temp;}

        if(temp=="#To plot file"){getline(file,temp);while(temp==""){getline(file,temp);}  *prog2=temp;}
  }

  file.close();

}

else{

cout<<"Error loading programs config file"<<endl;
}

}

void LoadBuffers(int *maxAtomBuff){


//ifstream file("/usr/local/GroMod/include/buffers.conf");

ifstream file("./ConfigFiles/buffers.conf");


if(file)

{
  string temp;
  int N=100;


  for (int i=0;i<N;i++){

        getline(file,temp);

        if(temp=="#Buffer Max Atom"){ getline(file,temp); while(temp==""){getline(file,temp);} *maxAtomBuff = atof(temp.c_str());}

                       }

  file.close();

}

else{

cout<<"Error loading buffers config file"<<endl;
}






}



//LOAD HBOND FUNCTION CONFIG

void LoadHbondConfig(std::string *opath, std::string *ipath, float *CritRadius, std::string *oname, float *CritAngle,std::string *oxy_name, std::string *h1_name,std::string *h2_name){



//ifstream file("/usr/local/GroMod/include/hbond.conf");

ifstream file("./ConfigFiles/hbond.conf");

if(file)

{

  string temp;
  int N=100;
  for (int i=0;i<N;i++){

        getline(file,temp);

        if(temp=="#Default output Path"){ getline(file,temp); while(temp==""){getline(file,temp);} *opath=temp;}

        if(temp=="#Default input Path"){ getline(file,temp); while(temp==""){getline(file,temp);} *ipath=temp;}

        if(temp=="#Default output Name"){ getline(file,temp); while(temp==""){getline(file,temp);} *oname=temp;}

        if(temp=="#Critical radius (nm)"){ getline(file,temp); while(temp==""){getline(file,temp);} *CritRadius = atof(temp.c_str());}

        if(temp=="#Bond angle (degres)"){ getline(file,temp); while(temp==""){getline(file,temp);} *CritAngle = atof(temp.c_str());}

        if(temp=="#Default Water-Oxygen name"){ getline(file,temp); while(temp==""){getline(file,temp);} *oxy_name = temp; }

        if(temp=="#Default Hydrogen1 name"){ getline(file,temp); while(temp==""){getline(file,temp);} *h1_name = temp; }

        if(temp=="#Default Hydrogen2 name"){ getline(file,temp); while(temp==""){getline(file,temp);} *h2_name = temp; }


                       }

  file.close();

}

else{

cout<<"Error loading hbond config file"<<endl;
}





}

//LOAD BOX FUNCTION


void LoadBoxSizeAngles(string gropath, string boxpath, float *a, float *b, float *c, float *alpha, float *beta, float *teta,int nbAtoms)

{

        //Box file
ifstream file(boxpath.c_str());

if(file)

{

  string temp;
  int N=100;


  for (int i=0;i<N;i++){

        getline(file,temp);

        if(temp=="#alpha"){getline(file,temp);while(temp==""){getline(file,temp);} *alpha=atof(temp.c_str()) ; }

        if(temp=="#beta"){getline(file,temp);while(temp==""){getline(file,temp);}  *beta=atof(temp.c_str()) ; }

        if(temp=="#gamma"){getline(file,temp);while(temp==""){getline(file,temp);}  *teta=atof(temp.c_str()) ; }

  }

  file.close();

}

else{

cout<<"Error loading box file"<<endl;

}


 //GRO FILE

ifstream file2(gropath.c_str());

if(file2)

{

  string temp;

        getline(file2,temp);
        getline(file2,temp);

  for (int i = 0 ;i<nbAtoms;i++){getline(file2,temp);}

        file2>>*a;file2>>*b;file2>>*c;

        file2.close();

}

else{

cout<<"Error loading gro file"<<endl;

}


}





///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////
/////////////                    COORDINATES TRANSFORMATION
/////////////
//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


void CalcTrigo(float alpha, float beta, float gamma, float *cos_alpha, float *sin_alpha, float *cos_beta, float *sin_beta, float *cos_gamma, float *sin_gamma, float *v){

    const float pi = 3.141592653589;

    float alphaTrig = alpha*pi/180.0 ;
    float betaTrig = beta*pi/180.0 ;
    float gammaTrig = gamma*pi/180.0 ;


    *cos_alpha=cos(alphaTrig);
    *cos_beta=cos(betaTrig);
    *cos_gamma=cos(gammaTrig);
    *sin_alpha=sin(alphaTrig);
    *sin_beta=sin(betaTrig);
    *sin_gamma=sin(gammaTrig);

    float ac= *cos_alpha;
    float bc= *cos_beta;
    float gc= *cos_gamma;

    *v = sqrt( 1.0 - ( ac * ac ) - ( bc * bc ) - ( gc * gc ) + ( 2.0*ac*bc*gc ) );

    }

void CalcMatrixElements(float a, float b, float c, float v, float cos_alpha, float sin_alpha, float cos_beta, float sin_beta, float cos_gamma, float sin_gamma, float *m11, float *m12, float *m13, float *m22, float *m23, float *m33){

    *m11=a;

    *m12=b*cos_gamma;

    *m22= b*sin_gamma;

    *m23=(   c * (cos_alpha - (cos_beta*cos_gamma)  ) ) / sin_gamma;

    *m33=(c*v) / sin_gamma;

    *m13=  c*cos_beta;

}


void CalcInverseMatrixElements(float a, float b, float c, float v, float cos_alpha, float sin_alpha, float cos_beta, float sin_beta, float cos_gamma, float sin_gamma, float *m11, float *m12, float *m13, float *m22, float *m23, float *m33){

    *m11=1.0 / a;

    *m12= -cos_gamma/(a*sin_gamma);

    *m22= 1.0 / (b*sin_gamma);

    *m23=-(cos_alpha-cos_beta*cos_gamma)/(v*b*sin_gamma);

    *m33=sin_gamma/(c*v);


    //very complicated expression so we split it out

    float bb=b*cos_gamma;
    float dd= (   c * (cos_alpha - (cos_beta*cos_gamma)  ) ) / sin_gamma;
    float cc= c*cos_beta;
    float tt= b*sin_gamma;
    float ee= (c*v) / sin_gamma;

//    *m13= ( (bb*dd) - (cc*tt) ) / (a*ee*tt) ;

     *m13= (cos_alpha*cos_gamma-cos_beta)/(a*v*sin_gamma)  ;


}



void TransformCoordinates(float *x, float *y, float *z , float m11, float m12, float m13, float m22, float m23, float m33){

    float old_x = *x;
    float old_y = *y;
    float old_z = *z;

    *x= (m11*old_x) + (m12*old_y) + (m13*old_z) ;
    *y= (m22*old_y) + (m23*old_z) ;
    *z= (m33*old_z) ;

}

void TransformAllCoordinates(int nbAtoms,float x[],float y[],float z[],float m11, float m12, float m13, float m22, float m23, float m33){


        float xtemp=0,ytemp=0,ztemp=0;

        for(int i=0;i<nbAtoms;i++)  {

            xtemp=x[i];
            ytemp=y[i];
            ztemp=z[i];

            TransformCoordinates(&xtemp,&ytemp,&ztemp,m11,m12,m13,m22,m23,m33);

            x[i]=xtemp;
            y[i]=ytemp;
            z[i]=ztemp;

                     }
}





/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////
/////////////                    FRAMES LOADING
/////////////
//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


//load nbFrames


void LoadnbFrames(std::string path , int  *nbFrames, int nbAtoms)

{

ifstream file(path.c_str());

int i=0;
int om=nbAtoms;

if(file)
{

    char temp;
    string temp2;
    // !file.eof()

     while ( file.get(temp) ) {

         if(temp == '='){ i++; if(i%15==0)cout<<"\r"<< i << flush;}

    }
        cout<<"\r"<<flush;
    *nbFrames=i;

    file.close();
}

else{cout<<"Error loading .gro file"<<endl;}

}


//LOAD FRAME

void Loadframe(string path, float time[], int nbAtoms, float *dtframe){

ifstream file(path.c_str());

int i=0;

if(file)
{

   char temp;

   while(!file.eof()){

         file.get(temp);

         if(temp == '='){file>>time[i]; i++; if(i%15==0)cout<<"\r"<< "frame -> "<< i << flush;}

         }

        cout<<"\r"<<flush;
    file.close();

    *dtframe=time[1]-time[0];
}

else{cout<<"Error loading .gro file"<<endl;}


}




///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////
/////////////                    FILE .GRO LOADING
/////////////
//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


//LOAD NB OF ATOMS

void LoadnbAtoms(string path ,int *nbAtoms){

ifstream file(path.c_str());

if(file)
{

   string temp;
  getline(file,temp);
  file>>*nbAtoms;
  file.close();
}

else{cout<<"Error loading .gro file"<<endl;}
}



//LOAD THE CONf.GRO FILE


void LoadGroFile(fstream& file, string path, int nbAtoms, float x[], float y[], float z[], float vx[], float vy[], float vz[], string names[], int indexAtom[], float frametoload,float dt)

{


//ifstream file(path.c_str());

if(file)
{


   char temp; int i=0; float t=0;
   string temp2;
   bool continuee=true;


    while(continuee){

        while(temp != '='){ file.get(temp);}

                    file>>t;

                    if(t != frametoload){

                        t=t+dt;
                        file.get(temp);

                        }


                    else if(t==frametoload){

                    file>>temp2;

                    for(i=0;i<nbAtoms;i++){

                    file>>temp2;//1ALB
                    file>>names[i]; //ALB
                    file>>indexAtom[i]; //numero atome
                    file>>x[i];//x
                    file>>y[i]; //y
                    file>>z[i]; //z
                    file>>vx[i]; //vx
                    file>>vy[i]; //vy
                    file>>vz[i]; //vz
                                           }

                    continuee=false; }

}

}

else{cout<<"Error loading .gro file"<<endl;}



}




//LOAD nbAtoms, nbFrames, dtFrames in THE NEIGHBOR DATA FILE

void LoadNeighborDataFileConstants(string path, int *nbAtoms,int *nbFrames,float *dtFrames){


ifstream file(path.c_str());

if(file)
{


    string temp;
    getline(file,temp);
    file>> *nbAtoms  ; file >>temp;
    file>> *nbFrames ; file >>temp;
    file>> *dtFrames ; file >>temp;

    file.close();
}

else{cout<<"Error loading the NEIHGBOR DATA file"<<endl;}

}


//LOAD THE NEIGHBOR DATA FILE

void LoadNeighborDataFile(fstream& file, string path, int nbAtoms, string atomName[], int indexAtom[],int ToNeighborIndex[], float ToNeighborLenght[] , string ToNeighborName[], int nbAtomsToLook, float dtFrames, float frameToLoad){

//ifstream file(path.c_str());

if(file)
{

    string temp;
    char ch;
    int i,k;
    int tempInt=0;
    float t;
    bool continuee=true;

    while(continuee){


        while(ch != '='){ file.get(ch);}

        file>>t;

        if(t!=frameToLoad){t= t + dtFrames; file>>temp; for(i=0;i<nbAtoms-2;i++){  getline(file,temp);    }  }


        else if(t==frameToLoad){

            for(int l=0 ; l<nbAtoms ; l++){

             file>>indexAtom[l]; //Atom number
             file>>atomName[l]; //atom name



                for( k=(0+l*nbAtomsToLook) ; k<nbAtomsToLook+nbAtomsToLook*l ; k++){

                    file>>ToNeighborIndex[k]; //Neighbor k index
                    file>>ToNeighborName[k]; //Neighbor Name
                    file>>ToNeighborLenght[k]; //distance to this neighbor

                                                         }

                                          }

                continuee=false;
                    }

}

}

else{cout<<"Error loading the NEIGHBOR DATA file"<<endl;}

}


void LoadNeighborFrameFromDataFile(string DataFilePath, int nbFrames,float frames[]){


    ifstream file(DataFilePath.c_str());
    char temp;
    bool continuee=true;
    int tempatom=0;
    float t=0;
    string temp2;


if(file){

    for(int m=0;m<nbFrames;m++){  if(m%15==0) cout<<"\r frame -> "<< m <<"      "<<flush;

    while(temp != '='){ file.get(temp); } file.get(temp);

    file>>t;
    frames[m]=t;

        } cout<<"\r"<<flush;

}

else{cout<<"fr"<<endl;}

}




///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////
/////////////                    ION PAIR FUNCTIONS
/////////////
//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


float CountPairFile(std::string pairPath, int pairCount){


    ifstream file(pairPath.c_str());

    if(file){

    }




}





