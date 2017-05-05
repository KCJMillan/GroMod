#include <string>
#include <cmath>
#include <fstream>


// Load or edit config file

void LoadConfig(std::string *opath, std::string *ipath, int *nbAtomsTolook, std::string *oname, std::string *boxpath);

void LoadPairConfig(std::string *opath, std::string *ipath, double *CritRadius, std::string *oname);

void LoadProgramConfig(std::string *prog1, std::string *prog2);

void LoadHbondConfig(std::string *opath, std::string *ipath, float *CritRadius, std::string *oname, float *CritAngle,std::string *oxy_name, std::string *h1_name,std::string *h2_name);

void LoadBuffers(int *maxAtomBuff);

// LOAD GRO FILE


void LoadGroFile(std::fstream& file, std::string path, int nbAtoms, float x[], float y[], float z[],float vx[],float vy[], float vz[], std::string names[],int indexAtom[], float frametoload,float dt);

void LoadnbAtoms(std::string path, int *nbAtoms);

void Loadframe(std::string path, float time[], int nbAtoms, float *dtframe);

void LoadnbFrames(std::string path , int  *nbFrames, int nbAtoms);


// LOAD THE BOX

void LoadBoxSizeAngles(std::string gropath, std::string boxpath, float *a, float *b, float *c, float *alpha, float *beta, float *teta, int nbAtoms);

//COORDINATES TRANSFORMATION

void CalcTrigo(float alpha, float beta, float gamma, float *cos_alpha, float *sin_alpha, float *cos_beta, float *sin_beta, float *cos_gamma, float *sin_gamma, float *v);

void CalcMatrixElements(float a, float b, float c, float v, float cos_alpha, float sin_alpha, float cos_beta, float sin_beta, float cos_gamma, float sin_gamma, float *m11, float *m12, float *m13, float *m22, float *m23, float *m33);

void CalcInverseMatrixElements(float a, float b, float c, float v, float cos_alpha, float sin_alpha, float cos_beta, float sin_beta, float cos_gamma, float sin_gamma, float *m11, float *m12, float *m13, float *m22, float *m23, float *m33);

void TransformCoordinates(float *x, float *y, float *z , float m11, float m12, float m13, float m22, float m23, float m33);

void TransformAllCoordinates(int nbAtoms,float x[],float y[],float z[],float m11, float m12, float m13, float m22, float m23, float m33);

void TransformBackCoordinates(float *x, float *y, float *z , float m11, float m12, float m13, float m22, float m23, float m33);


//LOAD BOND DATA FILE PRODUCED BY THIS ENGINE

void LoadNeighborDataFileConstants(std::string path, int *nbAtoms,int *nbFrames,float *dtFrames);

void LoadNeighborDataFile(std::fstream& file, std::string path, int nbAtoms, std::string atomName[], int indexAtom[], int ToNeighborIndex[], float ToNeighborLenght[] , std::string ToNeighborName[], int nbAtomsToLook, float dtFrames, float frameToLoad);

void LoadNeighborFrameFromDataFile(std::string DataFilePath, int nbFrames,float frames[]);

//ION PAIR FUNCTIONS

float CountPairFile(std::string pairPath, int pairCount);


