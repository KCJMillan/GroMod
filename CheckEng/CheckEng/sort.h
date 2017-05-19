
#include <vector>
#include <string>

//MENU FUNCTIONS
void sortRemove(std::string gropath,std::string groname, std::string opath,int nbAtoms,int nbFrames,float frames[], float dtframe);

//REMOVAL FUNCTIONS
void removeAll(std::string gropath,std::string opath,std::string groname,int nbAtoms,int nbFrames,float frames[], float dtframe, int nbToRemove, std::string NameToRemove);

void removeZinf(std::string gropath,std::string opath,std::string groname,int nbAtoms,int nbFrames,float frames[], float dtframe, int nbToRemove, std::string NameToRemove);

void removeZsup(std::string gropath,std::string opath,std::string groname, int nbAtoms,int nbFrames,float frames[], float dtframe, int nbToRemove, std::string NameToRemove);

//GET FUNCTIONS
std::string getFirstLineFromGro(std::string gropath);

int getNbAtomAtFrame(std::string gropath, float frameToLook, float dt);

void LoadGroLines(std::string gropath, std::vector<std::string> & lineIn, int nbAtoms, int nbFrames);

void loadNbAtomEachFrame(std::string gropath, int nbAtomAtFrame[], int nbFrames);



//OTHER STUFF
void changeGroFilePath(std::string *gropath, std::string *groname,std::string ipath);

void loadNewGroConstants(int *nbAtoms,int *nbFrames,float *dtframe,float frames[], std::string gropath);

void LoadGro(std::string path, std::vector<std::string> &resNumbName, int nbAtoms, std::vector<float> &z, std::vector<std::string> &names, std::vector<int> &indexAtom, float nbFrames);
