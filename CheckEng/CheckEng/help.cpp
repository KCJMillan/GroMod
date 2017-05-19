#include "help.h"
#include <iostream>
#include <string>

using namespace std;

void help(){

        bool b=true;
        string main;


        cout<<endl<<endl<<"   ------Gromacs analysis module help------"<<endl;
        cout<<"      Type 'exit' to exit the help"<<endl<<endl;
        cout<<"type help + <command> to have more details"<<endl<<endl;

        cout<<endl<<" *************"<<endl;
        cout<<"   FUNCTIONS"<<endl;
        cout<<" *************"<<endl<<endl;
        cout<<"  -neighbor :"<<endl<<"          Input a .gro file (for now) and will research the number of neighbors desired of each atoms. Output the file needed for all the following analysis functions. PBC is taken into acount when researching neighbor. Modifidy parameters of PBC in the file box.box. Modifidy parameters in the config file to ajust the research parameters."<<endl<<endl;
        cout<<"  -ion pair :"<<endl<<"          Take as Input the neighbor data file previously outputed. Look around user inputed ions and see if they are below a certain critical radius. PArameters can be changed in the ion pair configuration file. Output 3 files : Data file containing information about each ion and the analysis file showing the average of ion pairs, warnings, and more. Last file is a data file who can be plotted showing the average of ion pairs as a function of time."<<endl<<endl;
        cout<<"  -hydrogen bonds (still implementing):"<<endl<<"             .....    "<<endl<<endl;
        cout<<"  -sort:"<<endl<<"          The user can modify a .gro file before starting a neighbor research. Atoms can be removed, replace or add depending the user input condition (ex : remove all, remove only z<50...). Can be usefull to remove fixed atoms from a configuration or to compare ion pairing depending on z.   "<<endl<<endl;
        cout<<"     -sort --open <input>:"<<endl<<"          Subfunction of <sort>. Allows the user to open the .gro file before or after the sort."<<endl<<endl;
        cout<<"  -clear:"<<endl<<"        Clear the board. "<<endl<<endl;


        cout<<endl<<" *********************"<<endl;
        cout<<"  CONFIGURATION FILES"<<endl;
        cout<<" *********************"<<endl<<endl;
        cout<<"  -Config.conf : "<<endl<<"Main configuration file of the Neighbor Research. Containg Default paths and parameters such as number of atoms to look around. "<<endl;
        cout<<"      Location - /usr/local/GroMod/include/conf.gro "<<endl<<endl;
        cout<<"  -ionPairConfig.conf: Main configuration file for the ion pair function. Containg Default paths, parameters such as Critical Radius of ion pair, default output name.. "<<endl;
        cout<<"      Location - /usr/local/GroMod/include/ionPairConfig.conf"<<endl<<endl;
        cout<<"  -prog.conf : "<<endl<<"Configuration files for external softwares to use. For exemple : gnuplot to plot or gedit to read file, etc."<<endl;
        cout<<"      Location - /usr/local/GroMod/include/prog.conf"<<endl<<endl;
        cout<<"  -buffers.conf : "<<endl<<"Configuration files for machine effiency. Leave as default for most usages. However, this file as to be changed in some cases."<<endl;
        cout<<"      Location - /usr/local/GroMod/include/buffers.conf"<<endl<<endl;

        cout<<endl<<" *********************"<<endl;
        cout<<"     OTHER SETTINGS"<<endl;
        cout<<" *********************"<<endl<<endl;
        cout<<"  -boxConf.box :"<<endl<<"     The path of this file has to be provided in the Conf.conf file. It's used to correct the distances under PBC."<<endl;

        cout<<endl<<"***************"<<endl;
        cout<<">>> ";

         while(b){ getline(cin,main) ; if(main=="exit"){b=false;}else{cout<<">>> ";}  }



}
