#include <iostream>
#include <math.h>
#include "conf.h"
#include "const.h"
#include "help.h"
#include "sort.h"
#include "hbond.h"
#include "test.h"
#include "ionpair.h"

using namespace std;

int main()
{

    string main="hello";
    cout<<endl<<endl;
    cout << " --------- Gromacs Help Analysis Engine ----------"<<endl<<endl;
    cout << " -------Created at Temple University : ICMS --------"<<endl;
    cout << "     ---------Type help for help---------" <<endl<<endl;cout<<">>";


    while(main!="exit")


    {

        getline(cin,main);

     //MAIN LOOP, MENU SELECTION

     if(main=="neighbor"){bondlengthAnalysis();cout<<endl;cout<<">>";}

     else if(main=="ion pair"){/* ionPairAnalysis();cout<<">>";*/  test3(); }

     else if(main=="sort"){ sortConf();cout<<">>"; }

     else if(main=="hbond"){ hbond();cout<<">>"; }

     else if(main=="test"){ test3(); }

     else if(main=="clear"){for(int i=0;i<100;i++){cout<<endl;}cout<<">>";}

     else if (main=="help"){help();cout<<endl<<" --------- Gromacs Help Analysis Engine ----------"<<endl<<">>"; }

     else {if(main!="exit"){cout<<">>";}}




    }

    return 0;
}
