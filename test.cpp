#include <string>
#include <iostream>
#include "test.h"
#include <vector>

using namespace std;


void test(){

        int ol=5 ;
        string test[ol][ol];

        for(int i=0;i<5;i++){

            for(int g = 1; g<5; g++){

                    test[i][g]="ol";


            }

        }

        cout<<test[1][2];



}
