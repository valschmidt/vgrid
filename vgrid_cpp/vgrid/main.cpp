#include <iostream>
#include <string>
#include <vector>
#include "../lib_vgrid/vgrid.h"
#include <unistd.h>

using namespace std;

int main(int argc, char **argv) {

	int opt;
    bool testflag = 0;
	
	while ((opt = getopt(argc, argv, "th")) != EOF )
	{
	
		switch (opt)
		{
		case 't':   // test
			{
				cout << "Testing vgrid.\n";
                testflag = 1;
				break;
			}
		case 'h':
			{
				cout << "vgrid usage: \n";
				break;
			}
		
		
		}
	}
	 
    // Create test data.
    int N = 10;
    vector <double> x;
    vector <double> y;
    vector <double> z;
    vector <double> w;
    
    int i = 0;
    for (i=0; i < N; i++){
     x.push_back(double(i));
     y.push_back(double(i));
     z.push_back(double(i));
     w.push_back(1);
    }
    
    //vgrid G = vgrid(1.0,1.0,"mean");
    vgrid G = vgrid(1.0,3.0,"mean");
    G.m_verbose = 1;
    G.add(x,y,z,w);
	G.hello_world();
    /*vector <double> x2 = {0};
    vector <double> y2 = {15.0};
    vector <double> z2 = {0};
    vector <double> w2 = {0};
    G.add(x2,y2,z2,w2);*/
    G.print2Dvector("ww",G.m_ww);
    G.print2Dvector("zw",G.m_zw);
    G.print2Dvector("zz",G.zz());
    G.print2Dvector("std",G.std());
	return 0;
}
