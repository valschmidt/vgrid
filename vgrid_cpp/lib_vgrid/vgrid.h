/* vgrid library */

#include <vector>
#include <iostream>
#include <string>
#include <cmath>


#ifndef VGRID_H_
#define VGRID_H_

using namespace std;

class vgrid
{
    double m_cs;
	double m_cinf;
	string m_gridtype;

    vector <double> m_xx;            // x coordinates of grid.
    vector <double> m_yy;            // y coordinates of grid.

    
	public:
        vgrid();
		vgrid(double cs, double cinf, string gridtype);
		~vgrid() {}
        
        bool m_verbose = 0;
        vector<vector <double> > m_zw;   // sum of z*weights.
        vector<vector <double> > m_ww;   // sum of weights.
        vector<vector <double> > m_nn;      // number of points included
        vector<vector <double> > m_varw; // sum of (x-mu)^2*weights.

        
        // Core Functions
        void hello_world(){ cout<<"Hello World" << "\n";}
        void create_grid();
        void expand_grid();
        void add(vector <double> x, vector <double> y, vector <double> z, vector <double> w);
        bool gridsizesanitycheck(vector <double> XX);

        // Gridding functions
        void grid_mean();
        
        // Retrieving data functions:
        vector <vector <double>> zz();
        vector <vector <double>> var();
        vector <vector <double>> std();
        
        // Utility Functions
        std::vector< double > arange(double start_val, double end_val, double increment);
        vector <double> concatinate(vector <double> x1, vector <double> x2);
        vector <vector <double>> concatinate2DUpDown(vector <vector <double>> upper, vector <vector <double>> lower);
        vector <vector <double>> concatinate2DLeftRight(vector <vector <double>> left, vector <vector <double>> right);
        void printvector(string name, vector <double> v);
        void printintvector(string name, vector <int> v);
        void print2Dvector(string name, vector <vector <double>> v);

    private:

        int m_grows = 0;                // number of rows in grid.
        int m_gcols = 0;                // number of columns in grid.
        int m_row = 0;                  // row that's being operated on.
        int m_col = 0;                  // col that's being operated on.
        vector <int> m_indices2add;     // indices of new data that will be assymulated into node at (row,col)
        

    
        vector <double> m_x;             // x values to add.
        vector <double> m_y;             // y values to add.
        vector <double> m_z;             // z values to add.
        vector <double> m_w;             // weight value to add.
    
};

#endif /* VGRID_H_ */
