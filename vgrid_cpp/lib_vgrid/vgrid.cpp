/* vgrid Library */

#include "vgrid.h"
#include <algorithm>


using namespace std;

/* Default class constructor */
vgrid::vgrid(){
 
    m_cs = 1.0;
    m_cinf = 1.0;
    m_gridtype = "mean";
    
}
/* Class constructor */
vgrid::vgrid(double cs, double cinf, string gridtype){
    
    m_cs = cs;
    m_cinf = cinf;
    m_gridtype = gridtype;
    
}

/* Verify the grid won't be too large. */
bool vgrid::gridsizesanitycheck(vector <double> XX){
    if (XX.size() > 1e4){
        return 0;
    }else {
        return 1;
    }
}

/* A utilty function to get a range of values */
/* When increment > 0, result = [start_val, >= end_val] */
/* When increment < 0, result = [<= start , end_val] */
vector <double> vgrid::arange(double start_val, double end_val, double increment){
    
    vector <double> result;
    double val;
    
    /* In increment > 0, values return start at start value and end after end 
    value is exceeded. If increment < 0, values start at end value and decrement 
    until values are less than start_value. */
    if (increment > 0.0) {
        val = start_val;
        while (val <= (end_val)){
            result.push_back(val);
            val += increment;
        }
    }
    else {
        val = end_val;
        while (val > start_val){
            result.insert(result.begin(),val);
            val += increment;
        }
    }
    return result;
}

/* A utility funtion to return a concatenation of 2 vectors */
vector <double> vgrid::concatenate(vector <double> x1, vector <double> x2){
 
    vector <double> result(x1);
    result.insert(result.end(), x2.begin(), x2.end());
    return result;
}

vector <vector <double>> vgrid::concatenate2DUpDown(vector <vector <double>> upper, vector <vector <double>> lower){
 
    if (upper[0].size() != lower[0].size()) {
        cout << "ERROR attempted vertical concatenation of 2D vectors having differing numbers of columns!\n";
        exit(1);
    }
    vector <vector <double>> result(upper);
    result.insert(result.end(), lower.begin(), lower.end());
    return result;
}

vector <vector <double>> vgrid::concatenate2DLeftRight(vector <vector <double>> left, vector <vector <double>> right){
 
    if (left.size() != right.size()) {
        cout << "ERROR attempted horizontal concatenation of 2D vectors having differing numbers of rows.\n";
        exit(1);
    }
    vector <vector <double>> result(left);
    for(int i = 0; i < right.size(); i++){
        result[i].insert(result[i].end(), right[i].begin(), right[i].end());
    }
    return result;
}

/* Creates a new grid based on initial points added */
void vgrid::create_grid(){
    
    cout << "m_cs:" << m_cs << "\n";
    cout << "m_cimf:" << m_cinf << "\n";
    cout << "m_gridtype:" << m_gridtype << "\n";
    
    
    // Get the bounds of the data.
    double min_x = *min_element(m_x.begin(), m_x.end());
    double max_x = *max_element(m_x.begin(), m_x.end());
    double min_y = *min_element(m_y.begin(), m_y.end());
    double max_y = *max_element(m_y.begin(), m_y.end());
    
    // Create the extents of x and y.
    double val = min_x;
    m_xx = arange(min_x - m_cinf, max_x + m_cinf, m_cs);
    m_yy = arange(min_y - m_cinf, max_y + m_cinf, m_cs);
    /*while (val <= max_x + m_cs){
     m_xx.push_back(val);
     val += m_cs;
    }
    */
    
    m_gcols = m_xx.size();
    m_grows = m_yy.size();
    
    // Initialize memory for grid.
    m_zw.resize(m_grows, vector <double> (m_gcols, NAN));
    m_ww.resize(m_grows, vector <double> (m_gcols, NAN));    
    m_nn.resize(m_grows, vector <double> (m_gcols, 0));
    m_varw.resize(m_grows, vector <double> (m_gcols, NAN));
    
    /*
    for (int i = 0; i<m_xx.size(); i++) {
        cout << m_xx[i] << "\n";
    }
    */
    
    if (m_verbose){
        cout << "Created new grid.\n";
        
        cout << "Grid is " << m_xx.size(); 
        cout << "x" << m_yy.size() << ".\n";

        cout << "   Cell Spacing: " << m_cs << "\n";
        cout << "   Cell Influence: " << m_cinf << "\n";
        cout << "   Grid type:  " << m_gridtype << "\n";
    }
    
}


/* Expand the grid based on the bounds of the new data.*/
void vgrid::expand_grid(){
    
    // Get the bounds of the new data.
    double min_x = *min_element(m_x.begin(), m_x.end());
    double max_x = *max_element(m_x.begin(), m_x.end());
    double min_y = *min_element(m_y.begin(), m_y.end());
    double max_y = *max_element(m_y.begin(), m_y.end());
    
    vector <double> gridextension;
    
    // Expand to the west..
    if (min_x < m_xx[0]){
        gridextension = arange(min_x, m_xx.front() - m_cs, -m_cs);
        m_xx = concatenate(gridextension, m_xx);

        /* Expand grid...*/
        vector <vector <double>> tmp (m_zw.size(), vector <double> (gridextension.size(), NAN));
        m_zw = concatenate2DLeftRight(tmp, m_zw);
        m_ww = concatenate2DLeftRight(tmp, m_ww);
        m_varw = concatenate2DLeftRight(tmp, m_varw);
        m_nn = concatenate2DLeftRight(tmp, m_nn);


    }
    // Expand to the east...
    if (max_x > m_xx[m_xx.size()-1]) {
        gridextension = arange(m_xx.back() + m_cs, max_x, m_cs);
        m_xx = concatenate(m_xx, gridextension);
        /* Expand grid...*/
        vector <vector <double>> tmp (m_zw.size(), vector <double> (gridextension.size(), NAN));
        m_zw = concatenate2DLeftRight(m_zw, tmp);
        m_ww = concatenate2DLeftRight(m_ww, tmp);
        m_varw = concatenate2DLeftRight(m_varw, tmp);
        m_nn = concatenate2DLeftRight(m_nn, tmp);

    }
    // Expand to the south...
        if (min_y < m_yy[0]){
        gridextension = arange(min_y, m_yy.front() - m_cs, -m_cs);
        m_yy = concatenate(gridextension, m_yy);
        /* Expand grid here...*/
        vector <vector <double>> tmp (gridextension.size(), vector <double> (m_zw[0].size(), NAN));
        m_zw = concatenate2DUpDown(m_zw,tmp);
        m_ww = concatenate2DUpDown(m_ww,tmp);
        m_varw = concatenate2DUpDown(m_varw,tmp);
        m_nn = concatenate2DUpDown(m_nn,tmp);

    }
    // Expand to the north...
    if (max_y > m_yy[m_yy.size()-1]) {
        gridextension = arange(m_yy.back() + m_cs, max_y, m_cs);
        m_yy = concatenate(m_yy, gridextension);
        //printvector("m_yy",m_yy);
        vector <vector <double>> tmp (gridextension.size(), vector <double> (m_zw[0].size(), NAN));
        m_zw = concatenate2DUpDown(tmp,m_zw);
        m_ww = concatenate2DUpDown(tmp,m_ww);
        m_varw = concatenate2DUpDown(tmp,m_varw);
        m_nn = concatenate2DUpDown(tmp,m_nn);


    }
    
    m_grows = m_yy.size();
    m_gcols = m_xx.size();
    

    if (m_verbose){
        cout << "Expanded grid.\n";
        cout << "Grid is " << m_yy.size(); 
        cout << "x" << m_xx.size() << ".\n";   
    }    
    
}

void vgrid::printvector(string name, vector <double> v){
    cout << "vector: " << name << "\n";
    for(int ii = 0; ii < v.size(); ii++){
        cout << v[ii] << "\n";
    }
}

void vgrid::printintvector(string name, vector <int> v){
    cout << "vector: " << name << "\n";
    for(int ii = 0; ii < v.size(); ii++){
        cout << v[ii] << "\n";
    }
}

void vgrid::print2Dvector(string name, vector <vector <double>> v){
    cout << "vector: " << name << "\n";
    
    for(int ii = 0; ii < v.size(); ii++){
        for(int jj = 0; jj < v[0].size(); jj++){
            cout << v[ii][jj] << "   ";
        }
        cout << "\n";
    }
}

void vgrid::printGridInfo(){
    /*"Print basic info about the grid."*/

    cout << "Grid Info:" << "\n";
    cout << "  Grid Cell Size: " << m_cs << "\n";
    cout << "  Grid Cell Influence:" << m_cinf << "\n";
    cout << "  Grid Size: " << m_grows << "x" << m_gcols << "\n";
    cout << "  Type: " << m_gridtype << "\n";
    cout << "  Grid Bounds: X: " << m_xx[0] << "-" << m_xx[m_xx.size()-1] 
                        << "\tY: " << m_yy[0] << "-" << m_yy[m_yy.size()-1] << "\n"; 
}

void vgrid::add(vector <double> x, vector <double> y, vector <double> z, vector <double> w){

    int idx, jdx;
    int i = 0;
    m_x = x;
    m_y = y;
    m_z = z;
    m_w = w;
    int N = m_x.size();
    bool foundinrow = 0;
    double cinf2 = m_cinf * m_cinf;
    double ddy;
    vector <int> yidx, xidx;
    vector <double> xtest;
    
    // Check to make sure all have the same number of values.
    if (m_y.size() != N | m_z.size() != N | m_w.size() != N){
        cout << "Error! x,y and z are not the same size. Aborting add()!\n";
        return;
    }
    
    // Make sure weights are not 0.
    for (i = 0; i < m_w.size(); i++){
        if (m_w[i] < 1e-20){
            m_w[i] = 1e-20;
        }
    }
    
    // Create grid or expand it. 
    if (m_grows == 0) {
        create_grid();
    }else{
        expand_grid();
    }
    

    // Loop through the rows of the grid...
    for(m_row = 0; m_row < m_grows; m_row++){
        
        // This next line handles the fact that "y" indexing is + down 
        // c++ but 
        //m_row = (m_grows - 1) - idx;   
        foundinrow = 0;
        xtest.clear();
        yidx.clear();
        
        // Check all points against the N/S test for closeness...
        for(i = 0; i < N; i++){
            ddy = (m_y[i] - m_yy[m_row]) * (m_y[i] - m_yy[m_row]);
            if (ddy < cinf2){ 
                yidx.push_back(i);
                xtest.push_back(cinf2 - ddy);
                foundinrow = 1;
            }
        }
        
        // If none are found, skip this row.
        if (!foundinrow){
            continue;
        }
        

        // Loop through the cols of the grid.
        for(m_col = 0; m_col < m_gcols; m_col++){

            m_indices2add.clear();
 
            // For values that passed the N/S test, conduct the E/W test...
            for(i = 0; i< yidx.size(); i++) {
                if ( ((m_x[yidx[i]] - m_xx[m_col]) * (m_x[yidx[i]] - m_xx[m_col])) < xtest[i]) {
                    m_indices2add.push_back(yidx[i]);
                }
            }
            
            m_nn[m_row][m_col] += m_indices2add.size();

            //printintvector("m_indices2add",m_indices2add);
            // Do the grid calculations here.
            // Note we don't have to pass the grid node, because we've 
            // set member variables m_row and m_col already.
            if (m_gridtype == "mean"){
                grid_mean();
            }
            
        } // End loop through cols (xx)
    } // End loop through rows (yy)
}

/* Calculate grided values on the fly when required */
vector <vector <double>> vgrid::zz(){
        vector <vector <double>> result (m_grows, vector <double> (m_gcols, NAN));

        for(int idx = 0; idx < m_grows; idx++){
            for(int jdx = 0; jdx < m_gcols; jdx++){
             result[idx][jdx] = m_zw[idx][jdx] / m_ww[idx][jdx];
            }
        }
        return result;
}

vector <vector <double>> vgrid::var(){
        vector <vector <double>> result (m_grows, vector <double> (m_gcols, NAN));
        for(int idx = 0; idx < m_grows; idx++){
            for(int jdx = 0; jdx < m_gcols; jdx++){
             result[idx][jdx] = m_varw[idx][jdx] / m_ww[idx][jdx];
            }
        }
        return result;
    
}

vector <vector <double>> vgrid::std(){
        vector <vector <double>> result (m_grows, vector <double> (m_gcols, NAN));
        for(int idx = 0; idx < m_grows; idx++){
            for(int jdx = 0; jdx < m_gcols; jdx++){
             result[idx][jdx] = sqrt(m_varw[idx][jdx] / m_ww[idx][jdx]);
            }
        }
        return result;
}

void vgrid::grid_mean(){
    for (int i = 0; i< m_indices2add.size(); i++){
                
        // Calculate sum of z*w
        if (isnan(m_zw[m_row][m_col])){
            m_zw[m_row][m_col] = m_z[m_indices2add[i]] * m_w[m_indices2add[i]];
        }else{
            m_zw[m_row][m_col] += m_z[m_indices2add[i]] * m_w[m_indices2add[i]];
        }
        
        // Calculate sum of w
        if (isnan(m_ww[m_row][m_col])){
            m_ww[m_row][m_col] = m_w[m_indices2add[i]];
        }else{
            m_ww[m_row][m_col] += m_w[m_indices2add[i]];
        }
        
        // Calculate sum of (z - mu[i][j])^2 * w for weighted variance.
        if (isnan(m_varw[m_row][m_col])){
            m_varw[m_row][m_col] = (m_z[m_indices2add[i]] - m_zw[m_row][m_col]) * 
                                    (m_z[m_indices2add[i]] - m_zw[m_row][m_col]) 
                                    * m_w[m_indices2add[i]];
        }else{
            m_varw[m_row][m_col] += (m_z[m_indices2add[i]] - m_zw[m_row][m_col]) * 
                                    (m_z[m_indices2add[i]] - m_zw[m_row][m_col]) 
                                    * m_w[m_indices2add[i]];
        }
  
        
    }
    
}
