#ifndef M2DO_FEA_INPUTPARSER_H
#define M2DO_FEA_INPUTPARSER_H
using namespace std;

struct INPUTS{
    double lxy[2];
    double E, nu, rho;
    
    int n_holes;
    vector<vector<double> > hole_locs;
    vector<double> hole_sizes;

    int n_BC;
    vector<vector<double> > BC_locs;
    vector<vector<double> > BC_tols;
    vector<int> BC_dircs;

    int n_loads;
    vector<vector<double> > load_locs;
    vector<vector<double> > load_tols;
    vector<int> load_dircs;
    vector<double> load_vals;

    double move_limit, max_area;      
};

INPUTS read_inputs(char * filename, bool echoFlag = true);

#endif