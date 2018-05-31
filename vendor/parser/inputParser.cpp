#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "inputParser.h"

inline vector<double> Locs(double x, double y)
{
    vector<double> tmp{x, y};
    return tmp;
};

void print_inputs(INPUTS output)
{
    printf("lx = %1.2f, ly = %1.2f\n", output.lxy[0], output.lxy[1]);
    printf("E = %1.2f, nu = %1.2f, rho = %1.2f\n", output.E, output.nu, output.rho);

    printf("Hole configurations (# %d) >> \n", output.n_holes);
    for (int cc = 0; cc < output.n_holes; cc++)
        printf("x = %1.2f, y = %1.2f, d = %1.2f\n", output.hole_locs[cc][0], output.hole_locs[cc][1], output.hole_sizes[cc]);

    printf("BC (# %d) >> \n", output.n_BC);
    for (int cc = 0; cc < output.n_BC; cc++)
        printf("x = %1.2f, y = %1.2f, tolx = %1.2e, toly = %1.2e, dir = %i\n",
            output.BC_locs[cc][0], output.BC_locs[cc][1],
            output.BC_tols[cc][0], output.BC_tols[cc][1],
            output.BC_dircs[cc]);

    printf("Loads (# %d) >> \n", output.n_loads);
    for (int cc = 0; cc < output.n_loads; cc++)
        printf("x = %1.2f, y = %1.2f, tolx = %1.2f, toly = %1.2f, dir = %i, val = %1.2f\n",
            output.load_locs[cc][0], output.load_locs[cc][1],
            output.load_tols[cc][0], output.load_tols[cc][1],
            output.load_dircs[cc], output.load_vals[cc]);
    printf("movelimit = %1.2f\n", output.move_limit);
    printf("max_area = %1.2f\n", output.max_area);
}

INPUTS read_inputs(char *filename, bool echoFlag)
{
    // FILE *fid = fopen(filename, "r");
    INPUTS output;

    string line;
    ifstream inpfile(filename);

    if (inpfile.is_open())
    {
        while (!inpfile.eof())
        {
            streampos oldpos = inpfile.tellg();
            if (getline(inpfile, line))
            {
                if (line.find("*") == 0)
                {
                    // DESIGN DOMAIN ================================
                    if (line.compare(0, 3, "*Lx") == 0)
                    {
                        getline(inpfile, line, ',');
                        output.lxy[0] = stoi(line);
                        getline(inpfile, line, '\n');
                        output.lxy[1] = stoi(line);
                    }
                    if (line.compare(0, 6, "*E, nu") == 0)
                    {
                        getline(inpfile, line, ',');
                        output.E = stod(line);
                        getline(inpfile, line, ',');
                        output.nu = stod(line);
                        getline(inpfile, line, '\n');
                        output.rho = stod(line);
                    }

                    // HOLES ================================
                    if (line.compare(0, 10, "*num_holes") == 0)
                    {
                        getline(inpfile, line, '\n');
                        output.n_holes = stoi(line);
                        output.hole_locs.resize(output.n_holes, vector<double>(2, 0.0));
                        output.hole_sizes.resize(output.n_holes, 0.0);
                    }
                    if (line.compare(0, 9, "*hole_pos") == 0)
                    {
                        for (unsigned int cc = 0; cc < output.n_holes; cc++)
                        {
                            getline(inpfile, line, ',');
                            output.hole_locs[cc][0] = stod(line);
                            getline(inpfile, line, ',');
                            output.hole_locs[cc][1] = stod(line);
                            getline(inpfile, line, '\n');
                            output.hole_sizes[cc] = stod(line);
                        }
                    }

                    // Dirichlet BC ============================
                    if (line.compare(0, 7, "*num_BC") == 0)
                    {
                        getline(inpfile, line, '\n');
                        output.n_BC = stoi(line);
                        output.BC_locs.resize(output.n_BC, vector<double>(2, 0.0));
                        output.BC_tols.resize(output.n_BC, vector<double>(2, 0.0));
                        output.BC_dircs.resize(output.n_BC, 0);
                    }
                    if (line.compare(0, 7, "*BC_pos") == 0)
                    {
                        for (unsigned int cc = 0; cc < output.n_BC; cc++)
                        {
                            getline(inpfile, line, ',');
                            output.BC_locs[cc][0] = stod(line);
                            getline(inpfile, line, ',');
                            output.BC_locs[cc][1] = stod(line);
                            getline(inpfile, line, ',');
                            output.BC_tols[cc][0] = stod(line);
                            getline(inpfile, line, ',');
                            output.BC_tols[cc][1] = stod(line);
                            getline(inpfile, line, '\n');
                            output.BC_dircs[cc] = stoi(line);
                        }
                    }

                    // LOADS =====================================
                    if (line.compare(0, 10, "*num_loads") == 0)
                    {
                        getline(inpfile, line, '\n');
                        output.n_loads = stoi(line);
                        output.load_locs.resize(output.n_loads, vector<double>(2, 0.0));
                        output.load_tols.resize(output.n_loads, vector<double>(2, 0.0));
                        output.load_dircs.resize(output.n_loads, 0);
                        output.load_vals.resize(output.n_loads, 0.0);
                    }
                    if (line.compare(0, 9, "*load_pos") == 0)
                    {
                        for (unsigned int cc = 0; cc < output.n_loads; cc++)
                        {
                            getline(inpfile, line, ',');
                            output.load_locs[cc][0] = stod(line);
                            getline(inpfile, line, ',');
                            output.load_locs[cc][1] = stod(line);
                            getline(inpfile, line, ',');
                            output.load_tols[cc][0] = stod(line);
                            getline(inpfile, line, ',');
                            output.load_tols[cc][1] = stod(line);
                            getline(inpfile, line, ',');
                            output.load_dircs[cc] = stoi(line);
                            getline(inpfile, line, '\n');
                            output.load_vals[cc] = stod(line);
                        }
                    }

                    // Optimization problem
                    if (line.compare(0, 11, "*move_limit") == 0)
                    {
                        getline(inpfile, line, '\n');
                        output.move_limit = stod(line);
                    }
                    if (line.compare(0, 7, "*volCon") == 0)
                    {
                        getline(inpfile, line, '\n');
                        output.max_area = stod(line);
                    }
                }
            }
        }
        if (echoFlag)
            print_inputs(output);
        return output;
    }
    else
    {
        cout << "There is no such file"  << endl;
        cout << "Return to default cantilever...\n";

        // DOMAIN
        output.lxy[0] = 160;
        output.lxy[1] = 80;
        output.E = 1.0;
        output.nu = 0.3;
        output.rho = 1.0;

        // HOLES
        output.n_holes = 23;
        output.hole_sizes.resize(23, 5);

        output.hole_locs.clear();
        output.hole_locs.push_back(Locs(16, 14));
        output.hole_locs.push_back(Locs(48, 14));
        output.hole_locs.push_back(Locs(80, 14));
        output.hole_locs.push_back(Locs(112, 14));
        output.hole_locs.push_back(Locs(144, 14));

        output.hole_locs.push_back(Locs(32, 27));
        output.hole_locs.push_back(Locs(64, 27));
        output.hole_locs.push_back(Locs(96, 27));
        output.hole_locs.push_back(Locs(128, 27));

        output.hole_locs.push_back(Locs(16, 40));
        output.hole_locs.push_back(Locs(48, 40));
        output.hole_locs.push_back(Locs(80, 40));
        output.hole_locs.push_back(Locs(112, 40));
        output.hole_locs.push_back(Locs(144, 40));

        output.hole_locs.push_back(Locs(32, 53));
        output.hole_locs.push_back(Locs(64, 53));
        output.hole_locs.push_back(Locs(96, 53));
        output.hole_locs.push_back(Locs(128, 53));

        output.hole_locs.push_back(Locs(16, 66));
        output.hole_locs.push_back(Locs(48, 66));
        output.hole_locs.push_back(Locs(80, 66));
        output.hole_locs.push_back(Locs(112, 66));
        output.hole_locs.push_back(Locs(144, 66));

        for (int ii = 0; ii < output.hole_locs.size(); ii++)
        {
            cout << output.hole_locs[ii][0] << "\t" << output.hole_locs[ii][1] << endl;
        }

        // BC
        output.n_BC = 1;
        output.BC_locs.clear();
        output.BC_locs.push_back(Locs(0.0, 0.0));
        output.BC_tols.clear();
        output.BC_tols.push_back(Locs(1e-12, 1e10));
        output.BC_dircs.clear();
        output.BC_dircs.push_back(2);

        // Loads
        output.n_loads = 1;
        output.load_locs.clear();
        output.load_locs.push_back(Locs(160, 40));
        output.load_tols.clear();
        output.load_tols.push_back(Locs(1e-12, 1e-12));
        output.load_dircs.clear();
        output.load_dircs.push_back(1);
        output.load_vals.clear();
        output.load_vals.push_back(-0.5);

        // misc.
        output.max_area = 0.4;
        output.move_limit = 0.5;
        
        if (echoFlag)
            print_inputs(output);
        return output;
    }
}
