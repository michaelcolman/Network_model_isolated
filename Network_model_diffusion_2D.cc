// Simplest implementation for ease of including in other code

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <math.h>

using namespace std;

// values for geo 50_50_1
#define NX 50
#define NY 50

// Reused utility functions =================================================\\|
double rush_larsen(double y, double ss, double tau, double dt)
{
    double gate;
    gate = ss - (ss-y)*exp(-dt/tau);
    return gate;
}

double sigmoid(double V, double V_half, double k)
{
    double value;
    value = 1/(1 + exp((V - V_half)/k) );
    return value;
}
// End reused utility functions ==============================================//|


int main()
{
    // Files
    FILE * single_cell_out;
    single_cell_out = fopen("Vm_cell.txt", "wt");

    FILE * vtk_out;
    char *string = (char*)malloc(500);

    char string_geo[100]              = "Fibrosis_geometry_50_50_1.txt";
    char string_orientation[100]      = "Fibrosis_geometry_50_50_1_orientation.txt";

    // Global simulation varaibles
    double t, dt, total_time;
    int iteration_count = 0;
    total_time = 300;
    dt = 0.02; // ms

    // Global model and state variables
    double Vm[NX][NY];
    double Vm_next[NX][NY];

    double OX[NX][NY]; // orientation, x-component
    double OY[NX][NY]; // orientation, y-component
    double OZ[NX][NY]; // orientation, z-component

    int     geo[NX][NY];

    // NETWORK MODEL variables =============================================================================\\|
    double dx, dy, dz;
    dx = dy = dz = 0.25; // mm

    // Only used to define conductances
    double D1, D_AR;
    D1 = 0.2;
    D_AR = 5;

    double Gl, Gt;          // long(axial) and transverse junction conductance
    Gl = D1/dx;       // Can define from D1/dx if D1 is known already
    Gt = Gl/D_AR;

    double  gGap_node_xx[NX][NY]; // conductance in direction xx (based on orientation)
    double  gGap_node_yy[NX][NY];
    double  gGap_node_xypp[NX][NY];
    double  gGap_node_xypm[NX][NY];

    int Njunc;  // Number of junctions in the model
    // END NETWORK MODEL variables =========================================================================//|

    // Cell model / reaction variables (just needed to test actual implementation) =========================\\|
    // conductance |
    double gIp0d             = 16;       // s/mF
    double gIp3r            = 0.05;     // s/mF
    double gIp4r            = 0.3;      // s/mF

    // variables 
    double Ip0d;                // The current itself
    double Ip0d_va_ss;          // voltage activation, steady state
    double Ip0d_va_tau;         // voltage activation, time constant
    double Ip0d_vi_1_ss;        // voltage inactivation, steady state
    double Ip0d_vi_1_tau;       // voltage inactivation, time constant
    double Ip0d_vi_2_ss;        // voltage inactivation, steady state
    double Ip0d_vi_2_tau;       // voltage inactivation, time constant
    double Ip0d_va_al;          // voltage activation, alpha transition rate (1-y-> y)
    double Ip0d_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    double Ip0d_vi_1_al;        // voltage activation, alpha transition rate (1-y -> y)
    double Ip0d_vi_1_bet;       // voltage activation, beta transition rate  (y -> 1-y)
    double Ip0d_vi_2_al;        // voltage activation, alpha transition rate (1-y -> y)
    double Ip0d_vi_2_bet;       // voltage activation, beta transition rate  (y -> 1-y)

    double Ip3r;                // The current itself
    double Ip3r_va_ss;          // voltage activation, steady state
    double Ip3r_va_tau;         // voltage activation, time constant
    double Ip3r_va_al;          // voltage activation, alpha transition rate (1-y -> y)
    double Ip3r_va_bet;         // voltage activation, beta transition rate  (y -> 1-y)
    double Ip3r_vi_ti;          // Time-independent inactivation gate

    double Ip4r;                // The current itself

    double Istim;
    double Itot;

    // State variables
    double Ip0d_va[NX][NY];             // voltage activation
    double Ip0d_vi_1[NX][NY];           // voltage inactivation 1
    double Ip0d_vi_2[NX][NY];           // voltage inactivation 2
    double Ip2r_va[NX][NY];             // voltage activation
    double Ip2r_vi[NX][NY];             // voltage activation
    double Ip3r_va[NX][NY];             // voltage activation
    // End Cell model / reaction variables  ================================================================//|

    // Initial conditions of state variables ===============================================================\\|
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            Vm[i][j]               = -85;
            Ip0d_va[i][j]         = 0;
            Ip0d_vi_1[i][j]        = 1;
            Ip0d_vi_2[i][j]        = 1;
            Ip3r_va[i][j]          = 0;
        }
    }
    // End Initial conditions of state variables ===========================================================//|

    // Tissue geometry and fibres ==========================================================================\\|
    FILE * geo_file;
    FILE * orientation_file;

    geo_file = fopen(string_geo, "r");
    orientation_file = fopen(string_orientation, "r");

    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            /*geo[i][j] = 1;
              OX[i][j] = 0;
              OY[i][j] = 1;
              OZ[i][j] = 0;*/
            fscanf(geo_file, "%d ", &geo[i][j]);
            fscanf(orientation_file, "%lf %lf %lf ", &OX[i][j], &OY[i][j], &OZ[i][j]);
        }
    }
    fclose(geo_file);
    fclose(orientation_file);
    // End tissue geometry and fibres ======================================================================//|

    // Setup NETWORK MODEL =================================================================================\\|
    // Setup gGap for each node and direction based on orientation ================================\\|
    double theta_x_y;   // angle in x-y plane
    double W_axis_x_y; // Weight towards principal (axis; away from diagonal) in x-y plane
    double W_axis;     // Weight towards principal direction
    double W_diag_plane; // Diagonal, in plane
    double W_xx, W_yy;
    double W_xypp, W_xypm;
    double x, y; // local myocyte orientation components
    double x_t, y_t; // transverse componetns of myocyte orientation

    // "Pointing" parameters (1 if pointing in said direction, 0 otherwise)
    double Px, Py;
    double Pxyp;

    // Weights for transverse connections
    double W_t_xx, W_t_yy;
    double W_t_xypp, W_t_xypm;

    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            if (geo[i][j] > 0)
            {
                // local copy of orientation components
                x = OX[i][j];
                y = OY[i][j];

                // Default all weights to 0, so only need to add to relevant ones
                W_axis = W_diag_plane = 0;
                W_xx = W_yy = 0;
                W_xypp = W_xypm = 0;

                W_t_xx = W_t_yy = 0;
                W_t_xypp = W_t_xypm =  0;

                // Default gGap to zero so also can just be added where relevant
                gGap_node_xx[i][j] = gGap_node_yy[i][j] = 0;
                gGap_node_xypp[i][j] = gGap_node_xypm[i][j] = 0;

                // Default all pointing parameters to 0
                Px = Py = Pxyp = 0;

                if (x*x + y*y > 0) // do fibres exist?
                {
                    // First, calculate theta and weights in each plane
                    theta_x_y       = asin(fabs(x)/(sqrt(x*x + y*y)));
                    if (x == 0 && y == 0)  theta_x_y = 0; // override any nans from a /0
                    W_axis_x_y    = fabs(theta_x_y - 0.25*M_PI)/(0.25*M_PI);

                    // Set pointing parameters
                    // Main axes (only one can be non-zero at a time)
                    if (x*x >= y*y)         Px = 1;
                    else                    Py = 1;

                    // Diagonals (more than one can be non-zero)
                    if (x*y >= 0)          Pxyp = 1;

                    // Assign quadrant weights
                    W_axis            = W_axis_x_y;
                    W_diag_plane        = (1-W_axis_x_y);

                    // Assign weights to directions
                    // Primary (axial)
                    W_xx                = Px        * W_axis;    // only non-zero if pointing primarily towards x
                    W_yy                = Py        * W_axis;

                    W_xypp              = Pxyp      * W_diag_plane;
                    W_xypm              = (1-Pxyp)  * W_diag_plane;

                    // Transverse
                    W_t_xx              = Py        * W_axis;
                    W_t_yy              = Px        * W_axis;

                    W_t_xypp            = (1-Pxyp)  * W_diag_plane;
                    W_t_xypm            = Pxyp      * W_diag_plane;
                } // end checking if fibre exists
                else // if no fibre
                {
                    printf("ERROR: No fibre found for node %d %d || model will work but is desgined for anisotropic simulations\n", i, j);
                }

                // Set all below to 1.0 if just want standard implementation; these factors just help ensure symmetries
                double symm_fac_diag_long           = 1.0; 
                double symm_fac_diag_trans          = 1.0;

                // Assign gGap based on weights
                gGap_node_xx[i][j]        = Gt*W_t_xx     + Gl*W_xx; // simply = weight in transverse * gt + weight along fibre * Gl
                gGap_node_yy[i][j]        = Gt*W_t_yy     + Gl*W_yy; // note that both weights can be 0, giving no contribution to said direction

                gGap_node_xypp[i][j]      = (Gt*W_t_xypp)*symm_fac_diag_trans   + Gl*W_xypp * symm_fac_diag_long;
                gGap_node_xypm[i][j]      = (Gt*W_t_xypm)*symm_fac_diag_trans   + Gl*W_xypm * symm_fac_diag_long;

                // Geometry scaling
                // if you want this indepnendet of dx, then remove all dx, dy, dz below, and include 1/sqrt(2) and 1/sqrt(3) factors in diagonals and corners
                gGap_node_xx[i][j]     *= 1.0/dx;
                gGap_node_yy[i][j]     *= 1.0/dy;

                gGap_node_xypp[i][j]   *= 1.0/sqrt(dx*dx + dy*dy);
                gGap_node_xypm[i][j]   *= 1.0/sqrt(dx*dx + dy*dy);

            } // end if geo is > 1 loop
        } // end x
    } //end y
    // End Setup gGap for each node and direction based on orientation ============================//|

    // Calc Njunctions =============================================================================\\|
    Njunc = 0;
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            if (geo[i][j] > 0) // if it is an actual cell/node
            {
                if (i < NX-1 &&              geo[i+1][j] > 0)         Njunc++;
                if (j < NY-1 &&              geo[i][j+1] > 0)         Njunc++;
                if (i < NX-1 && j < NY-1 &&  geo[i+1][j+1] > 0)       Njunc++;
                if (i > 0    && j < NY-1 &&  geo[i-1][j+1] > 0)       Njunc++;
            } // end if
        } // end x
    } // end y
    printf("Number of junctions = %d\n", Njunc);
    // End calc Njunctions =========================================================================//|

    // Allocate arrays which are of size Njunc =====================================================\\|
    double  IDiff;     // junction current (not needed as an array in this implementation, but can be if desired)
    double  *gGap_jn;  // junction conductance
    int     *jn_map_minus_x;  // returns Ncell array index of the minus cell for IDiff[i][j]
    int     *jn_map_minus_y;  // returns Ncell array index of the minus cell for IDiff[i][j]
    int     *jn_map_plus_x;
    int     *jn_map_plus_y;
    int     *connection_type_jn; // whether the jucntion is comprised transverse, long, or a mix

    gGap_jn               = new double [Njunc];
    //IDiff               = new double [Njunc];
    jn_map_minus_x        = new int [Njunc];
    jn_map_plus_x         = new int [Njunc];
    jn_map_minus_y        = new int [Njunc];
    jn_map_plus_y         = new int [Njunc];
    connection_type_jn  = new int [Njunc];
    // End Allocate arrays which are of size Njunc =================================================//|

    // Create junctions and relate to nodes ========================================================\\|
    bool xx;
    bool yy;
    bool xypp;
    bool xypm;

    int count_junc = 0;

    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            if (geo[i][j] > 0) // if it is an actual cell/node
            {
                xx      = false;
                yy      = false;
                xypp    = false;
                xypm    = false;

                if (i < NX-1 &&                                         geo[i+1][j] > 0)         xx      = true;
                if (j < NY-1 &&                                         geo[i][j+1] > 0)         yy      = true;
                if (i < NX-1 && j < NY-1 &&                             geo[i+1][j+1] > 0)       xypp    = true;
                if (i > 0        && j < NY-1 &&                         geo[i-1][j+1] > 0)       xypm    = true;

                // Calculate junction g as an average of the contribution from the two nodes that comprise that junction
                // Principal directions
                if (xx == true) // i.e., if there is a connection in the x direction
                {
                    gGap_jn[count_junc]        = (gGap_node_xx[i][j] + gGap_node_xx[i+1][j])/2;

                    // "minus" is always current cell
                    jn_map_minus_x[count_junc]  = i; // seems redundant here, but later we loop over junctions so need this
                    jn_map_minus_y[count_junc]  = j; 

                    jn_map_plus_x[count_junc]   = i+1;
                    jn_map_plus_y[count_junc]   = j;

                    count_junc++; // now increments so that next junction is unique
                } 
                if (yy == true) // i.e., if there is a connection in the y direction
                {
                    gGap_jn[count_junc]        = (gGap_node_yy[i][j] + gGap_node_yy[i][j+1])/2;

                    jn_map_minus_x[count_junc]  = i; 
                    jn_map_minus_y[count_junc]  = j;

                    jn_map_plus_x[count_junc]   = i;
                    jn_map_plus_y[count_junc]   = j+1;

                    count_junc++;
                }
                if (xypp == true) 
                {
                    gGap_jn[count_junc]        = (gGap_node_xypp[i][j] + gGap_node_xypp[i+1][j+1])/2;

                    jn_map_minus_x[count_junc]  = i;
                    jn_map_minus_y[count_junc]  = j;

                    jn_map_plus_x[count_junc]   = i+1;
                    jn_map_plus_y[count_junc]   = j+1;

                    count_junc++; 
                }
                if (xypm == true)
                {
                    gGap_jn[count_junc]        = (gGap_node_xypm[i][j] + gGap_node_xypm[i-1][j+1])/2;

                    jn_map_minus_x[count_junc]  = i;
                    jn_map_minus_y[count_junc]  = j;

                    jn_map_plus_x[count_junc]   = i-1;
                    jn_map_plus_y[count_junc]   = j+1;

                    count_junc++; 
                }
            } // end if
        } // end x
    } // end y
    printf("second count of juncs = %d\n", count_junc);
    // End Set g for actual junctions and maps to relate junctions to nodes ========================//|
    // End Setup NETWORK MODEL =============================================================================//|

    // Time loop ===========================================================================================\\|
    for (t = 0; t < total_time; t+=dt)
    {
        // Spatial loop ===============================================================================\\|
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                // Cell/reaction model ==============================================================\\|
                // Set Activation gate alpha and beta
                Ip0d_va_al                 = 0.32*(Vm[i][j]+47.13)/(1-exp(-0.1*(Vm[i][j]+47.13)));
                Ip0d_va_bet                = 0.08*exp(-Vm[i][j]/11);

                // Set inactivation gates alphas and betas
                if (Vm[i][j] < -40.0)
                {
                    Ip0d_vi_1_al           = 0.135*exp((80+Vm[i][j])/-6.8);
                    Ip0d_vi_1_bet          = 3.56*exp(0.079*Vm[i][j])+310000*exp(0.35*Vm[i][j]);
                    Ip0d_vi_2_al           = (-127140*exp(0.2444*Vm[i][j])-0.00003474*exp(-0.04391*Vm[i][j]))*((Vm[i][j]+37.78)/(1+exp(0.311*(Vm[i][j]+79.23))));
                    Ip0d_vi_2_bet          = (0.1212*exp(-0.01052*Vm[i][j]))/(1+exp(-0.1378*(Vm[i][j]+40.14)));
                }
                else
                {
                    Ip0d_vi_1_al           = 0;
                    Ip0d_vi_1_bet          = 1/(0.13*(1+exp((Vm[i][j]+10.66)/-11.1)));
                    Ip0d_vi_2_al           = 0;
                    Ip0d_vi_2_bet          = (0.3*exp(-0.0000002535*Vm[i][j]))/(1+exp(-0.1*(Vm[i][j]+32)));
                }

                // Set tau and SS from alpha and beta
                Ip0d_va_tau                = 1/(Ip0d_va_al + Ip0d_va_bet); // 1/(a+b)
                Ip0d_vi_1_tau              = 1/(Ip0d_vi_1_al + Ip0d_vi_1_bet);
                Ip0d_vi_2_tau              = 1/(Ip0d_vi_2_al + Ip0d_vi_2_bet);
                Ip0d_va_ss                 = Ip0d_va_al * Ip0d_va_tau; // a*tau
                Ip0d_vi_1_ss               = Ip0d_vi_1_al * Ip0d_vi_1_tau;
                Ip0d_vi_2_ss               = Ip0d_vi_2_al * Ip0d_vi_2_tau;

                Ip0d_va[i][j]                      = rush_larsen(Ip0d_va[i][j], Ip0d_va_ss, Ip0d_va_tau, dt); // lib/Membrane.c
                Ip0d_vi_1[i][j]                    = rush_larsen(Ip0d_vi_1[i][j], Ip0d_vi_1_ss, Ip0d_vi_1_tau, dt);
                Ip0d_vi_2[i][j]                    = rush_larsen(Ip0d_vi_2[i][j], Ip0d_vi_2_ss, Ip0d_vi_2_tau, dt);

                Ip0d       = gIp0d * pow(Ip0d_va[i][j], 3) * Ip0d_vi_1[i][j] * Ip0d_vi_2[i][j] * (Vm[i][j] - 76);

                // 3r
                // Alpha
                if (fabs(Vm[i][j]+14.1)< 1e-10)        Ip3r_va_al  = 0.0015; // Denominator = 0 clause
                else                                Ip3r_va_al  = 0.0003*(Vm[i][j]+14.1)/(1-exp((Vm[i][j]+14.1)/-5));

                // Beta
                if (fabs(Vm[i][j]-3.3328) < 1e-10)     Ip3r_va_bet = 3.7836118e-4;
                else                                Ip3r_va_bet = 0.000073898*(Vm[i][j]-3.3328)/(exp((Vm[i][j]-3.3328)/5.1237)-1);

                // time constant and steady state
                Ip3r_va_tau        =  1/(Ip3r_va_al+Ip3r_va_bet);
                Ip3r_va_ss         = sigmoid(Vm[i][j], -14.10, -6.5); // V, V1/2, k

                // Time-independent inactivation gate
                Ip3r_vi_ti         = sigmoid(Vm[i][j], -15, 22.4);

                Ip3r_va[i][j]            = rush_larsen(Ip3r_va[i][j], Ip3r_va_ss, Ip3r_va_tau, dt);

                Ip3r               = gIp3r * Ip3r_va[i][j] * Ip3r_vi_ti * (Vm[i][j] - (-88));

                Ip4r               = gIp4r * (Vm[i][j] - (-88))/((1+exp(0.07*(Vm[i][j] - (-80)))));

                Itot = Ip0d + Ip3r + Ip4r;
                // End cell/reaction model ==========================================================//|

                // Istim
                if (t < 2.0 && i < 10 && j < 10) Istim = -20;
                else Istim = 0;

                // Update voltage
                Vm_next[i][j] = Vm[i][j] - dt*(Itot + Istim);

                // Output single cell file
                if (iteration_count % (int)(1/dt) == 0)
                {
                    if (i == 3 && j == 3) fprintf(single_cell_out, "%f %f\n", t, Vm_next[i][j]); 
                    if (i == 3 && j == 3) printf("%f %f\n", t, Vm_next[i][j]);
                }
            }
        }
        // End Spatial loop ===========================================================================//|

        // NETWORK MODEL coupling =====================================================================\\|
        int i1, i2, j1, j2, k1, k2;
        for (int n = 0; n < Njunc; n++) // looping over junctions (not nodes)
        {
            i1 = jn_map_minus_x[n];
            j1 = jn_map_minus_y[n];

            i2 = jn_map_plus_x[n];
            j2 = jn_map_plus_y[n];

            IDiff    = gGap_jn[n] * (Vm[i2][j2] - Vm[i1][j1]);

            // update voltages of nodes connected to that junction
            Vm_next[i2][j2] = Vm_next[i2][j2] + dt*(-IDiff);
            Vm_next[i1][j1] = Vm_next[i1][j1] + dt*(IDiff);
        }
        // End NETWORK MODEL coupling =================================================================//|

        // assign Vm to updated Vm_next
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                Vm[i][j] = Vm_next[i][j];
            }
        }

        // Output spatial vtk
        if ((iteration_count/5) % (int)(1/dt) == 0)
        {
            sprintf(string, "Vm_%.0f.vtk", t);
            vtk_out = fopen(string, "wt");

            fprintf(vtk_out, "# vtk DataFile Version 3.0\n");
            fprintf(vtk_out, "vtk output\n");
            fprintf(vtk_out, "ASCII\n");
            fprintf(vtk_out, "DATASET STRUCTURED_POINTS\n");
            fprintf(vtk_out, "DIMENSIONS %d %d %d\n", NX, NY, 1);
            fprintf(vtk_out, "SPACING 1 1 1\n");
            fprintf(vtk_out, "ORIGIN 0 0 0\n");
            fprintf(vtk_out, "POINT_DATA %d\n", NX*NY*1);
            fprintf(vtk_out, "SCALARS geometry float 1\n");
            fprintf(vtk_out, "LOOKUP_TABLE default\n");

            for (int j = 0; j < NY; j++)
            {
                for (int i = 0; i < NX; i++)
                {
                    fprintf(vtk_out, "%f ", Vm[i][j]);
                }
                fprintf(vtk_out, "\n");
            }
            fclose(vtk_out);
        }

        iteration_count++;
    }
    // End time loop =======================================================================================//|

    fclose(single_cell_out);

    // Deallocate dynamic arrays
    delete []   gGap_jn;
    //delete []   IDiff;
    delete []   jn_map_minus_x;
    delete []   jn_map_plus_x;
    delete []   jn_map_minus_y;
    delete []   jn_map_plus_y;
    delete []   connection_type_jn;

} // end main
