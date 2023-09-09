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
#define NZ 1

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
    double Vm[NX][NY][NZ];
    double Vm_next[NX][NY][NZ];

    double OX[NX][NY][NZ]; // orientation, x-component
    double OY[NX][NY][NZ]; // orientation, y-component
    double OZ[NX][NY][NZ]; // orientation, z-component

    int     geo[NX][NY][NZ];

    // NETWORK MODEL variables =============================================================================\\|
    double dx, dy, dz;
    dx = dy = dz = 0.25; // mm
    
    // Only used to define conductances
    double D1, D_AR;
    D1 = 0.2;
    D_AR = 5;

    double Gl, Gt, Gt2;          // long(axial) and transverse junction conductance
    Gl = D1/dx;       // Can define from D1/dx if D1 is known already
    Gt = Gl/D_AR;
    Gt2 = Gt; // set to smaller if three values are used; set to Gt to make it primary myocyte orientation only

    double  gGap_node_xx[NX][NY][NZ]; // conductance in direction xx (based on orientation)
    double  gGap_node_yy[NX][NY][NZ];
    double  gGap_node_zz[NX][NY][NZ];
    double  gGap_node_xypp[NX][NY][NZ];
    double  gGap_node_xypm[NX][NY][NZ];
    double  gGap_node_xzpp[NX][NY][NZ];
    double  gGap_node_xzpm[NX][NY][NZ];
    double  gGap_node_yzpp[NX][NY][NZ];
    double  gGap_node_yzpm[NX][NY][NZ];
    double  gGap_node_xyzppp[NX][NY][NZ];
    double  gGap_node_xyzppm[NX][NY][NZ];
    double  gGap_node_xyzpmp[NX][NY][NZ];
    double  gGap_node_xyzmpp[NX][NY][NZ];

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
    double Ip0d_va[NX][NY][NZ];             // voltage activation
    double Ip0d_vi_1[NX][NY][NZ];           // voltage inactivation 1
    double Ip0d_vi_2[NX][NY][NZ];           // voltage inactivation 2
    double Ip2r_va[NX][NY][NZ];             // voltage activation
    double Ip2r_vi[NX][NY][NZ];             // voltage activation
    double Ip3r_va[NX][NY][NZ];             // voltage activation
    // End Cell model / reaction variables  ================================================================//|

    // Initial conditions of state variables ===============================================================\\|
    for (int k = 0; k < NZ; k++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                Vm[i][j][k]               = -85;
                Ip0d_va[i][j][k]         = 0;
                Ip0d_vi_1[i][j][k]        = 1;
                Ip0d_vi_2[i][j][k]        = 1;
                Ip3r_va[i][j][k]          = 0;
            }
        }
    }
    // End Initial conditions of state variables ===========================================================//|

    // Tissue geometry and fibres ==========================================================================\\|
    FILE * geo_file;
    FILE * orientation_file;

    geo_file = fopen(string_geo, "r");
    orientation_file = fopen(string_orientation, "r");

    for (int k = 0; k < NZ; k++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                /*geo[i][j][k] = 1;
                OX[i][j][k] = 0;
                OY[i][j][k] = 1;
                OZ[i][j][k] = 0;*/
                fscanf(geo_file, "%d ", &geo[i][j][k]);
                fscanf(orientation_file, "%lf %lf %lf ", &OX[i][j][k], &OY[i][j][k], &OZ[i][j][k]);
            }
        }
    }
    fclose(geo_file);
    fclose(orientation_file);
    // End tissue geometry and fibres ======================================================================//|
    
    // Setup NETWORK MODEL =================================================================================\\|
    // Setup gGap for each node and direction based on orientation ================================\\|
    double theta_x_y;   // angle in x-y plane
    double theta_z_x;
    double theta_z_y;
    double W_axis_x_y; // Weight towards principal (axis; away from diagonal) in x-y plane
    double W_axis_x_z; // Weight towards principal (axis; away from diagonal) in x-z plane
    double W_axis_y_z; // Weight towards principal (axis; away from diagonal) in y-z plane
    double W_axis;     // Weight towards principal direction
    double W_diag_plane; // Diagonal, in plane
    double W_diag_elevation; // Diagonal, elevated
    double W_diag_diag; // corner
    double W_xx, W_yy, W_zz;
    double W_xypp, W_xypm, W_xzpp, W_xzpm, W_yzpp, W_yzpm;
    double W_xyzppp, W_xyzppm, W_xyzpmp, W_xyzmpp;
    double x, y, z; // local myocyte orientation components
    double x_t, y_t, z_t; // transverse 1
    double x_t2, y_t2, z_t2; // transverse 2
    double theta, phi;
    double theta_t, phi_t;
    double theta_t2, phi_t2;

    // "Pointing" parameters (1 if pointing in said direction, 0 otherwise)
    double Px, Py, Pz;
    double Pxyp, Pxzp, Pyzp;
    double Px_t, Py_t, Pz_t;
    double Pxyp_t, Pxzp_t, Pyzp_t;
    double Px_t2, Py_t2, Pz_t2;
    double Pxyp_t2, Pxzp_t2, Pyzp_t2;

    // Weights for transverse connections
    double W_t_xx, W_t_yy, W_t_zz;
    double W_t_xypp, W_t_xypm, W_t_xzpp, W_t_xzpm, W_t_yzpp, W_t_yzpm;
    double W_t_xyzppp, W_t_xyzppm, W_t_xyzpmp, W_t_xyzmpp;
    double W_t2_xx, W_t2_yy, W_t2_zz;
    double W_t2_xypp, W_t2_xypm, W_t2_xzpp, W_t2_xzpm, W_t2_yzpp, W_t2_yzpm;
    double W_t2_xyzppp, W_t2_xyzppm, W_t2_xyzpmp, W_t2_xyzmpp;

    for (int k = 0; k < NZ; k++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                if (geo[i][j][k] > 0)
                {
                    // local copy of orientation components
                    x = OX[i][j][k];
                    y = OY[i][j][k];
                    z = OZ[i][j][k];

                    // Transverse orientation components
                    // NOTE: if your data has the three eigenvectors, simply set x_t, x_t2 etc from those components as above
                    // The below is for cases where we need to calculate the transverse vectors
                    if (z == 1) // special case as theta undefined. But easy as set to other axes
                    {
                        x_t     = 1;
                        y_t     = 0;
                        z_t     = 0;

                        x_t2    = 0;
                        y_t2    = 1;
                        z_t2    = 0;
                    }
                    else // general case
                    {
                        // First, transform by 90 degree rotation in phi
                        phi = asin(z);
                        phi_t2 = phi + M_PI/2;
                        x_t2 = (x/cos(phi))*cos(phi_t2);
                        y_t2 = (y/cos(phi))*cos(phi_t2);
                        z_t2 = sin(phi_t2);

                        // Now, find normal to that vector and the original (from cross product of the 2)
                        x_t = y*z_t2 - z*y_t2;
                        y_t = z*x_t2 - x*z_t2;
                        z_t = x*y_t2 - y*x_t2;
                    }

                    // Default all weights to 0, so only need to add to relevant ones
                    W_axis = W_diag_plane = W_diag_elevation = W_diag_diag = 0;
                    W_xx = W_yy = W_zz = 0;
                    W_xypp = W_xypm = W_xzpp = W_xzpm = W_yzpp = W_yzpm = 0;
                    W_xyzppp = W_xyzppm = W_xyzpmp = W_xyzmpp = 0;

                    W_t_xx = W_t_yy = W_t_zz = 0;
                    W_t_xypp = W_t_xypm = W_t_xzpp = W_t_xzpm = W_t_yzpp = W_t_yzpm = 0;
                    W_t_xyzppp = W_t_xyzppm = W_t_xyzpmp = W_t_xyzmpp = 0;
                    W_t2_xx = W_t2_yy = W_t2_zz = 0;
                    W_t2_xypp = W_t2_xypm = W_t2_xzpp = W_t2_xzpm = W_t2_yzpp = W_t2_yzpm = 0;
                    W_t2_xyzppp = W_t2_xyzppm = W_t2_xyzpmp = W_t2_xyzmpp = 0;

                    // Default gGap to zero so also can just be added where relevant
                    gGap_node_xx[i][j][k] = gGap_node_yy[i][j][k] = gGap_node_zz[i][j][k] = 0;
                    gGap_node_xypp[i][j][k] = gGap_node_xypm[i][j][k] = gGap_node_xzpp[i][j][k] = 0;
                    gGap_node_xzpm[i][j][k] = gGap_node_yzpp[i][j][k] = gGap_node_yzpm[i][j][k] = 0;
                    gGap_node_xyzppp[i][j][k] = gGap_node_xyzppm[i][j][k] = gGap_node_xyzpmp[i][j][k] = gGap_node_xyzmpp[i][j][k] = 0;

                    // Default all pointing parameters to 0
                    Px = Py = Pz = Pxyp = Pxzp = Pyzp = 0;
                    Px_t = Py_t = Pz_t = Pxyp_t = Pxzp_t = Pyzp_t = 0;
                    Px_t2 = Py_t2 = Pz_t2 = Pxyp_t2 = Pxzp_t2 = Pyzp_t2 = 0;

                    if (x*x + y*y + z*z > 0) // do fibres exist?
                    {
                        // Primary eigenvector ===========================================================\\|
                        // First, calculate theta and weights in each plane
                        theta_x_y       = asin(fabs(x)/(sqrt(x*x + y*y)));
                        theta_z_x       = asin(fabs(z)/(sqrt(z*z + x*x)));
                        theta_z_y       = asin(fabs(z)/(sqrt(z*z + y*y)));
                        if (x == 0 && y == 0)  theta_x_y = 0; // override any nans from a /0
                        if (x == 0 && z == 0)  theta_z_x = 0;
                        if (y == 0 && z == 0)  theta_z_y = 0;
                        W_axis_x_y    = fabs(theta_x_y - 0.25*M_PI)/(0.25*M_PI);
                        W_axis_x_z    = fabs(theta_z_x - 0.25*M_PI)/(0.25*M_PI);
                        W_axis_y_z    = fabs(theta_z_y - 0.25*M_PI)/(0.25*M_PI);

                        //printf("%f %f %f %f %f %f\n", theta_x_y, theta_z_x, theta_z_y, W_axis_x_y, W_axis_x_z, W_axis_y_z);

                        // Set pointing parameters
                        // Main axes (only one can be non-zero at a time)
                        if (x*x >= y*y && x*x >= z*z)   Px = 1;
                        else if (y*y >= z*z)            Py = 1;
                        else                            Pz = 1;

                        // Diagonals (more than one can be non-zero)
                        if (x*y >= 0)                   Pxyp = 1;
                        if (x*z >= 0)                   Pxzp = 1;
                        if (y*z >= 0)                   Pyzp = 1;

                        // Assign quadrant weights
                        W_axis            = Px*W_axis_x_y*W_axis_x_z          + Py*W_axis_x_y*W_axis_y_z          + Pz*W_axis_x_z*W_axis_y_z;
                        W_diag_plane        = Px*(1-W_axis_x_y)*(W_axis_x_z)    + Py*(1-W_axis_x_y)*(W_axis_y_z)    + Pz*(1-W_axis_x_z)*W_axis_y_z;
                        W_diag_elevation    = Px*W_axis_x_y*(1-W_axis_x_z)      + Py*W_axis_x_y*(1-W_axis_y_z)      + Pz*(W_axis_x_z)*(1-W_axis_y_z);
                        W_diag_diag         = Px*(1-W_axis_x_y)*(1-W_axis_x_z)  + Py*(1-W_axis_x_y)*(1-W_axis_y_z)  + Pz*(1-W_axis_x_z)*(1-W_axis_y_z);

                        //printf("W %f %f %f %f\n", W_axis, W_diag_plane, W_diag_elevation, W_diag_diag);

                        // Assign weights to directions
                        W_xx                = Px * W_axis;    // only non-zero if pointing primarily towards x
                        W_yy                = Py * W_axis;
                        W_zz                = Pz * W_axis;

                        W_xypp              = (Px + Py) * Pxyp      * W_diag_plane;
                        W_xypm              = (Px + Py) * (1-Pxyp)  * W_diag_plane;

                        W_xzpp              = Px * Pxzp     * W_diag_elevation      + Pz * Pxzp     * W_diag_plane;
                        W_xzpm              = Px * (1-Pxzp) * W_diag_elevation      + Pz * (1-Pxzp) * W_diag_plane;

                        W_yzpp              = (Py + Pz) * Pyzp      * W_diag_elevation;
                        W_yzpm              = (Py + Pz) * (1-Pyzp)  * W_diag_elevation;

                        W_xyzppp            = Pxzp      *   Pyzp    * W_diag_diag;
                        W_xyzpmp            = (1-Pxyp)  *   Pxzp    * W_diag_diag;
                        W_xyzmpp            = (1-Pxyp)  *  (1-Pxzp) * W_diag_diag;
                        W_xyzppm            = Pxyp      *  (1-Pxzp) * W_diag_diag;
                        // End Primary eigenvector =======================================================//|

                        // Transverse directions =========================================================\\|
                        // transverse 1 =================================================\\|
                        theta_x_y       = asin(fabs(x_t)/(sqrt(x_t*x_t + y_t*y_t)));
                        theta_z_x       = asin(fabs(z_t)/(sqrt(z_t*z_t + x_t*x_t)));
                        theta_z_y       = asin(fabs(z_t)/(sqrt(z_t*z_t + y_t*y_t)));
                        if (x_t == 0 && y_t == 0)  theta_x_y = 0; // override any nans from a /0
                        if (x_t == 0 && z_t == 0)  theta_z_x = 0;
                        if (y_t == 0 && z_t == 0)  theta_z_y = 0;
                        W_axis_x_y    = fabs(theta_x_y - 0.25*M_PI)/(0.25*M_PI);
                        W_axis_x_z    = fabs(theta_z_x - 0.25*M_PI)/(0.25*M_PI);
                        W_axis_y_z    = fabs(theta_z_y - 0.25*M_PI)/(0.25*M_PI);

                        //printf("%f %f %f %f %f %f\n", theta_x_y, theta_z_x, theta_z_y, W_axis_x_y, W_axis_x_z, W_axis_y_z);

                        if (x_t*x_t >= y_t*y_t && x_t*x_t >= z_t*z_t)   Px_t = 1;
                        else if (y_t*y_t >= z_t*z_t)                    Py_t = 1;
                        else                                            Pz_t = 1;

                        // Diagonals (more than one can be non-zero)
                        if (x_t*y_t >= 0)                   Pxyp_t = 1;
                        if (x_t*z_t >= 0)                   Pxzp_t = 1;
                        if (y_t*z_t >= 0)                   Pyzp_t = 1;

                        // Assign quadrant weights
                        W_axis            = Px_t*W_axis_x_y*W_axis_x_z          + Py_t*W_axis_x_y*W_axis_y_z          + Pz_t*W_axis_x_z*W_axis_y_z;
                        W_diag_plane        = Px_t*(1-W_axis_x_y)*(W_axis_x_z)    + Py_t*(1-W_axis_x_y)*(W_axis_y_z)    + Pz_t*(1-W_axis_x_z)*W_axis_y_z;
                        W_diag_elevation    = Px_t*W_axis_x_y*(1-W_axis_x_z)      + Py_t*W_axis_x_y*(1-W_axis_y_z)      + Pz_t*(W_axis_x_z)*(1-W_axis_y_z);
                        W_diag_diag         = Px_t*(1-W_axis_x_y)*(1-W_axis_x_z)  + Py_t*(1-W_axis_x_y)*(1-W_axis_y_z)  + Pz_t*(1-W_axis_x_z)*(1-W_axis_y_z);

                        //printf("W %f %f %f %f\n", W_axis, W_diag_plane, W_diag_elevation, W_diag_diag);

                        // Assign weights to directions
                        W_t_xx                = Px_t * W_axis;    // only non-zero if pointing primarily towards x
                        W_t_yy                = Py_t * W_axis;
                        W_t_zz                = Pz_t * W_axis;

                        W_t_xypp              = (Px_t + Py_t) * Pxyp_t      * W_diag_plane;
                        W_t_xypm              = (Px_t + Py_t) * (1-Pxyp_t)  * W_diag_plane;

                        W_t_xzpp              = Px_t * Pxzp_t     * W_diag_elevation      + Pz_t * Pxzp_t     * W_diag_plane;
                        W_t_xzpm              = Px_t * (1-Pxzp_t) * W_diag_elevation      + Pz_t * (1-Pxzp_t) * W_diag_plane;

                        W_t_yzpp              = (Py_t + Pz_t) * Pyzp_t      * W_diag_elevation;
                        W_t_yzpm              = (Py_t + Pz_t) * (1-Pyzp_t)  * W_diag_elevation;

                        W_t_xyzppp            = Pxzp_t      *   Pyzp_t    * W_diag_diag;
                        W_t_xyzpmp            = (1-Pxyp_t)  *   Pxzp_t    * W_diag_diag;
                        W_t_xyzmpp            = (1-Pxyp_t)  *  (1-Pxzp_t) * W_diag_diag;
                        W_t_xyzppm            = Pxyp_t      *  (1-Pxzp_t) * W_diag_diag;
                        // End transverse 1 =============================================//|

                        // Transverse 2 =================================================\\|
                        theta_x_y       = asin(fabs(x_t2)/(sqrt(x_t2*x_t2 + y_t2*y_t2)));
                        theta_z_x       = asin(fabs(z_t2)/(sqrt(z_t2*z_t2 + x_t2*x_t2)));
                        theta_z_y       = asin(fabs(z_t2)/(sqrt(z_t2*z_t2 + y_t2*y_t2)));
                        if (x_t2 == 0 && y_t2 == 0)  theta_x_y = 0; // override any nans from a /0
                        if (x_t2 == 0 && z_t2 == 0)  theta_z_x = 0;
                        if (y_t2 == 0 && z_t2 == 0)  theta_z_y = 0;
                        W_axis_x_y    = fabs(theta_x_y - 0.25*M_PI)/(0.25*M_PI);
                        W_axis_x_z    = fabs(theta_z_x - 0.25*M_PI)/(0.25*M_PI);
                        W_axis_y_z    = fabs(theta_z_y - 0.25*M_PI)/(0.25*M_PI);

                        //printf("%f %f %f %f %f %f\n\n", theta_x_y, theta_z_x, theta_z_y, W_axis_x_y, W_axis_x_z, W_axis_y_z);

                        if (x_t2*x_t2 >= y_t2*y_t2 && x_t2*x_t2 >= z_t2*z_t2)   Px_t2 = 1;
                        else if (y_t2*y_t2 >= z_t2*z_t2)                        Py_t2 = 1;
                        else                                                    Pz_t2 = 1;

                        // Diagonals (more than one can be non-zero)
                        if (x_t2*y_t2 >= 0)                   Pxyp_t2 = 1;
                        if (x_t2*z_t2 >= 0)                   Pxzp_t2 = 1;
                        if (y_t2*z_t2 >= 0)                   Pyzp_t2 = 1;

                        //printf("Points t %f %f %f\n", Px_t, Py_t, Pz_t);
                        //printf("Points t2 %f %f %f\n", Px_t2, Py_t2, Pz_t2);

                        // Assign quadrant weights
                        W_axis            = Px_t2*W_axis_x_y*W_axis_x_z          + Py_t2*W_axis_x_y*W_axis_y_z          + Pz_t2*W_axis_x_z*W_axis_y_z;
                        W_diag_plane        = Px_t2*(1-W_axis_x_y)*(W_axis_x_z)    + Py_t2*(1-W_axis_x_y)*(W_axis_y_z)    + Pz_t2*(1-W_axis_x_z)*W_axis_y_z;
                        W_diag_elevation    = Px_t2*W_axis_x_y*(1-W_axis_x_z)      + Py_t2*W_axis_x_y*(1-W_axis_y_z)      + Pz_t2*(W_axis_x_z)*(1-W_axis_y_z);
                        W_diag_diag         = Px_t2*(1-W_axis_x_y)*(1-W_axis_x_z)  + Py_t2*(1-W_axis_x_y)*(1-W_axis_y_z)  + Pz_t2*(1-W_axis_x_z)*(1-W_axis_y_z);

                        //printf("W %f %f %f %f\n", W_axis, W_diag_plane, W_diag_elevation, W_diag_diag);

                        // Assign weights to directions
                        W_t2_xx                = Px_t2 * W_axis;    // only non-zero if pointing primarily towards x
                        W_t2_yy                = Py_t2 * W_axis;
                        W_t2_zz                = Pz_t2 * W_axis;

                        W_t2_xypp              = (Px_t2 + Py_t2) * Pxyp_t2      * W_diag_plane;
                        W_t2_xypm              = (Px_t2 + Py_t2) * (1-Pxyp_t2)  * W_diag_plane;

                        W_t2_xzpp              = Px_t2 * Pxzp_t2     * W_diag_elevation      + Pz_t2 * Pxzp_t2     * W_diag_plane;
                        W_t2_xzpm              = Px_t2 * (1-Pxzp_t2) * W_diag_elevation      + Pz_t2 * (1-Pxzp_t2) * W_diag_plane;

                        W_t2_yzpp              = (Py_t2 + Pz_t2) * Pyzp_t2      * W_diag_elevation;
                        W_t2_yzpm              = (Py_t2 + Pz_t2) * (1-Pyzp_t2)  * W_diag_elevation;

                        W_t2_xyzppp            = Pxzp_t2      *   Pyzp_t2    * W_diag_diag;
                        W_t2_xyzpmp            = (1-Pxyp_t2)  *   Pxzp_t2    * W_diag_diag;
                        W_t2_xyzmpp            = (1-Pxyp_t2)  *  (1-Pxzp_t2) * W_diag_diag;
                        W_t2_xyzppm            = Pxyp_t2      *  (1-Pxzp_t2) * W_diag_diag;
                        // End transverse 2 =============================================//|
                        // End Transverse directions =====================================================//|
                    } // end checking if fibre exists
                    else // if no fibre
                    {
                        printf("ERROR: No fibre found for node %d %d %d || model will work but is desgined for anisotropic simulations\n", i, j, k);
                    }

                    // Set all below to 1.0 if just want standard implementation; these factors just help ensure symmetries
                    double symm_fac_diag_long           = 1.0;
                    double symm_fac_diag_trans          = 1.0;
                    double symm_fac_corner              = 1.0;

                    // Assign gGap based on weights
                    gGap_node_xx[i][j][k]        = Gt*W_t_xx     + Gt2*W_t2_xx     + Gl*W_xx; // simply = weight in transverse * gt + weight along fibre * Gl
                    gGap_node_yy[i][j][k]        = Gt*W_t_yy     + Gt2*W_t2_yy     + Gl*W_yy; // note that both weights can be 0, giving no contribution to said direction
                    gGap_node_zz[i][j][k]        = Gt*W_t_zz     + Gt2*W_t2_zz     + Gl*W_zz;

                    gGap_node_xypp[i][j][k]      = (Gt*W_t_xypp   + Gt2*W_t2_xypp)*symm_fac_diag_trans   + Gl*W_xypp * symm_fac_diag_long;
                    gGap_node_xypm[i][j][k]      = (Gt*W_t_xypm   + Gt2*W_t2_xypm)*symm_fac_diag_trans   + Gl*W_xypm * symm_fac_diag_long;
                    gGap_node_xzpp[i][j][k]      = (Gt*W_t_xzpp   + Gt2*W_t2_xzpp)*symm_fac_diag_trans   + Gl*W_xzpp * symm_fac_diag_long;
                    gGap_node_xzpm[i][j][k]      = (Gt*W_t_xzpm   + Gt2*W_t2_xzpm)*symm_fac_diag_trans   + Gl*W_xzpm * symm_fac_diag_long;
                    gGap_node_yzpp[i][j][k]      = (Gt*W_t_yzpp   + Gt2*W_t2_yzpp)*symm_fac_diag_trans   + Gl*W_yzpp * symm_fac_diag_long;
                    gGap_node_yzpm[i][j][k]      = (Gt*W_t_yzpm   + Gt2*W_t2_yzpm)*symm_fac_diag_trans   + Gl*W_yzpm * symm_fac_diag_long;

                    gGap_node_xyzppp[i][j][k]    = Gt*W_t_xyzppp + Gt2*W_t2_xyzppp + Gl*W_xyzppp * symm_fac_corner;
                    gGap_node_xyzppm[i][j][k]    = Gt*W_t_xyzppm + Gt2*W_t2_xyzppm + Gl*W_xyzppm * symm_fac_corner;
                    gGap_node_xyzpmp[i][j][k]    = Gt*W_t_xyzpmp + Gt2*W_t2_xyzpmp + Gl*W_xyzpmp * symm_fac_corner;
                    gGap_node_xyzmpp[i][j][k]    = Gt*W_t_xyzmpp + Gt2*W_t2_xyzmpp + Gl*W_xyzmpp * symm_fac_corner;

                    // Geometry scaling
                    // if you want this indepnendet of dx, then remove all dx, dy, dz below, and include 1/sqrt(2) and 1/sqrt(3) factors in diagonals and corners
                    gGap_node_xx[i][j][k]     *= 1.0/dx;
                    gGap_node_yy[i][j][k]     *= 1.0/dy;
                    gGap_node_zz[i][j][k]     *= 1.0/dz;

                    gGap_node_xypp[i][j][k]   *= 1.0/sqrt(dx*dx + dy*dy);
                    gGap_node_xypm[i][j][k]   *= 1.0/sqrt(dx*dx + dy*dy);
                    gGap_node_xzpp[i][j][k]   *= 1.0/sqrt(dx*dx + dz*dz);
                    gGap_node_xzpm[i][j][k]   *= 1.0/sqrt(dx*dx + dz*dz);
                    gGap_node_yzpp[i][j][k]   *= 1.0/sqrt(dy*dy + dz*dz);
                    gGap_node_yzpm[i][j][k]   *= 1.0/sqrt(dy*dy + dz*dz);

                    gGap_node_xyzppp[i][j][k] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);
                    gGap_node_xyzppm[i][j][k] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);
                    gGap_node_xyzpmp[i][j][k] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);
                    gGap_node_xyzmpp[i][j][k] *= 1.0/sqrt(dx*dx + dy*dy + dz*dz);

                } // end if geo is > 1 loop
            } // end x
        } //end y
    } // end spatial loop
    // End Setup gGap for each node and direction based on orientation ============================//|

    // Calc Njunctions =============================================================================\\|
    Njunc = 0;
    for (int k = 0; k < NZ; k++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                if (geo[i][j][k] > 0) // if it is an actual cell/node
                {
                    if (i < NX-1 &&                                         geo[i+1][j][k] > 0)         Njunc++;
                    if (j < NY-1 &&                                         geo[i][j+1][k] > 0)         Njunc++;
                    if (k < NZ-1 &&                                         geo[i][j][k+1] > 0)         Njunc++;
                    if (i < NX-1 && j < NY-1 &&                             geo[i+1][j+1][k] > 0)       Njunc++;
                    if (i > 0        && j < NY-1 &&                         geo[i-1][j+1][k] > 0)       Njunc++;
                    if (i < NX-1 && k < NZ-1 &&                             geo[i+1][j][k+1] > 0)       Njunc++;
                    if (i > 0        && k < NZ-1 &&                         geo[i-1][j][k+1] > 0)       Njunc++;
                    if (j < NY-1 && k < NZ-1 &&                             geo[i][j+1][k+1] > 0)       Njunc++;
                    if (j > 0        && k < NZ-1 &&                         geo[i][j-1][k+1] > 0)       Njunc++;
                    if (i < NX-1 && j < NY-1 && k < NZ-1 &&                 geo[i+1][j+1][k+1] > 0)     Njunc++;
                    if (i > 0        && j < NY-1 && k < NZ-1 &&             geo[i-1][j+1][k+1] > 0)     Njunc++;
                    if (i < NX-1 && j > 0        && k < NZ-1 &&             geo[i+1][j-1][k+1] > 0)     Njunc++;
                    if (i > 0        && j > 0        && k < NZ-1 &&         geo[i-1][j-1][k+1] > 0)     Njunc++; // ppm is same as mmp
                } // end if
            } // end x
        } // end y
    } // end spatial loop
    printf("Number of junctions = %d\n", Njunc);
    // End calc Njunctions =========================================================================//|

    // Allocate arrays which are of size Njunc =====================================================\\|
    double  IDiff;     // junction current (not needed as an array in this implementation, but can be if desired)
    double  *gGap_jn;  // junction conductance
    int     *jn_map_minus_x;  // returns Ncell array index of the minus cell for IDiff[i][j][k]
    int     *jn_map_minus_y;  // returns Ncell array index of the minus cell for IDiff[i][j][k]
    int     *jn_map_minus_z;  // returns Ncell array index of the minus cell for IDiff[i][j][k]
    int     *jn_map_plus_x;
    int     *jn_map_plus_y;
    int     *jn_map_plus_z;
    int     *connection_type_jn; // whether the jucntion is comprised transverse, long, or a mix

    gGap_jn               = new double [Njunc];
    //IDiff               = new double [Njunc];
    jn_map_minus_x        = new int [Njunc];
    jn_map_plus_x         = new int [Njunc];
    jn_map_minus_y        = new int [Njunc];
    jn_map_plus_y         = new int [Njunc];
    jn_map_minus_z        = new int [Njunc];
    jn_map_plus_z         = new int [Njunc];
    connection_type_jn  = new int [Njunc];
    // End Allocate arrays which are of size Njunc =================================================//|

    // Create junctions and relate to nodes ========================================================\\|
    bool xx;
    bool yy;
    bool zz;
    bool xypp;
    bool xypm;
    bool xzpp;
    bool xzpm;
    bool yzpp;
    bool yzpm;
    bool xyzppp;
    bool xyzppm;
    bool xyzpmp;
    bool xyzmpp;

    int count_junc = 0;

    for (int k = 0; k < NZ; k++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int i = 0; i < NX; i++)
            {
                if (geo[i][j][k] > 0) // if it is an actual cell/node
                {
                    xx      = false;
                    yy      = false;
                    zz      = false;
                    xypp    = false;
                    xypm    = false;
                    xzpp    = false;
                    xzpm    = false;
                    yzpp    = false;
                    yzpm    = false;
                    xyzppp  = false;
                    xyzppm  = false;
                    xyzpmp  = false;
                    xyzmpp  = false;

                    if (i < NX-1 &&                                         geo[i+1][j][k] > 0)         xx      = true;
                    if (j < NY-1 &&                                         geo[i][j+1][k] > 0)         yy      = true;
                    if (k < NZ-1 &&                                         geo[i][j][k+1] > 0)         zz      = true;
                    if (i < NX-1 && j < NY-1 &&                             geo[i+1][j+1][k] > 0)       xypp    = true;
                    if (i > 0        && j < NY-1 &&                         geo[i-1][j+1][k] > 0)       xypm    = true;
                    if (i < NX-1 && k < NZ-1 &&                             geo[i+1][j][k+1] > 0)       xzpp    = true;
                    if (i > 0        && k < NZ-1 &&                         geo[i-1][j][k+1] > 0)       xzpm    = true;
                    if (j < NY-1 && k < NZ-1 &&                             geo[i][j+1][k+1] > 0)       yzpp    = true;
                    if (j > 0        && k < NZ-1 &&                         geo[i][j-1][k+1] > 0)       yzpm    = true;
                    if (i < NX-1 && j < NY-1 && k < NZ-1 &&                 geo[i+1][j+1][k+1] > 0)     xyzppp  = true;
                    if (i > 0        && j < NY-1 && k < NZ-1 &&             geo[i-1][j+1][k+1] > 0)     xyzmpp  = true;
                    if (i < NX-1 && j > 0        && k < NZ-1 &&             geo[i+1][j-1][k+1] > 0)     xyzpmp  = true;
                    if (i > 0        && j > 0        && k < NZ-1 &&         geo[i-1][j-1][k+1] > 0)     xyzppm  = true; // ppm is same as mmp

                    // Calculate junction g as an average of the contribution from the two nodes that comprise that junction
                    // Principal directions
                    if (xx == true) // i.e., if there is a connection in the x direction
                    {
                        gGap_jn[count_junc]        = (gGap_node_xx[i][j][k] + gGap_node_xx[i+1][j][k])/2;

                        // "minus" is always current cell
                        jn_map_minus_x[count_junc]  = i; // seems redundant here, but later we loop over junctions so need this
                        jn_map_minus_y[count_junc]  = j; 
                        jn_map_minus_z[count_junc]  = k; 

                        jn_map_plus_x[count_junc]   = i+1;
                        jn_map_plus_y[count_junc]   = j;
                        jn_map_plus_z[count_junc]   = k;

                        count_junc++; // now increments so that next junction is unique
                    } 
                    if (yy == true) // i.e., if there is a connection in the y direction
                    {
                        gGap_jn[count_junc]        = (gGap_node_yy[i][j][k] + gGap_node_yy[i][j+1][k])/2;

                        jn_map_minus_x[count_junc]  = i; 
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i;
                        jn_map_plus_y[count_junc]   = j+1;
                        jn_map_plus_z[count_junc]   = k;

                        count_junc++;
                    }
                    if (zz == true) // i.e., if there is a connection in the z direction
                    {
                        gGap_jn[count_junc]        = (gGap_node_zz[i][j][k] + gGap_node_zz[i][j][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i;
                        jn_map_plus_y[count_junc]   = j;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (xypp == true) 
                    {
                        gGap_jn[count_junc]        = (gGap_node_xypp[i][j][k] + gGap_node_xypp[i+1][j+1][k])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i+1;
                        jn_map_plus_y[count_junc]   = j+1;
                        jn_map_plus_z[count_junc]   = k;

                        count_junc++; 
                    }
                    if (xypm == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_xypm[i][j][k] + gGap_node_xypm[i-1][j+1][k])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i-1;
                        jn_map_plus_y[count_junc]   = j+1;
                        jn_map_plus_z[count_junc]   = k;

                        count_junc++; 
                    }
                    if (xzpp == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_xzpp[i][j][k] + gGap_node_xzpp[i+1][j][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i+1;
                        jn_map_plus_y[count_junc]   = j;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (xzpm == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_xzpm[i][j][k] + gGap_node_xzpm[i-1][j][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i-1;
                        jn_map_plus_y[count_junc]   = j;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (yzpp == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_yzpp[i][j][k] + gGap_node_yzpp[i][j+1][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i;
                        jn_map_plus_y[count_junc]   = j+1;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (yzpm == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_yzpm[i][j][k] + gGap_node_yzpm[i][j-1][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i;
                        jn_map_plus_y[count_junc]   = j-1;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (xyzppp == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_xyzppp[i][j][k] + gGap_node_xyzppp[i+1][j+1][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i+1;
                        jn_map_plus_y[count_junc]   = j+1;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (xyzmpp == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_xyzmpp[i][j][k] + gGap_node_xyzmpp[i-1][j+1][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i-1;
                        jn_map_plus_y[count_junc]   = j+1;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (xyzpmp == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_xyzpmp[i][j][k] + gGap_node_xyzpmp[i+1][j-1][k+1])/2;

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i+1;
                        jn_map_plus_y[count_junc]   = j-1;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                    if (xyzppm == true)
                    {
                        gGap_jn[count_junc]        = (gGap_node_xyzppm[i][j][k] + gGap_node_xyzppm[i-1][j-1][k+1])/2; // ppm is same as mmp

                        jn_map_minus_x[count_junc]  = i;
                        jn_map_minus_y[count_junc]  = j;
                        jn_map_minus_z[count_junc]  = k;

                        jn_map_plus_x[count_junc]   = i-1;
                        jn_map_plus_y[count_junc]   = j-1;
                        jn_map_plus_z[count_junc]   = k+1;

                        count_junc++;
                    }
                } // end if
            } // end x
        } // end y
    } // end spatial loop
    printf("second count of juncs = %d\n", count_junc);
    // End Set g for actual junctions and maps to relate junctions to nodes ========================//|
    // End Setup NETWORK MODEL =============================================================================//|

    // Time loop ===========================================================================================\\|
    for (t = 0; t < total_time; t+=dt)
    {
        // Spatial loop ===============================================================================\\|
        for (int k = 0; k < NZ; k++)
        {
            for (int j = 0; j < NY; j++)
            {
                for (int i = 0; i < NX; i++)
                {
                    // Cell/reaction model ==============================================================\\|
                    // Set Activation gate alpha and beta
                    Ip0d_va_al                 = 0.32*(Vm[i][j][k]+47.13)/(1-exp(-0.1*(Vm[i][j][k]+47.13)));
                    Ip0d_va_bet                = 0.08*exp(-Vm[i][j][k]/11);

                    // Set inactivation gates alphas and betas
                    if (Vm[i][j][k] < -40.0)
                    {
                        Ip0d_vi_1_al           = 0.135*exp((80+Vm[i][j][k])/-6.8);
                        Ip0d_vi_1_bet          = 3.56*exp(0.079*Vm[i][j][k])+310000*exp(0.35*Vm[i][j][k]);
                        Ip0d_vi_2_al           = (-127140*exp(0.2444*Vm[i][j][k])-0.00003474*exp(-0.04391*Vm[i][j][k]))*((Vm[i][j][k]+37.78)/(1+exp(0.311*(Vm[i][j][k]+79.23))));
                        Ip0d_vi_2_bet          = (0.1212*exp(-0.01052*Vm[i][j][k]))/(1+exp(-0.1378*(Vm[i][j][k]+40.14)));
                    }
                    else
                    {
                        Ip0d_vi_1_al           = 0;
                        Ip0d_vi_1_bet          = 1/(0.13*(1+exp((Vm[i][j][k]+10.66)/-11.1)));
                        Ip0d_vi_2_al           = 0;
                        Ip0d_vi_2_bet          = (0.3*exp(-0.0000002535*Vm[i][j][k]))/(1+exp(-0.1*(Vm[i][j][k]+32)));
                    }

                    // Set tau and SS from alpha and beta
                    Ip0d_va_tau                = 1/(Ip0d_va_al + Ip0d_va_bet); // 1/(a+b)
                    Ip0d_vi_1_tau              = 1/(Ip0d_vi_1_al + Ip0d_vi_1_bet);
                    Ip0d_vi_2_tau              = 1/(Ip0d_vi_2_al + Ip0d_vi_2_bet);
                    Ip0d_va_ss                 = Ip0d_va_al * Ip0d_va_tau; // a*tau
                    Ip0d_vi_1_ss               = Ip0d_vi_1_al * Ip0d_vi_1_tau;
                    Ip0d_vi_2_ss               = Ip0d_vi_2_al * Ip0d_vi_2_tau;

                    Ip0d_va[i][j][k]                      = rush_larsen(Ip0d_va[i][j][k], Ip0d_va_ss, Ip0d_va_tau, dt); // lib/Membrane.c
                    Ip0d_vi_1[i][j][k]                    = rush_larsen(Ip0d_vi_1[i][j][k], Ip0d_vi_1_ss, Ip0d_vi_1_tau, dt);
                    Ip0d_vi_2[i][j][k]                    = rush_larsen(Ip0d_vi_2[i][j][k], Ip0d_vi_2_ss, Ip0d_vi_2_tau, dt);

                    Ip0d       = gIp0d * pow(Ip0d_va[i][j][k], 3) * Ip0d_vi_1[i][j][k] * Ip0d_vi_2[i][j][k] * (Vm[i][j][k] - 76);

                    // 3r
                    // Alpha
                    if (fabs(Vm[i][j][k]+14.1)< 1e-10)        Ip3r_va_al  = 0.0015; // Denominator = 0 clause
                    else                                Ip3r_va_al  = 0.0003*(Vm[i][j][k]+14.1)/(1-exp((Vm[i][j][k]+14.1)/-5));

                    // Beta
                    if (fabs(Vm[i][j][k]-3.3328) < 1e-10)     Ip3r_va_bet = 3.7836118e-4;
                    else                                Ip3r_va_bet = 0.000073898*(Vm[i][j][k]-3.3328)/(exp((Vm[i][j][k]-3.3328)/5.1237)-1);

                    // time constant and steady state
                    Ip3r_va_tau        =  1/(Ip3r_va_al+Ip3r_va_bet);
                    Ip3r_va_ss         = sigmoid(Vm[i][j][k], -14.10, -6.5); // V, V1/2, k

                    // Time-independent inactivation gate
                    Ip3r_vi_ti         = sigmoid(Vm[i][j][k], -15, 22.4);

                    Ip3r_va[i][j][k]            = rush_larsen(Ip3r_va[i][j][k], Ip3r_va_ss, Ip3r_va_tau, dt);

                    Ip3r               = gIp3r * Ip3r_va[i][j][k] * Ip3r_vi_ti * (Vm[i][j][k] - (-88));

                    Ip4r               = gIp4r * (Vm[i][j][k] - (-88))/((1+exp(0.07*(Vm[i][j][k] - (-80)))));

                    Itot = Ip0d + Ip3r + Ip4r;
                    // End cell/reaction model ==========================================================//|

                    // Istim
                    if (t < 2.0 && i < 10 && j < 10) Istim = -20;
                    else Istim = 0;

                    // Update voltage
                    Vm_next[i][j][k] = Vm[i][j][k] - dt*(Itot + Istim);

                    // Output single cell file
                    if (iteration_count % (int)(1/dt) == 0)
                    {
                        if (i == 3 && j == 3 && k == 0) fprintf(single_cell_out, "%f %f\n", t, Vm_next[i][j][k]); 
                        if (i == 3 && j == 3 && k == 0) printf("%f %f\n", t, Vm_next[i][j][k]);
                    }
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
            k1 = jn_map_minus_z[n];

            i2 = jn_map_plus_x[n];
            j2 = jn_map_plus_y[n];
            k2 = jn_map_plus_z[n];

            IDiff    = gGap_jn[n] * (Vm[i2][j2][k2] - Vm[i1][j1][k1]);

            // update voltages of nodes connected to that junction
            Vm_next[i2][j2][k2] = Vm_next[i2][j2][k2] + dt*(-IDiff);
            Vm_next[i1][j1][k1] = Vm_next[i1][j1][k1] + dt*(IDiff);
        }
        // End NETWORK MODEL coupling =================================================================//|

        // assign Vm to updated Vm_next
        for (int k = 0; k < NZ; k++)
        {
            for (int j = 0; j < NY; j++)
            {
                for (int i = 0; i < NX; i++)
                {
                    Vm[i][j][k] = Vm_next[i][j][k];
                }
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
            fprintf(vtk_out, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
            fprintf(vtk_out, "SPACING 1 1 1\n");
            fprintf(vtk_out, "ORIGIN 0 0 0\n");
            fprintf(vtk_out, "POINT_DATA %d\n", NX*NY*NZ);
            fprintf(vtk_out, "SCALARS geometry float 1\n");
            fprintf(vtk_out, "LOOKUP_TABLE default\n");

            for (int k = 0; k < NZ; k++)
            {
                for (int j = 0; j < NY; j++)
                {
                    for (int i = 0; i < NX; i++)
                    {
                        fprintf(vtk_out, "%f ", Vm[i][j][k]);
                    }
                    fprintf(vtk_out, "\n");
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
    delete []   jn_map_minus_z;
    delete []   jn_map_plus_z;
    delete []   connection_type_jn;

} // end main
