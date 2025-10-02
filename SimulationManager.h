/**
 * Road Simulator
 * @name Holden Holzer
 * @email holden.holzer@ucalgary.ca
 *
 * Modified from CPSC-453 assignment template files
 */

#pragma once

#include <glm/glm.hpp>
#include <vector>

#define BLUE_1 0.05f
#define BLUE_2 0.95f
#define GREEN_1 0.5f
#define RED_1 0.5f

class GraphicObject;

enum CellType
{
    Normal, // normal internal cell
    WallN, // v = 0, u = k, dp/dy = 0
    WallW,
    WallS,
    WallE,
    SlipN,
    SlipW,
    SlipS,
    SlipE,
    ExitN,
    ExitW,
    ExitS,
    ExitE,
    InletN,
    InletW,
    InletS,
    InletE
};

struct Node
{
    double x_pos = 0.0; // the x position of the center of the node
    double y_pos = 0.0; // the y position of the center of the node

    double u__ = 0.0; // u at time n - 1
    double u_  = 0.0; // u at time n
    double u_m = 0.0; // u at iteration m-1
    double u   = 0.0; // u at time n + 1 (current m)
    double u_err = 0.0; // used to measuring convergence
    double u_c = 0.0; // u correction

    double v__ = 0.0; // v at time n - 1
    double v_  = 0.0; // v at time n 
    double v_m = 0.0; // v at iteration m-1
    double v   = 0.0; // v at time n + 1 (current m)
    double v_err = 0.0; // used to measure convergence
    double v_c = 0.0; // v correction

    double p   = 0.0; // pressure at time n + 1 (current m)
    double p_c = 0.0; // pressure correction
    double p_err = 0.0; // used to measure convergence
    double p_m = 0.0; // the pressure at time m - 1
    
    double rho = 0.0; // the density at iteration m
    double rho_m = 0.0; // the density at iteration m - 1

    double phi = 0.0;
    double phi_err = 0.0;

    double m_x = 0.0;
    double m_y = 0.0;


    // coefficients for the x direction
    // format A * V + Q = 0
    double Ap_x = 0.0;  // the diagonal coefficient
    double An_ux = 0.0; // the coefficient of the north velocity
    double Anw_vx = 0.0;
    double Aw_ux = 0.0;
    double Asw_vx = 0.0;
    double As_ux = 0.0;
    double Ase_vx = 0.0;
    double Ae_ux = 0.0;
    double Ane_vx = 0.0;
    double Q_x = 0.0; // the x component of the constant vector

    double Apn = 0.0;
    double Apnn = 0.0;
    double Apw = 0.0;
    double Apww = 0.0;
    double Aps = 0.0;
    double Apss = 0.0;
    double Ape = 0.0;
    double Apee = 0.0;
    double Ap_p = 0.0;
    double Qp = 0.0;

    // coefficients for the y direction
    double Ap_y = 0.0; // the diagonal coefficient
    double An_vy = 0.0; // the coefficient of the north velocity
    double Anw_uy = 0.0;
    double Aw_vy = 0.0;
    double Asw_uy = 0.0;
    double As_vy = 0.0;
    double Ase_uy = 0.0;
    double Ae_vy = 0.0;
    double Ane_uy = 0.0;
    double Q_y = 0.0; // the y component of the constant vector

    CellType type = CellType::Normal;
    bool is_edge = false;
    bool is_source = false;
};

class SimulationManager
{
private: 

    GraphicObject* mesh_plot = nullptr;
    GraphicObject* vector_feild = nullptr;
    GraphicObject* result_plot = nullptr;

    glm::vec3 mesh_color = glm::vec3(0.3f, 0.5f, 0.5f); // the color of a surface
    glm::vec3 vec_color = glm::vec3(1.0f, 0.3f, 0.3f); // the color of a normal cell

    // the constants to be calculated for the color calculations
    float red_const = 0;
    float blue_const = 0;
    float green_const = 0;

    int time_step = 0;
    double time = 0.0;

    const float VEC_SCALE = 0.02;
    const double CON_SCALE = 2.0;
    double Mu = 0.5; // diffusive constant
    const double MAX_U = 10.0;
    const double MAX_P = 10.0;

    double Gamma = 0.3;
    double w_phi = 0.01;

    // STUFF FOR THE SIMULATION
    const int MESH_X = 200;
    const int MESH_Y = 25;

    const int RES_X_VEL = 200; // the resolution in the x direction of the vector plot
    const int RES_Y_VEL = RES_X_VEL * MESH_Y / MESH_X;
    const double rho_o = 1.000; // the density of the fluid
    const double H = 0.2; // the side length of a cell
    const double H_inv = 1.0 / H;
    const double V = H*H; // the volume of a cell
    double delta_time = 0.1; // the size of the time step
    double w_v = 0.01; // the under-relaxation factor for the velocity approximation solver
    double w_p = 0.4; // the under-relaxation factor for the pressure correction solver
    double a = 0.4; // the update factor for the pressure

    const double ERRP_GOAL = 0.0001; // the error goal
    const double ERRV_GOAL = 0.0001;
    const double ERRT_GOAL = 0.0001;
    const double ERRI_GOAL = 0.0001;

    const int MAXIT_T = 20;
    const int MAXIT_P = 20;
    const int MAXIT_V = 20;
    const int MAXIT_PHI = 100;

    int p_count = 0; // the number of iterations taken by the pressure solver
    int v_count = 0; // the number of iterations taken by the velocity solver
    int t_count = 0; // the number of iterations taken by the time solver
    int phi_count = 0;

    double max_p = MAX_P;
    double min_p = 0.0;
    double max_v = MAX_U;
    double min_v = 0.0;
    double max_d = 0.0;
    double range_v = MAX_U;
    double range_p = MAX_P;

    std::vector<std::vector<Node>> mesh; // the mesh that describes the simulation (x,y)
    
public:
    /**
     * constructor for the simulation manager
     */
    SimulationManager();

    /**
     *  initialize the simulation
     */
    void InitializeSimulation(GraphicObject* _mesh, GraphicObject* _vec, GraphicObject* _res);

    /**
     *  update the simulation 
     */
    void UpdateSimulation();

    /**
     * used to reset the simulation
     */
    void ResetSimulation();

private:

    /**
    * updates the time step and increments the velocity and pressure values to the next time step
    * 
    */    
    void UpdateOuterIteration();

    void SolveTimeStep();
    void AdvanceInnerIteration();

    /**
     * discritizes a cell for an inner iteration 
     */
    void DiscretizeCellInnerIteration(int x, int y);
    void DiscretizePressureCorrection(int x, int y);
    void DiscretizeConcentration(int x, int y);

    void FindVelocityApproximation();
    void VelocityFeildUpdate(int x, int y);

    void FindPressureCorrection();
    void PressureCorrectionUpdate(int x, int y);

    void FindConcentration();
    void ConcentrationUpdate(int x, int y);
    

    void PlotVectorFeild();
    void PlotVectorFeild2();
    void PlotVelocity();
    void PlotPressure();
    void PlotConcentration();
    void PlotDivergence();
    void Plotdphidy();
    void Plotdphidx();
    void PlotMesh();
    void PlotBoundaryTypes();

    // get the color of a cell
    glm::vec3 GetBoundaryColor(int x, int y);
   
    /**
     * maps a normalized value (0-1) to a unique color
     */   
    glm::vec3 MapColor(double val);
    double GetDivergenceI(int x, int y);
    glm::vec3 GetVelocity(double x, double y);
    glm::vec3 GetVNormal(double x, double y);
    double GetSpeedI(int x, int y);
    double GetSpeed(double x, double y);
    double GetdphidyI(int x, int y);
    double GetdphidxI(int x, int y);    
};