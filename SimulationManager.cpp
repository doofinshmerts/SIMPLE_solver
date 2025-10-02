
#include "SimulationManager.h" // replace this with you simulation driver
#include "GraphicObject.h" // replace this with your data structure for storing plot results

#include <omp.h>

#include <stdio.h>
SimulationManager::SimulationManager()
{

    // calculate the color constants
    red_const = 1.0f / (1.0f - RED_1);
    blue_const = 4.0f / (BLUE_1 * BLUE_1 - 2.0f * BLUE_1 * BLUE_2 + BLUE_2 * BLUE_2);
    green_const = 1.0f / GREEN_1;

    // pre-allocate the mesh
    mesh.resize(MESH_X);
    for(int i = 0; i < MESH_X; i++)
    {
        mesh[i].resize(MESH_Y);
    }

    // assign the position and velocity of each cell
    for(int x_i = 0; x_i < MESH_X; x_i++)
    {
        for(int y_i = 0; y_i < MESH_Y; y_i++)
        {
            Node temp;
            temp.x_pos = H*0.5 + double(x_i)*H;
            temp.y_pos = H*0.5 + double(y_i)*H;
            temp.u = 0.0;
            temp.u_ = 0.0;
            temp.u_m = 0.0;
            temp.u__ = 0.0;

            temp.v = 0.0;
            temp.v_ = 0.0;
            temp.v_m = 0.0;
            temp.v__ = 0.0;
            
            temp.p_m = 0.0;
            temp.p = 0.0;
            temp.type = CellType::Normal;
            temp.rho = rho_o;
            temp.rho_m = rho_o;

            if((x_i > 1) && (x_i < (MESH_X - 2)) && (y_i > 1) && (y_i < (MESH_Y - 2)))
            {
                temp.is_edge = false;
            }
            else
            {
                temp.is_edge = true;
            }
            mesh[x_i][y_i] = temp;
        }
    }

    // assign the top and bottom walls
    for(int x_i = 0; x_i < MESH_X; x_i++)
    {
        mesh[x_i][0].type = CellType::WallS;
        mesh[x_i][MESH_Y-1].type = CellType::WallN;
        mesh[x_i][MESH_Y-1].u = 0.0;
        //mesh[x_i][MESH_Y-1].phi = 1.0;
    }
    // assign the exit
    for (int y_i = 0; y_i < MESH_Y; y_i++)
    {
        mesh[MESH_X - 1][y_i].type = CellType::ExitE;
    }

    // assign the inlet
    for (int y_i = 0; y_i < MESH_Y - 1; y_i++)
    {
        mesh[0][y_i].type = CellType::InletW;
        mesh[0][y_i].u = 10.0;
        mesh[0][y_i].u_m = 0.0;
        mesh[0][y_i].v = 0.0;
        mesh[0][y_i].v_m = 0.0;
        mesh[0][y_i].p = 0.0;
    }

    //mesh[MESH_X *0.5][MESH_Y * 0.5].type = CellType::ExitE;

    
    int y_start = MESH_Y * 0.2;
    int y_end = MESH_Y * 0.25;

    for (int y = y_start; y < y_end; y++)
    {
        //mesh[0][y].phi = 1.0;
    }
    /*
    for (int y = y_start; y < y_end; y++)
    {
        mesh[0][y].type = CellType::InletW;
        mesh[0][y].u = 10.0;
    }
    */

    /*
    double r = H * double(MESH_Y) * 0.01;
    r = r*r;
    double cx = H * double(MESH_X)*0.1;
    double cy = H * double(MESH_Y)*0.2;

    
    for(int x = 0; x < MESH_X; x++)
    {
        for(int y = 0; y < MESH_Y; y++)
        {
            Node* node = &mesh[x][y];

            double del_x = node->x_pos - cx;
            double del_y = node->y_pos - cy;
            double dist = del_x * del_x + del_y * del_y;

            if( dist < r)
            {
                mesh[x][y].phi = 1.0;
                mesh[x][y].is_source = true;
            }
        }
    }

    cx = H * double(MESH_X) * 0.3;
    cy = H * double(MESH_Y) * 0.5;

    for (int x = 0; x < MESH_X; x++)
    {
        for (int y = 0; y < MESH_Y; y++)
        {
            Node *node = &mesh[x][y];

            double del_x = node->x_pos - cx;
            double del_y = node->y_pos - cy;
            double dist = del_x * del_x + del_y * del_y;

            if (dist < r)
            {
                mesh[x][y].phi = 1.0;
                mesh[x][y].is_source = true;
            }
        }
    }

    int x_start = MESH_X * 0.2;
    int x_end = MESH_X * 0.3;

    y_start = MESH_Y * 0.2;
    y_end = MESH_Y * 0.3;

    for(int x = x_start; x < x_end; x++)
    {
        for(int y = y_start; y < y_end; y++)
        {
            mesh[x][y].type = CellType::WallE;
        } 
    }

    

    for(int y = y_start; y < y_end; y++)
    {
        mesh[x_start][y].type = CellType::WallE;
        mesh[x_end][y].type = CellType::WallW;
        //mesh[x_end][y].phi = 1.0;
        mesh[x_start+1][y].type = CellType::WallE;
        mesh[x_end-1][y].type = CellType::WallW;
    }

    for (int x = x_start; x < x_end; x++)
    {
        mesh[x][y_start].type = CellType::WallN;
        mesh[x][y_start + 1].type = CellType::WallN;
        mesh[x][y_end].type = CellType::WallS;
        mesh[x][y_end - 1].type = CellType::WallS;
    }
    */
    return;
}

void SimulationManager::InitializeSimulation(GraphicObject *_mesh, GraphicObject *_vec, GraphicObject *_res)
{
    mesh_plot = _mesh;
    vector_feild = _vec;
    result_plot = _res;

    mesh_plot->vertices.resize(MESH_X * MESH_Y * 8); // lines
    mesh_plot->colors.resize(MESH_X * MESH_Y * 8);

    double d = 1.0 / MESH_X;
    double h_d = 0.5 * d;
    // set the position of the vertices
    for (int xi = 0; xi < MESH_X; xi++)
    {
        double x = (0.5 + (double(xi))) * d;
        int base_idx = xi * MESH_Y * 8;

        for (int yi = 0; yi < MESH_Y; yi++)
        {
            double y = (0.5 + (double(yi))) * d;

            glm::vec3 v1(x - h_d, y + h_d, 0.0);
            glm::vec3 v2(x + h_d, y + h_d, 0.0);
            glm::vec3 v3(x + h_d, y - h_d, 0.0);
            glm::vec3 v4(x - h_d, y - h_d, 0.0);
            // top edge
            mesh_plot->vertices[base_idx + yi * 8] = v1;
            mesh_plot->vertices[base_idx + yi * 8 + 1] = v2;

            // right edge
            mesh_plot->vertices[base_idx + yi * 8 + 2] = v2;
            mesh_plot->vertices[base_idx + yi * 8 + 3] = v3;

            // bottom edge
            mesh_plot->vertices[base_idx + yi * 8 + 4] = v3;
            mesh_plot->vertices[base_idx + yi * 8 + 5] = v4;

            // left edge
            mesh_plot->vertices[base_idx + yi * 8 + 6] = v4;
            mesh_plot->vertices[base_idx + yi * 8 + 7] = v1;

            glm::vec3 color;

            switch (mesh[xi][yi].type)
            {
            case CellType::Normal:
                color = glm::vec3(1.0f, 1.0f, 1.0f);
                break;
            case CellType::WallN:
                // set the pressure at the wall to that of the south boundary
                color = glm::vec3(1.0f, 0.0f, 0.0f);
                break;
            case CellType::WallW:
                // set the pressure at the wall to that of the east boundary
                color = glm::vec3(1.0f, 0.0f, 0.0f);
                break;
            case CellType::WallS:
                // set the pressure at the wall to that of the north boundary
                color = glm::vec3(1.0f, 0.0f, 0.0f);
                break;
            case CellType::WallE:
                // set the pressure at the wall to that of the west boundary
                color = glm::vec3(1.0f, 0.0f, 0.0f);
                break;
            case CellType::SlipN:
                /// set the pressure at the wall to that of the south boundary
                color = glm::vec3(0.0f, 0.0f, 1.0f);
                break;
            case CellType::SlipW:
                // set the pressure at the wall to that of the east boundary
                color = glm::vec3(0.0f, 0.0f, 1.0f);
                break;
            case CellType::SlipS:
                // set the pressure at the wall to that of the north boundary
                color = glm::vec3(0.0f, 0.0f, 1.0f);
                break;
            case CellType::SlipE:
                // set the pressure at the wall to that of the west boundary
                color = glm::vec3(0.0f, 0.0f, 1.0f);
                break;
            case CellType::ExitN:
                // set the pressure at the wall to that of the south boundary
                color = glm::vec3(0.0f, 1.0f, 0.0f);
                break;
            case CellType::ExitW:
                // set the pressure at the wall to that of the east boundary
                color = glm::vec3(0.0f, 1.0f, 0.0f);
                break;
            case CellType::ExitS:
                // set the pressure at the wall to that of the north boundary
                color = glm::vec3(0.0f, 1.0f, 0.0f);
                break;
            case CellType::ExitE:
                // set the pressure at the wall to that of the west boundary
                color = glm::vec3(0.0f, 1.0f, 0.0f);
                break;
            default:
                color = glm::vec3(1.0f, 0.0f, 1.0f);
                break;
            }

            mesh_plot->colors[base_idx + yi * 8] = color;
            mesh_plot->colors[base_idx + yi * 8 + 1] = color;
            mesh_plot->colors[base_idx + yi * 8 + 2] = color;
            mesh_plot->colors[base_idx + yi * 8 + 3] = color;
            mesh_plot->colors[base_idx + yi * 8 + 4] = color;
            mesh_plot->colors[base_idx + yi * 8 + 5] = color;
            mesh_plot->colors[base_idx + yi * 8 + 6] = color;
            mesh_plot->colors[base_idx + yi * 8 + 7] = color;
        }
    }

    mesh_plot->UpdateBuffers();



    vector_feild->vertices.resize(RES_X_VEL * RES_Y_VEL* 2); // lines (individual)
    vector_feild->colors.resize(RES_X_VEL * RES_Y_VEL * 2);

    vector_feild->UpdateBuffers();

    result_plot->vertices.resize(MESH_X * MESH_Y * 6); // triangle mesh
    result_plot->colors.resize(MESH_X * MESH_Y * 6);

    // set the position of the vertices
   

    // set the position of the vertices
    for (int xi = 0; xi < MESH_X; xi++)
    {
        double x = (0.5 + (double(xi))) * d;
        int base_idx = xi * MESH_Y * 6;

        for (int yi = 0; yi < MESH_Y; yi++)
        {
            double y = (0.5 + (double(yi))) * d;

            glm::vec3 v1(x - h_d, y + h_d, 0.0);
            glm::vec3 v2(x + h_d, y + h_d, 0.0);
            glm::vec3 v3(x + h_d, y - h_d, 0.0);
            glm::vec3 v4(x - h_d, y - h_d, 0.0);
            // first triangle
            result_plot->vertices[base_idx + yi * 6] = v4;
            result_plot->vertices[base_idx + yi * 6 + 1] = v3;
            result_plot->vertices[base_idx + yi * 6 + 2] = v2;
            // second triangle
            result_plot->vertices[base_idx + yi * 6 + 3] = v2;
            result_plot->vertices[base_idx + yi * 6 + 4] = v1;
            result_plot->vertices[base_idx + yi * 6 + 5] = v4;
        }
    }

    result_plot->UpdateBuffers();
    PlotVectorFeild();

    //UpdateOuterIteration();
}

void SimulationManager::UpdateOuterIteration()
{
    // advance the values
    // u__ = u_
    // u_ = u
    // u_m = u

    // v__ = v_
    // v_ = v
    // v_m = v

    // p_m = p
    // rho_m = rho

    time_step++;
    time += delta_time;

    for (int x = 0; x < MESH_X; x++)
    {
        for (int y = 0; y < MESH_Y; y++)
        {
            if (mesh[x][y].type == CellType::Normal)
            {
                mesh[x][y].u__ = mesh[x][y].u_;
                mesh[x][y].u_ = mesh[x][y].u_m * 0.6 + 0.4 * mesh[x][y].u;
                mesh[x][y].u_m = mesh[x][y].u_m * 0.6 + 0.4 * mesh[x][y].u;

                mesh[x][y].v__ = mesh[x][y].v_;
                mesh[x][y].v_ = mesh[x][y].v_m * 0.6 + 0.4 * mesh[x][y].v;
                mesh[x][y].v_m = mesh[x][y].v_m * 0.6 + 0.4 * mesh[x][y].v;

                mesh[x][y].p_m = mesh[x][y].p;
                mesh[x][y].rho_m = mesh[x][y].rho;
            }
        }
    }

    // update inner iteration
    SolveTimeStep();

    max_p = -10000000000.0;
    min_p =  10000000000.0;

    max_v = -10000000000.0;
    min_v =  10000000000.0;

    for (int i = 0; i < MESH_X; i++)
    {
        for (int j = 0; j < MESH_Y; j++)
        {
            if (mesh[i][j].type == CellType::Normal)
            {
                double v = mesh[i][j].u * mesh[i][j].u + mesh[i][j].v * mesh[i][j].v;
                v = std::sqrt(v);

                double p = mesh[i][j].p;
                if (v > max_v)
                {
                    max_v = v;
                }
                else if(v < min_v)
                {
                    min_v = v;
                }

                if (p > max_p)
                {
                    max_p = p;
                }
                else if (p < min_p)
                {
                    min_p = p;
                }
            }
        }
    }
    range_v = max_v - min_v;
    range_p = max_p - min_p;

    double p_start = mesh[0.5*MESH_X][0.5*MESH_Y].p;
    printf("t_count: %d, v_count %d, p_count: %d \n", t_count, v_count, p_count);
    printf("Time step: %d, Time: %10.3f [s], Velocity range: %10.3f - %10.3f [m/s], Pressure range: %10.3f - %10.3f [Pa], max divergence %10.3f [Kg/s]\n",time_step, time, max_v, min_v, max_p, min_p, max_d);

    int x_start = MESH_X * 0.3;
    int x_end = MESH_X * 0.4;

    int y_start = MESH_Y * 0.4;
    int y_end = MESH_Y * 0.6;

    double front = 0.0;
    double back = 0.0;
    for(int y = y_start; y < y_end; y++)
    {
        front += mesh[x_start][y].p;
        back += mesh[x_end][y].p;
    }

    double drag = H * (front - back);
    printf("pressure: %10.3f\n", p_start);
}

void SimulationManager::SolveTimeStep()
{
    t_count = 0;
    double err = ERRT_GOAL*100;

    while(err > ERRT_GOAL)
    {
        err = 0.0;
        AdvanceInnerIteration();

        if(t_count > MAXIT_T)
        {
            break;
        }

        for(int x = 0; x < MESH_X; x++)
        {
            for(int y = 0; y < MESH_Y; y++)
            {
                if(mesh[x][y].type == CellType::Normal)
                {
                    double e = mesh[x][y].u - mesh[x][y].u_m;
                    err += e * e;
                }
            }   
        }

        err = std::sqrt(err / ( (MESH_X-2) * (MESH_Y - 2))) / range_v;
        //printf("count: %5d, error Time Step: %12.10f\n", count, err);
        t_count++;
    }

    FindConcentration();
}

void SimulationManager::AdvanceInnerIteration()
{
    // step 0: update the previous iteration values
    for(int x = 0; x < MESH_X; x++)
    {
        for(int y = 0; y < MESH_Y; y++)
        {
            double v;
            double u;
            switch (mesh[x][y].type)
            {
            case CellType::Normal:
                mesh[x][y].rho_m = mesh[x][y].rho;
                mesh[x][y].u_m = mesh[x][y].u;
                mesh[x][y].v_m = mesh[x][y].v;
                mesh[x][y].p_m = mesh[x][y].p;
                mesh[x][y].p_c = 0.0;
                break;
            case CellType::WallN:
                //set the pressure at the wall to that of the south boundary 
                mesh[x][y].p_m = mesh[x][y-1].p;
                mesh[x][y].p = mesh[x][y-1].p;
                break;
            case CellType::WallW:
                //set the pressure at the wall to that of the east boundary
                mesh[x][y].p_m = mesh[x+1][y].p;
                mesh[x][y].p= mesh[x+1][y].p;
                break;
            case CellType::WallS:
                // set the pressure at the wall to that of the north boundary
                mesh[x][y].p_m = mesh[x][y+1].p;
                mesh[x][y].p = mesh[x][y+1].p;
                break;
            case CellType::WallE:
                // set the pressure at the wall to that of the west boundary
                mesh[x][y].p_m = mesh[x-1][y].p;
                mesh[x][y].p = mesh[x-1][y].p;
                break;
            case CellType::SlipN:
                /// set the pressure at the wall to that of the south boundary
                mesh[x][y].p_m = mesh[x][y - 1].p;
                mesh[x][y].p = mesh[x][y-1].p;
                mesh[x][y].u_m = mesh[x][y-1].u;
                break;
            case CellType::SlipW:
                // set the pressure at the wall to that of the east boundary
                mesh[x][y].p_m = mesh[x + 1][y].p;
                mesh[x][y].p = mesh[x+1][y].p;
                mesh[x][y].v_m = mesh[x + 1][y].v;
                break;
            case CellType::SlipS:
                // set the pressure at the wall to that of the north boundary
                mesh[x][y].p_m = mesh[x][y + 1].p;
                mesh[x][y].p = mesh[x][y+1].p;
                mesh[x][y].u_m = mesh[x][y + 1].u;
                break;
            case CellType::SlipE:
                // set the pressure at the wall to that of the west boundary
                mesh[x][y].p_m = mesh[x - 1][y].p;
                mesh[x][y].p = mesh[x-1][y].p;
                mesh[x][y].v_m = mesh[x - 1][y].v;
                break;
            case CellType::ExitN:
                // set the pressure at the wall to that of the south boundary
                //mesh[x][y].p_m = mesh[x][y - 1].p;
                //mesh[x][y].p = mesh[x][y-1].p;
                mesh[x][y].u_m = mesh[x][y - 1].u;
                v = mesh[x][y-1].v;
                if(v > 0.0)
                {
                    mesh[x][y].v_m = v;
                }
                else
                {
                    mesh[x][y].v_m = 0.0;
                }
                break;
            case CellType::ExitW:
                // set the pressure at the wall to that of the east boundary
                //mesh[x][y].p_m = mesh[x + 1][y].p;
                //mesh[x][y].p = mesh[x + 1][y].p;
                u = mesh[x + 1][y].u;
                if (u < 0.0)
                {
                    mesh[x][y].u_m = u;
                }
                else
                {
                    mesh[x][y].u_m = 0.0;
                }
                mesh[x][y].v_m = mesh[x + 1][y].v;
                break;
            case CellType::ExitS:
                // set the pressure at the wall to that of the north boundary
                //mesh[x][y].p_m = mesh[x][y + 1].p;
                //mesh[x][y].p = mesh[x][y + 1].p;
                mesh[x][y].u_m = mesh[x][y + 1].u;
                v = mesh[x][y + 1].v;
                if (v < 0.0)
                {
                    mesh[x][y].v_m = v;
                }
                else
                {
                    mesh[x][y].v_m = 0.0;
                }
                break;
            case CellType::ExitE:
                // set the pressure at the wall to that of the west boundary
                //mesh[x][y].p_m = mesh[x - 1][y].p;
                //mesh[x][y].p = mesh[x - 1][y].p;
                u = mesh[x - 1][y].u;
                if (u > 0.0)
                {
                    mesh[x][y].u_m = u;
                }
                else
                {
                    mesh[x][y].u_m = 0.0;
                }
                mesh[x][y].v_m = mesh[x - 1][y].v;
                break;
            case CellType::InletN:
                mesh[x][y].p = mesh[x][y-1].p;
                mesh[x][y].p_m = mesh[x][y-1].p;
                break;
            case CellType::InletW:
                mesh[x][y].p = mesh[x+1][y].p;
                mesh[x][y].p_m = mesh[x+1][y].p;
                break;
            case CellType::InletS:
                mesh[x][y].p = mesh[x][y+1].p;
                mesh[x][y].p_m = mesh[x][y+1].p;
                break;
            case CellType::InletE:
                mesh[x][y].p = mesh[x-1][y].p;
                mesh[x][y].p_m = mesh[x-1][y].p;
                break;
            default:
                break;
            }
        }
    }

    // step 1: Discritize for iteration m
    for(int x = 0; x < MESH_X; x++)
    {
        for(int y = 0; y < MESH_Y; y++)
        {
            DiscretizeCellInnerIteration(x, y);
        }
    }


    // step 2: solve for velocity field approximations
    FindVelocityApproximation();

    /*
    // step 3: solve for pressure correction
    for(int x = 0; x < MESH_RES; x++)
    {
        for(int y = 0; y < MESH_RES; y++)
        {
            GetAdvectionEstimate(x, y);
        }
    }
    */
    // step 3: solve for pressure correction
    for (int x = 0; x < MESH_X; x++)
    {
        for (int y = 0; y < MESH_Y; y++)
        {
            DiscretizePressureCorrection(x, y);
        }
    }

    FindPressureCorrection();
    
    // step 4: solve for the velocity correction
    // step 5: use corrections to find next values
    // step 6: update the density

   

    
}

void SimulationManager::FindVelocityApproximation()
{
    v_count = 0;
    // calculate the error
    double err = ERRV_GOAL*100.0;
    
   
    while(ERRV_GOAL < err)
    {
        v_count++;
        err = 0.0;

        if(v_count > MAXIT_V)
        {
            break;
        }
        for (int x = 0; x < MESH_X; x++)
        {
            for (int y = 0; y < MESH_Y; y++)
            {
                VelocityFeildUpdate(x, y);
                if(mesh[x][y].type == CellType::Normal)
                {
                    double e = mesh[x][y].u_err;
                    err += e * e;
                }
            }
        }

        err = std::sqrt(err / ((MESH_X - 2) * (MESH_Y - 2))) / range_v;
        //printf("count: %5d, vel approx error: %12.10f\n", count, err);
    }

   
}

void SimulationManager::VelocityFeildUpdate(int x, int y)
{
    Node* node = &mesh[x][y];
    
    double ODx;
    double ODy;
    double u_new;
    double v_new;

    switch(node->type)
    {
        case CellType::Normal:
            // results of the off diagonal calculation
            ODx = node->An_ux * mesh[x][y + 1].u + node->Aw_ux * mesh[x - 1][y].u + node->As_ux * mesh[x][y - 1].u + node->Ae_ux * mesh[x + 1][y].u;
            ODx += node->Anw_vx * mesh[x - 1][y + 1].v + node->Asw_vx * mesh[x - 1][y - 1].v + node->Ase_vx * mesh[x + 1][y - 1].v + node->Ane_vx * mesh[x + 1][y + 1].v;

            ODy = node->An_vy * mesh[x][y + 1].v + node->Aw_vy * mesh[x - 1][y].v + node->As_vy * mesh[x][y - 1].v + node->Ae_vy * mesh[x + 1][y].v;
            ODy += node->Anw_uy * mesh[x - 1][y + 1].u + node->Asw_uy * mesh[x - 1][y - 1].u + node->Ase_uy * mesh[x + 1][y - 1].u + node->Ane_uy * mesh[x + 1][y + 1].u;
            u_new = (-node->Q_x - ODx) / node->Ap_x;
            v_new = (-node->Q_y - ODy) / node->Ap_y;
            u_new = u_new - node->u;
            v_new = v_new - node->v;
            node->u = node->u + w_v * u_new;
            node->v = node->v + w_v * v_new;
            node->u_err = u_new;
            node->v_err = v_new;
            break;

        case CellType::SlipN:
            // set the u velocity to that of the south node
            node->u = mesh[x][y-1].u;
            break;
            
        case CellType::SlipW:
            // set the v velocity to that of the east node
            node->v = mesh[x+1][y].v;
            break;

        case CellType::SlipS:
            // set the u velocity to that of the north node
            node->u = mesh[x][y+1].u;
            break;

        case CellType::SlipE:
            // set the v velocity to that of the west node
            node->v = mesh[x - 1][y].v;
            break;

        case CellType::ExitN:
            // set the u and v velocity to that of the south node
            node->u = mesh[x][y-1].u;
            node->v = mesh[x][y-1].v;
            break;

        case CellType::ExitW:
            // set the u and v velocity to that of the east node
            node->u = mesh[x-1][y].u;
            node->v = mesh[x-1][y].v;
            break;

        case CellType::ExitS:
            // set the u and v velocity to that of the north node
            node->u = mesh[x][y+1].u;
            node->v = mesh[x][y+1].v;
            break;

        case CellType::ExitE:

            // set the u and v velocity to that of the west node
            node->u = mesh[x-1][y].u;
            node->v = mesh[x-1][y].v;
            break;

        default:
    return;
    }
}

void SimulationManager::FindPressureCorrection()
{
    // calculate the error
    double err = ERRP_GOAL * 100.0;
    p_count = 0;

    while (err > ERRP_GOAL)
    {
        p_count++;
        if(p_count > MAXIT_P)
        {
            break;
        }
        err = 0.0;
        for (int x = 0; x < MESH_X; x++)
        {
            for (int y = 0; y < MESH_Y; y++)
            {
                PressureCorrectionUpdate(x, y);
                if(mesh[x][y].type == CellType::Normal)
                {
                    double e = mesh[x][y].p_err;
                    err += e * e;
                }
            }
        }

        err = std::sqrt(err / ( (MESH_X - 2) * (MESH_Y) - 2)) / range_p;
        printf("count: %5d, p_c error: %12.10f\n", p_count, err);
    }

   
    // apply the correction
    for (int x = 0; x < MESH_X; x++)
    {
        for (int y = 0; y < MESH_Y; y++)
        {
            Node* node = &mesh[x][y];

            if(node->type == CellType::Normal)
            {
                node->u_c = -0.5 * H * (mesh[x+1][y].p_c - mesh[x-1][y].p_c) / node->Ap_x;
                node->v_c = -0.5 * H * (mesh[x][y+1].p_c - mesh[x][y-1].p_c) / node->Ap_y;
                node->p = node->p_m + a * node->p_c;
                node->v = node->v + a * node->v_c;
                node->u = node->u + a * node->u_c; 
            }
        }
    }
}

void SimulationManager::FindConcentration()
{
    for(int x = 0; x < MESH_X; x++)
    {
        for(int y = 0; y < MESH_Y; y++)
        {
            DiscretizeConcentration(x, y);
        }
    }

    // calculate the error
    double err = ERRI_GOAL * 100.0;
    phi_count = 0;

    while (err > ERRP_GOAL)
    {
        phi_count++;
        if (phi_count > MAXIT_PHI)
        {
            break;
        }
        err = 0.0;
        for (int x = 0; x < MESH_X; x++)
        {
            for (int y = 0; y < MESH_Y; y++)
            {
                ConcentrationUpdate(x, y);

                if(mesh[x][y].type == CellType::Normal)
                {
                    double e = mesh[x][y].phi_err;
                    err += e*e;
                }
            }
        }

        err = std::sqrt(err / ((MESH_X - 2) * (MESH_Y)-2));
        //printf("count: %5d, p_c error: %12.10f\n", p_count, err);
    }
}

void SimulationManager::ConcentrationUpdate(int x, int y)
{
    Node *node = &mesh[x][y];

    double ODx;
    double phi_new;

    if(node->type == CellType::Normal && !node->is_source)
    {
        // results of the off diagonal calculation
        ODx = node->An_ux * mesh[x][y + 1].phi + node->Aw_ux * mesh[x - 1][y].phi + node->As_ux * mesh[x][y - 1].phi + node->Ae_ux * mesh[x + 1][y].phi;

        phi_new = ODx * node->Ap_x;
        phi_new = phi_new - node->phi;
        node->phi_err = phi_new;
        node->phi = node->phi + w_phi * phi_new;
    }
}

void SimulationManager::PressureCorrectionUpdate(int x, int y)
{
    Node *node = &mesh[x][y];

    double ODx;
    double pc_new;

    switch (node->type)
    {
    case CellType::Normal:
        // results of the off diagonal calculation
        ODx = node->Apn * mesh[x][y + 1].p_c + node->Apw * mesh[x - 1][y].p_c + node->Aps * mesh[x][y - 1].p_c + node->Ape * mesh[x + 1][y].p_c;
        if(!node->is_edge)
        {
            ODx += node->Apnn * mesh[x][y+2].p_c + node->Apww * mesh[x-2][y].p_c + node->Apss * mesh[x][y-2].p_c + node->Apee * mesh[x+2][y].p_c;
        }
        pc_new = (node->Qp + ODx) / node->Ap_p;
        pc_new = pc_new - node->p_c;
        node->p_err = pc_new;
        node->p_c = node->p_c + w_p * pc_new;
        break;

    case CellType::WallN:
        // set the pressure at the wall to that of the south boundary
        mesh[x][y].p_c = mesh[x][y - 1].p_c;
        break;
    case CellType::WallW:
        // set the pressure at the wall to that of the east boundary
        mesh[x][y].p_c = mesh[x + 1][y].p_c;
        break;
    case CellType::WallS:
        // set the pressure at the wall to that of the north boundary
        mesh[x][y].p_c = mesh[x][y + 1].p_c;
        break;
    case CellType::WallE:
        // set the pressure at the wall to that of the west boundary
        mesh[x][y].p_c = mesh[x - 1][y].p_c;
        break;
    case CellType::SlipN:
        /// set the pressure at the wall to that of the south boundary
        mesh[x][y].p_c = mesh[x][y - 1].p_c;
        break;
    case CellType::SlipW:
        // set the pressure at the wall to that of the east boundary
        mesh[x][y].p_c = mesh[x + 1][y].p_c;
        break;
    case CellType::SlipS:
        // set the pressure at the wall to that of the north boundary
        mesh[x][y].p_c = mesh[x][y + 1].p_c;
        break;
    case CellType::SlipE:
        // set the pressure at the wall to that of the west boundary
        mesh[x][y].p_c = mesh[x - 1][y].p_c;
        break;
    case CellType::ExitN:
        // set the pressure at the wall to that of the south boundary
        break;
    case CellType::ExitW:
        // set the pressure at the wall to that of the east boundary
        break;
    case CellType::ExitS:
        // set the pressure at the wall to that of the north boundary
        break;
    case CellType::ExitE:
        // set the pressure at the wall to that of the west boundary
        break;
    case CellType::InletN:
        mesh[x][y].p_c = mesh[x][y-1].p_c;
        break;
    case CellType::InletW:
        mesh[x][y].p_c = mesh[x+1][y].p_c;
        break;
    case CellType::InletS:
        mesh[x][y].p_c = mesh[x][y+1].p_c;
        break;
    case CellType::InletE:
        mesh[x][y].p_c = mesh[x-1][y].p_c;
        break;
    default:
        break;
    }
}

void SimulationManager::DiscretizeCellInnerIteration(int x, int y)
{

    // A * V + Q = 0

    // get a reference to the current node
    Node* node = &mesh[x][y];

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

    double Ap_y = 0.0;  // the diagonal coefficient
    double An_vy = 0.0; // the coefficient of the north velocity
    double Anw_uy = 0.0;
    double Aw_vy = 0.0;
    double Asw_uy = 0.0;
    double As_vy = 0.0;
    double Ase_uy = 0.0;
    double Ae_vy = 0.0;
    double Ane_uy = 0.0;
    double Q_y = 0.0; // the y component of the constant vector

    // contributions from time ########################################################################################################3

    // x-direction
    Ap_x = 3.0 * V * rho_o / (2.0 * delta_time);                          // 3 rho *delV / 2 * delT
    Q_x = V * rho_o * (-4.0 * node->u_ + node->u__) / (2.0 * delta_time); // (rho * delV / 2 * delT) * (-4u_n + u_(n-1))

    // y-direction
    Ap_y = 3.0 * V * rho_o / (2.0 * delta_time);
    Q_y = V * rho_o * (-4.0 * node->v_ + node->v__) / (2.0 * delta_time);
    // contributions from momentum ######################################################################################################

    // only normal cells need to be discritized
    if(node->type != CellType::Normal)
    {
        node->Ap_x = Ap_x;
        node->Ap_y = Ap_y;
        return;
    }

    // use the previous iteration for estimating the flow rates
    double mn = H* 0.5 * ( node->v_m * node->rho_m + mesh[x][y+1].v_m * mesh[x][y+1].rho_m); // 1/2 * ( (rho v)p + (rho v)n)_(m-1)
    double mw = -0.5 * H * ( node->u_m * node->rho_m + mesh[x -1][y].u_m * mesh[x- 1][y].rho_m); // 1/2 * ( (rho v)p + (rho v)n)_(m-1)
    double ms = -0.5 * H * ( node->v_m * node->rho_m + mesh[x][y-1].v_m * mesh[x][y-1].rho_m); // 1/2 * ( (rho v)p + (rho v)n)_(m-1)
    double me = 0.5 * H * (node->u_m * node->rho_m + mesh[x+1][y].u_m * mesh[x+1][y].rho_m); // 1/2 * ( (rho v)p + (rho v)n)_(m-1)

    //printf("mn: %6.3f, mw: %6.3f, ms: %6.3f, me: %6.3f \n", mn, mw, ms, me);
    if(mn > 0)
    {
        // x direction
        Ap_x += mn;
        Q_x += 0.5 * mn * (mesh[x][y+1].u_m - node->u_m);

        // y direction
        Ap_y += mn;
        Q_y += 0.5 * mn * (mesh[x][y+1].v_m - node->v_m);
    }
    else
    {
        // x direction
        An_ux += mn;
        Q_x += 0.5 * mn * (node->u_m - mesh[x][y+1].u_m);

        // y direction
        An_vy += mn;
        Q_y += 0.5 * mn * (node->v_m - mesh[x][y+1].v_m);
    }

    if (mw > 0)
    {
        // x direction
        Ap_x += mw;
        Q_x += 0.5 * mw * (mesh[x-1][y].u_m - node->u_m);

        // y direction
        Ap_y += mw;
        Q_y += 0.5 * mw * (mesh[x-1][y].v_m - node->v_m);
    }
    else
    {
        // x direction
        Aw_ux += mw;
        Q_x += 0.5 * mw * (node->u_m - mesh[x-1][y].u_m);

        // y direction
        Aw_vy += mw;
        Q_y += 0.5 * mw * (node->v_m - mesh[x-1][y].v_m);
    }

    if (ms > 0)
    {
        // x direction
        Ap_x += ms;
        Q_x += 0.5 * ms * (mesh[x][y -1].u_m - node->u_m);

        // y direction
        Ap_y += ms;
        Q_y += 0.5 * ms * (mesh[x][y -1].v_m - node->v_m);
    }
    else
    {
        // x direction
        As_ux += ms;
        Q_x += 0.5 * ms * (node->u_m - mesh[x][y - 1].u_m);

        // y direction
        As_vy += ms;
        Q_y += 0.5 * ms * (node->v_m - mesh[x][y - 1].v_m);
    }

    if (me > 0)
    {
        // x direction
        Ap_x += me;
        Q_x += 0.5 * me * (mesh[x+1][y].u_m - node->u_m);

        // y direction
        Ap_y += me;
        Q_y += 0.5 * me * (mesh[x+1][y].v_m - node->v_m);
    }
    else
    {
        // x direction
        Ae_ux += me;
        Q_x += 0.5 * me * (node->u_m - mesh[x+1][y].u_m);

        // y direction
        Ae_vy += me;
        Q_y += 0.5 * me * (node->v_m - mesh[x+1][y].v_m);
    }

    // contributions from the shear forces ###################################################################################################################################

    Ap_x += Mu * 6.0;
    An_ux += -Mu;
    Aw_ux += -2.0 * Mu;
    As_ux += -Mu;
    Ae_ux += -2.0 * Mu;
    Anw_vx += Mu * 0.25;
    Ase_vx += Mu * 0.25;
    Ane_vx += -Mu * 0.25;
    Asw_vx += -Mu * 0.25;

    Ap_y += Mu*6.0;
    An_vy += -2.0 * Mu;
    Aw_vy += -Mu;
    As_vy += -2.0 * Mu;
    Ae_vy += -Mu;
    Anw_uy += Mu * 0.25;
    Ase_uy += Mu * 0.25;
    Ane_uy += -Mu * 0.25;
    Asw_uy += -Mu * 0.25;

    // contributions from the pressure forces ######################################################################################################################################

    Q_x += 0.5 * (mesh[x+1][y].p_m - mesh[x-1][y].p_m) * H;
    Q_y += 0.5 * (mesh[x][y+1].p_m - mesh[x][y-1].p_m) * H;

    // record the constants

    
    //printf("x: %d, y: %d, Ap_x = %6.3f, An_ux = %6.3f, Aw_ux = %6.3f, As_ux = %6.3f, Ae_ux = %6.3f, Anw_vx = %6.3f, Ase_vx = %6.3f, Ane_vx = %6.3f, Asw_vx = %6.3f, Q_x = %6.3f \n",x,y,Ap_x,An_ux,Aw_ux,As_ux,Ae_ux,Anw_vx,Ase_vx,Ane_vx,Asw_vx,Q_x);
    
    node->Ap_x = Ap_x;
    node->An_ux = An_ux;
    node->Aw_ux = Aw_ux;
    node->As_ux = As_ux;
    node->Ae_ux = Ae_ux;
    node->Anw_vx = Anw_vx;
    node->Ase_vx = Ase_vx;
    node->Ane_vx = Ane_vx;
    node->Asw_vx = Asw_vx;
    node->Q_x = Q_x;

    node->Ap_y = Ap_x;
    node->An_vy = An_ux;
    node->Aw_vy = Aw_ux;
    node->As_vy = As_ux;
    node->Ae_vy = Ae_ux;
    node->Anw_uy = Anw_vx;
    node->Ase_uy = Ase_vx;
    node->Ane_uy = Ane_vx;
    node->Asw_uy = Asw_vx;
    node->Q_y = Q_y;
}

void SimulationManager::DiscretizePressureCorrection(int x, int y)
{
    Node* node = &mesh[x][y];

    if(node->type != CellType::Normal)
    {
        return;
    }

    node->Qp = -(mesh[x + 1][y].u - mesh[x - 1][y].u + mesh[x][y + 1].v - mesh[x][y - 1].v) * H_inv;

    // this is used to determine if second order discritization can be used
    if(!node->is_edge)
    {
        double c = 0.1666666667;
        node->Ap_p = 1.0;
        node->Apn = c;
        node->Apw = c;
        node->Aps = c;
        node->Ape = c;
        c = c * 0.5;
        node->Apnn = c;
        node->Apww = c;
        node->Apss = c;
        node->Apee = c;
        node->Qp = node->Qp * 0.33333;
    }
    else
    {

        double c = 0.25;
        node->Ap_p = 1.0;
        node->Apn = c;
        node->Apw = c;
        node->Aps = c;
        node->Ape = c;
        node->Qp = node->Qp * c;
    }
}

void SimulationManager::DiscretizeConcentration(int x, int y)
{
    
    // get a reference to the current node
    Node *node = &mesh[x][y];

    if(node->type != CellType::Normal)
    {
        return;
    }

    // use the previous iteration for estimating the flow rates
    double mn = H * 0.5 * (node->v * node->rho + mesh[x][y + 1].v * mesh[x][y + 1].rho);  // 1/2 * ( (rho v)p + (rho v)n)_(m-1)
    double mw = -0.5 * H * (node->u * node->rho + mesh[x - 1][y].u * mesh[x - 1][y].rho); // 1/2 * ( (rho v)p + (rho v)n)_(m-1)
    double ms = -0.5 * H * (node->v * node->rho + mesh[x][y - 1].v * mesh[x][y - 1].rho); // 1/2 * ( (rho v)p + (rho v)n)_(m-1)
    double me = 0.5 * H * (node->u * node->rho + mesh[x + 1][y].u * mesh[x + 1][y].rho);  // 1/2 * ( (rho v)p + (rho v)n)_(m-1)

    double Ae = 0.5 * me - Gamma;
    double Aw = 0.5 * mw - Gamma;
    double An = 0.5 * mn - Gamma;
    double As = 0.5 * ms - Gamma;

    double Ap = -1.0 / (4.0 * Gamma);

    node->Ap_x = Ap;
    node->Aw_ux = Aw;
    node->Ae_ux = Ae;
    node->An_ux = An;
    node->As_ux = As;
}

void SimulationManager::UpdateSimulation()
{
    //AdvanceInnerIteration();
    UpdateOuterIteration();
    //PlotConcentration();
    //Plotdphidx();
    //Plotdphidy();
    //PlotVectorFeild();
    PlotVectorFeild2();
    PlotVelocity();
    //PlotConcentration();
    //PlotPressure();
    //PlotDivergence();
    return;
}

void SimulationManager::ResetSimulation()
{
    // assign the position and velocity of each cell
    for (int x_i = 0; x_i < MESH_X; x_i++)
    {
        for (int y_i = 0; y_i < MESH_Y; y_i++)
        {
            Node temp;
            temp.x_pos = H * 0.5 + double(x_i) * H;
            temp.y_pos = H * 0.5 + double(y_i) * H;
            temp.u = 0.0;
            temp.u_ = 0.0;
            temp.u_m = 0.0;
            temp.u__ = 0.0;

            temp.v = 0.0;
            temp.v_ = 0.0;
            temp.v_m = 0.0;
            temp.v__ = 0.0;

            temp.p_m = 0.0;
            temp.p = 0.0;
            temp.type = CellType::Normal;
            temp.rho = rho_o;
            temp.rho_m = rho_o;
            mesh[x_i][y_i] = temp;
        }
    }

    // assign the inlet
    for (int y_i = 0; y_i < MESH_Y; y_i++)
    {
        mesh[0][y_i].type = CellType::InletW;
        mesh[0][y_i].u = 10.0;
        mesh[0][y_i].u_m = 10.0;
        mesh[0][y_i].v = 0.0;
        mesh[0][y_i].v_m = 0.0;
        mesh[0][y_i].p = 0.0;
    }

    // assign the exit
    for (int y_i = 0; y_i < MESH_Y; y_i++)
    {
        mesh[MESH_X - 1][y_i].type = CellType::ExitE;
    }

    // assign the top and bottom walls
    for (int x_i = 0; x_i < MESH_X; x_i++)
    {
        mesh[x_i][0].type = CellType::WallS;
        mesh[x_i][MESH_Y - 1].type = CellType::WallN;
    }
    return;
}

glm::vec3 SimulationManager::MapColor(double val)
{
    // normalize the value
    float x = float(val);

    // if it is greater than 1 then make it red
    if (x >= 1.0f)
    {
        return glm::vec3(1.0f, 0.0f, 0.0f);
    }

    // calculate the color values
    float blue = (GREEN_1 - x) * green_const;

    float green = blue_const * (x - BLUE_1) * (BLUE_2 - x);

    float red = red_const * (x - RED_1);

    //printf("x: %5.2f, r: %5.2f, g: %5.2f, b: %5.2f\n",val, red, green, blue);

    return glm::vec3(red, green, blue);
}

double SimulationManager::GetSpeedI(int x, int y)
{
    double u = mesh[x][y].u;
    double v = mesh[x][y].v;
    return sqrt(u*u + v*v);
}

double SimulationManager::GetSpeed(double x, double y)
{
    int x_i = (int)((x * H_inv) - 0.5);
    double x_weight = (H * (0.5 + double(x_i + 1)) - x) * H_inv;

    int y_i = (int)((y * H_inv) - 0.5);
    double y_weight = (H * (0.5 + double(y_i + 1)) - y) * H_inv;

    double v_b = GetSpeedI(x_i, y_i);
    double v_x;
    double v_y;
    
    if(x_i < MESH_X-1)
    {
        v_x = GetSpeedI(x_i + 1, y_i);
    }
    else
    {
        v_x = v_b;
    }

    if(y_i < MESH_Y-1)
    {
        v_y = GetSpeedI(x_i, y_i + 1);
    }
    else
    {
        v_y = v_b;
    }

    double u = x_weight * v_x + y_weight * v_y + (1.0 - x_weight - y_weight) * v_b;
    return u;
}

void SimulationManager::PlotVectorFeild()
{

    float v_f = 1.0f / (max_v);

    double d = 1.0 / RES_X_VEL;
    double h_d = 0.5 * d;

    double map_f = double(MESH_X) / double(RES_X_VEL);

    for (int xi = 0; xi < RES_X_VEL; xi++)
    {
        double x = (0.5 + (double(xi))) * H * map_f;
        double xv = (0.5 + (double(xi))) * d;
        int base_idx = xi*RES_Y_VEL*2;

        for(int yi = 0; yi < RES_Y_VEL; yi++)
        {
            double y = (0.5 + (double(yi))) * H * map_f;
            double yv = (0.5 + (double(yi))) * d;
            glm::vec3 pos = glm::vec3(xv, yv, 0.1);

            glm::vec3 temp = GetVNormal(x,y) * VEC_SCALE + pos;
            vector_feild->vertices[base_idx + yi*2] = pos;
            vector_feild->vertices[base_idx + yi*2 + 1] = temp;

            double u = (GetSpeed(x, y)) * v_f;
            glm::vec3 color = MapColor(u);
            //glm::vec3 color = glm::vec3(1.0, 0.0, 0.0);
            vector_feild->colors[base_idx + yi*2] = color;
            vector_feild->colors[base_idx + yi*2 + 1] = color;
        }
    }

    vector_feild->UpdateBuffers();
}

void SimulationManager::PlotVectorFeild2()
{

    float v_f = 1.0f / (max_v);

    double d = 1.0 / RES_X_VEL;
    double h_d = 0.5 * d;

    double map_f = double(MESH_X) / double(RES_X_VEL);

    for (int xi = 0; xi < RES_X_VEL; xi += 20)
    {
        double x = (0.5 + (double(xi))) * H * map_f;
        double xv = (0.5 + (double(xi))) * d;
        int base_idx = xi * RES_Y_VEL * 2;

        for (int yi = 0; yi < RES_Y_VEL; yi++)
        {
            double y = (0.5 + (double(yi))) * H * map_f;
            double yv = (0.5 + (double(yi))) * d;
            glm::vec3 pos = glm::vec3(xv, yv, 0.1);
            double length = GetSpeed(x,y)*v_f;

            glm::vec3 temp = GetVNormal(x, y)*3.0f * float(length) * VEC_SCALE + pos;
            vector_feild->vertices[base_idx + yi * 2] = pos;
            vector_feild->vertices[base_idx + yi * 2 + 1] = temp;

            //glm::vec3 color = MapColor(length);
            glm::vec3 color = glm::vec3(0.0, 0.0, 0.0);
            vector_feild->colors[base_idx + yi * 2] = color;
            vector_feild->colors[base_idx + yi * 2 + 1] = color;
        }
    }

    vector_feild->UpdateBuffers();
}

void SimulationManager::PlotVelocity()
{
    // set the position of the vertices
    // set the position of the vertices

    double v_f = 1.0 / max_v;
    for (int xi = 0; xi < MESH_X; xi++)
    {
        int base_idx = xi * MESH_Y * 6;

        for (int yi = 0; yi < MESH_Y; yi++)
        {
            glm::vec3 c1 = MapColor(GetSpeedI(xi, yi) * v_f);
            // first triangle
            result_plot->colors[base_idx + yi * 6] = c1;
            result_plot->colors[base_idx + yi * 6 + 1] = c1;
            result_plot->colors[base_idx + yi * 6 + 2] = c1;

            result_plot->colors[base_idx + yi * 6 + 3] = c1;
            result_plot->colors[base_idx + yi * 6 + 4] = c1;
            result_plot->colors[base_idx + yi * 6 + 5] = c1;
        }
    }
    result_plot->UpdateBuffers();
}

void SimulationManager::PlotPressure()
{
    // set the position of the vertices
    // set the position of the vertices
    double p_f = 1.0 / (max_p - min_p);
    for (int xi = 0; xi < MESH_X; xi++)
    {
        int base_idx = xi * MESH_Y * 6;

        for (int yi = 0; yi < MESH_Y; yi++)
        {
            glm::vec3 c1 = MapColor((mesh[xi][yi].p - min_p) * p_f);

            // first triangle
            result_plot->colors[base_idx + yi * 6] = c1;
            result_plot->colors[base_idx + yi * 6 + 1] = c1;
            result_plot->colors[base_idx + yi * 6 + 2] = c1;

            result_plot->colors[base_idx + yi * 6 + 3] = c1;
            result_plot->colors[base_idx + yi * 6 + 4] = c1;
            result_plot->colors[base_idx + yi * 6 + 5] = c1;
        }
    }
    result_plot->UpdateBuffers();
}

void SimulationManager::PlotConcentration()
{
    // set the position of the vertices
    // set the position of the vertices
    for (int xi = 0; xi < MESH_X; xi++)
    {
        int base_idx = xi * MESH_Y * 6;

        for (int yi = 0; yi < MESH_Y; yi++)
        {
            glm::vec3 c1 = MapColor(mesh[xi][yi].phi*CON_SCALE);

            // first triangle
            result_plot->colors[base_idx + yi * 6] = c1;
            result_plot->colors[base_idx + yi * 6 + 1] = c1;
            result_plot->colors[base_idx + yi * 6 + 2] = c1;

            result_plot->colors[base_idx + yi * 6 + 3] = c1;
            result_plot->colors[base_idx + yi * 6 + 4] = c1;
            result_plot->colors[base_idx + yi * 6 + 5] = c1;
        }
    }
    result_plot->UpdateBuffers();
}

void SimulationManager::PlotDivergence()
{
    // set the position of the vertices
    // set the position of the vertices
    for (int xi = 0; xi < MESH_X; xi++)
    {
        int base_idx = xi * MESH_Y * 6;

        for (int yi = 0; yi < MESH_X; yi++)
        {
            glm::vec3 c1 = MapColor(GetDivergenceI(xi, yi) * CON_SCALE);

            // first triangle
            result_plot->colors[base_idx + yi * 6] = c1;
            result_plot->colors[base_idx + yi * 6 + 1] = c1;
            result_plot->colors[base_idx + yi * 6 + 2] = c1;

            result_plot->colors[base_idx + yi * 6 + 3] = c1;
            result_plot->colors[base_idx + yi * 6 + 4] = c1;
            result_plot->colors[base_idx + yi * 6 + 5] = c1;
        }
    }
    result_plot->UpdateBuffers();
}

void SimulationManager::Plotdphidy()
{
    // set the position of the vertices
    // set the position of the vertices
    for (int xi = 0; xi < MESH_X; xi++)
    {
        int base_idx = xi * MESH_Y * 6;

        for (int yi = 0; yi < MESH_Y; yi++)
        {
            glm::vec3 c1 = MapColor(GetdphidyI(xi, yi) * CON_SCALE);
            // first triangle
            result_plot->colors[base_idx + yi * 6] = c1;
            result_plot->colors[base_idx + yi * 6 + 1] = c1;
            result_plot->colors[base_idx + yi * 6 + 2] = c1;

            result_plot->colors[base_idx + yi * 6 + 3] = c1;
            result_plot->colors[base_idx + yi * 6 + 4] = c1;
            result_plot->colors[base_idx + yi * 6 + 5] = c1;
        }
    }
    result_plot->UpdateBuffers();
}

void SimulationManager::Plotdphidx()
{
    // set the position of the vertices
    // set the position of the vertices
    for (int xi = 0; xi < MESH_X; xi++)
    {
        int base_idx = xi * MESH_Y * 6;

        for (int yi = 0; yi < MESH_Y; yi++)
        {
            glm::vec3 c1 = MapColor(GetdphidxI(xi, yi) * CON_SCALE);
            // first triangle
            result_plot->colors[base_idx + yi * 6] = c1;
            result_plot->colors[base_idx + yi * 6 + 1] = c1;
            result_plot->colors[base_idx + yi * 6 + 2] = c1;

            result_plot->colors[base_idx + yi * 6 + 3] = c1;
            result_plot->colors[base_idx + yi * 6 + 4] = c1;
            result_plot->colors[base_idx + yi * 6 + 5] = c1;
        }
    }
    result_plot->UpdateBuffers();
}

void SimulationManager::PlotMesh()
{

}

// get the color of a cell
glm::vec3 SimulationManager::GetBoundaryColor(int x, int y)
{
    return glm::vec3(0.0f, 0.0f, 0.0f);
}

double SimulationManager::GetDivergenceI(int x, int y)
{

    double dvdy = 0;
    double dudx = 0;
    // special case if y == 0 or y == MESH_RES -1
    if (y == 0)
    {
        // us forward differencing
        double phi1 = mesh[x][y].v;
        double phi2 = mesh[x][y + 1].v;
        dvdy = (phi2 - phi1) * H_inv;
    }
    else if (y == (MESH_Y - 1))
    {
        // us backward differencing
        double phi1 = mesh[x][y - 1].v;
        double phi2 = mesh[x][y].v;
        dvdy = (phi2 - phi1) * H_inv;
    }
    else
    {
        // use centered differencing scheme for dphi / dy
        double phi1 = mesh[x][y - 1].v;
        double phi2 = mesh[x][y + 1].v;

        dvdy = (phi2 - phi1) * 0.5 * H_inv;
    }

    // special case if x == 0 or x == MESH_RES -1
    if (x == 0)
    {
        // us forward differencing
        double phi1 = mesh[x][y].u;
        double phi2 = mesh[x + 1][y].u;
        dudx = (phi2 - phi1) * H_inv;
    }
    else if (x == (MESH_X - 1))
    {
        // us backward differencing
        double phi1 = mesh[x - 1][y].u;
        double phi2 = mesh[x][y].u;
        dudx = (phi2 - phi1) * H_inv;
    }
    else
    {
        // use centered differencing scheme for dphi / dy
        double phi1 = mesh[x - 1][y].u;
        double phi2 = mesh[x + 1][y].u;

        dudx = (phi2 - phi1) * 0.5 * H_inv;
    }

    double div = dudx + dvdy;

    if(std::abs(div) > max_d && (x > 1) && (y > 1) && (x < MESH_X - 2) && (y < MESH_Y - 2))
    {
        max_d = std::abs(div);
    }

    return dudx + dvdy;
}

glm::vec3 SimulationManager::GetVNormal(double x, double y)
{
    int x_i = (int)((x * H_inv) - 0.5);
    double x_weight = (H * (0.5 + double(x_i + 1)) - x) * H_inv;

    int y_i = (int)((y * H_inv) - 0.5);
    double y_weight = (H * (0.5 + double(y_i + 1)) - y) * H_inv;

    double u_b = mesh[x_i][y_i].u;
    double v_b = mesh[x_i][y_i].v;
    double u_x;
    double v_x;
    double u_y;
    double v_y;

    if (x_i < MESH_X - 1)
    {
        u_x = mesh[x_i + 1][y_i].u;
        v_x = mesh[x_i + 1][y_i].v;
    }
    else
    {
        u_x = u_b;
        v_x = v_b;
    }

    if (y_i < MESH_Y - 1)
    {
        u_y = mesh[x_i][y_i + 1].u;
        v_y = mesh[x_i][y_i + 1].v;
    }
    else
    {
        v_y = v_b;
        u_y = u_b;
    }

    double u = x_weight * u_x + y_weight * u_y + (1.0 - x_weight - y_weight) * u_b;
    double v = x_weight * v_x + y_weight * v_y + (1.0 - x_weight - y_weight) * v_b;

    glm::vec3 temp = glm::vec3(u, v, 0.0);

    return glm::normalize(temp);
}

glm::vec3 SimulationManager::GetVelocity(double x, double y)
{
    int x_i = (int)((x * H_inv) - 0.5);
    double x_weight = (H * (0.5 + double(x_i + 1)) - x) * H_inv;

    int y_i = (int)((y * H_inv) - 0.5);
    double y_weight = (H * (0.5 + double(y_i + 1)) - y) * H_inv;

    double u_b = mesh[x_i][y_i].u;
    double v_b = mesh[x_i][y_i].v;
    double u_x;
    double v_x;
    double u_y;
    double v_y;

    if (x_i < MESH_X - 1)
    {
        u_x = mesh[x_i + 1][y_i].u;
        v_x = mesh[x_i + 1][y_i].v;
    }
    else
    {
        u_x = u_b;
        v_x = v_b;
    }

    if (y_i < MESH_Y - 1)
    {
        u_y = mesh[x_i][y_i + 1].u;
        v_y = mesh[x_i][y_i + 1].v;
    }
    else
    {
        v_y = v_b;
        u_y = u_b;
    }

    double u = x_weight * u_x + y_weight * u_y + (1.0 - x_weight - y_weight) * u_b;
    double v = x_weight * v_x + y_weight * v_y + (1.0 - x_weight - y_weight) * v_b;

    glm::vec3 temp = glm::vec3(u, v, 0.0);

    return temp;
}

double SimulationManager::GetdphidyI(int x, int y)
{
    // special case if y == 0 or y == MESH_RES -1
    if(y == 0)
    {
        // us forward differencing
        double phi1 = mesh[x][y].phi;
        double phi2 = mesh[x][y+1].phi;
        return (phi2 - phi1) * H_inv;
    }
    else if(y == (MESH_Y - 1))
    {
        // us backward differencing
        double phi1 = mesh[x][y-1].phi;
        double phi2 = mesh[x][y].phi;
        return (phi2 - phi1) * H_inv;
    }

    // use centered differencing scheme for dphi / dy
    double phi1 = mesh[x][y-1].phi;
    double phi2 = mesh[x][y+1].phi;
    
    return (phi2 - phi1) * 0.5 * H_inv;
}

double SimulationManager::GetdphidxI(int x, int y)
{
    // special case if x == 0 or x == MESH_RES -1
    if (x == 0)
    {
        // us forward differencing
        double phi1 = mesh[x][y].phi;
        double phi2 = mesh[x+1][y].phi;
        return (phi2 - phi1) * H_inv;
    }
    else if (x == (MESH_X - 1))
    {
        // us backward differencing
        double phi1 = mesh[x - 1][y].phi;
        double phi2 = mesh[x][y].phi;
        return (phi2 - phi1) * H_inv;
    }

    // use centered differencing scheme for dphi / dy
    double phi1 = mesh[x - 1][y].phi;
    double phi2 = mesh[x + 1][y].phi;

    return (phi2 - phi1) * 0.5 * H_inv;
}
