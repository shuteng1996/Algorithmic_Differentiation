#include <iostream>
#include <cmath>
#include <mpi.h>
#include "dco.hpp"
#include <array>
#include <stdio.h>
#include <math.h>
#include <chrono>

int I, N;

// For this task, we assumed that the # of input variables is the same as the
// # of output variables

// Two constant variables here are (1) heat diffusion constant (2) heat source term 
const double heat_diff_constant = 0.0005;
const double heat_source_term = 0;
const double delta_t = 0.001;
const double delta_x = 0.005;
const int x_span = 1;
const int N_x = (int)(x_span/delta_x)+1;
const int N_p = 2;
const int N_jac_x_send = N_x*N_x;
const int N_jac_p_send = N_x*N_p;
const int N_steps_per_processor = 50000;


template<typename T>
using vec_x = std::array<T, N_x>;


template<typename PT>
using vec_p = std::array<PT, N_p>;

template<typename T>
using mat_x_res = std::array<vec_x<T>, N_x>;

template<typename PT>
using mat_p_res = std::array<vec_p<PT>, N_x>;

template<typename T>
using jacobian_x_p = std::pair<mat_x_res<T>, mat_p_res<T>>;


// Here starts the finite element discretization of the heat diffusion equation

template<typename T, typename PT>
vec_x<T> f(vec_x<T> x, const vec_p<PT> p) {
    
    vec_x<T> y;

    for (int i=1; i< N_x-1; i++) {

        y[i] = x[i] + p[0] * (delta_t/pow(delta_x,2)) * (x[i+1] -2*x[i] + x[i-1]) + p[1] * delta_t;

    }

    // Boundary conditions fixed.
    y[0] = x[0];
    y[N_x - 1] = x[N_x - 1];

    return y;

}

// passive evolution for n times

template<typename T, typename PT>
vec_x<T> passive_evolution(int n, vec_x<T> x0, const vec_p<PT> p) {
    vec_x<T> x = x0;

    for (int i=0; i<n; i++) {

        x = f(x,p);

    }

    return x;
}


// From here, it follows the adjoint evolution for Multi-variate vector function
template<typename T, typename PT>
std::pair<vec_x<T>, jacobian_x_p<T>> adjoint_evolution(int N_proc, vec_x<T> x0, const vec_p<PT> p) {
    vec_x<T> x = x0;

    x = passive_evolution(N_steps_per_processor*I, x, p);   vec_x<T> y = passive_evolution(N_steps_per_processor, x, p);

    // Augment forward

    using DCO_M = typename dco::ga1s<T>; //Adjoint mode
    using DCO_T = typename DCO_M::type; //Adjoint type
    using DCO_TT = typename DCO_M::tape_t; //Tape type

    //activate

    using vec_x_dco = vec_x<DCO_T>;
    using vec_p_dco = vec_p<DCO_T>;

    vec_x_dco x_in; vec_p_dco p_in;

    std::copy(x.begin(), x.end(), x_in.begin());
    std::copy(p.begin(), p.end(), p_in.begin());

    // create tape
    DCO_M::global_tape = DCO_TT::create();

    jacobian_x_p<T> jac;
    //vec_x<T> y;

    mat_p_res<PT> p_res;


    if (I == N-1) {

        // contents
        vec_x<T> zeros; zeros.fill(0.0);
        // Step 1: Initialize jac.first
        for (int m = 0; m<N_x; m++) {
            
            jac.first[m] = zeros;
            jac.first[m][m] = 1.0;

        }


    }

    else {

        // contents to fill
        auto res_receive = (jac.first).data();

        MPI_Recv( &res_receive[0][0]  , N_jac_x_send ,  MPI_DOUBLE, I+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        
    }
    

    // The following steps are to calculate the derivative and do the overloading
    DCO_M::global_tape->register_variable(x_in.begin(), x_in.end());
    DCO_M::global_tape->register_variable(p_in.begin(), p_in.end());
    vec_x<DCO_T> output = passive_evolution(N_steps_per_processor, x_in, p_in);
    long l_bytes;
    for (int i=0; i<N_x; i++) {
     
            for (int j=0; j<N_x; j++) {

                dco::derivative(output[j]) = jac.first[i][j];

            }
            

	    l_bytes = dco::size_of(dco::ga1s<T>::global_tape); // tape size in bytes
	    
	    // Interpret adjoint
	    DCO_M::global_tape->interpret_adjoint();

            // update dy/d{p}
            for (int k1 = 0; k1<N_p; k1++) {

                (jac.second)[i][k1] = dco::derivative(p_in[k1]);

            }

            // update dy/d{x}
            for (int k2 = 0; k2<N_x; k2++) {

                (jac.first)[i][k2] = dco::derivative(x_in[k2]);

            }

            //DCO_M::global_tape->reset();
	    DCO_M::global_tape->zero_adjoints();
        }

    printf("Tape size in MB: %f\n", l_bytes / 1024.0 / 1024.0);
    

    if (I>0) {
        // Send to the previous processor

        // Question: How to combine the template with MPI_Type????????
        auto res_send = (jac.first).data();

        MPI_Send(&res_send[0][0], N_jac_x_send , MPI_DOUBLE, I-1, 0, MPI_COMM_WORLD);


    }
    auto res_allreduce = (jac.second).data();
    auto final_result = (p_res).data();

    MPI_Allreduce(&res_allreduce[0][0], &final_result[0][0] ,N_jac_p_send , MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

    return std::make_pair(y, std::make_pair(jac.first, p_res));

}


//


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &I);
    MPI_Comm_size(MPI_COMM_WORLD, &N);
    
    vec_x<double> zeros;
    zeros.fill(0.0);
    vec_x<double> x = zeros;
    x[0] = 2.0;   x[N_x-1] = 0.0; // Fix the boundary condition
    vec_p<double> p = {heat_diff_constant, heat_source_term};


    // Time it.
    auto t_start = std::chrono::high_resolution_clock::now();
    auto res = adjoint_evolution(N, x, p);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t_end - t_start;

    printf(" I am Processor # %d:  My time duration is: %f ms\n", I, ms_double);
    /*
    MPI_Barrier(MPI_COMM_WORLD);
    // obtain the results
    
    vec_x<double> y = res.first;
    mat_x_res<double> dydx = (res.second).first;
    mat_p_res<double> dydp = (res.second).second;
    
    if (I==N-1) {
        
        printf(" y = { ");

        for (int i=0; i<N_x; i++) {

            printf("%f  ",y[i]);

        }

        printf("}\n");

    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (I==0) {
        
        printf(" dydx = { ");

        for (int i=0; i<N_x; i++) {
            for (int j=0; j<N_x; j++) {

                printf(" %f ", dydx[i][j]);

            }
        }

        printf(" }\n");

        //printf("dydx = { %f , %f , %f , %f }\n", dydx[0][0], dydx[0][1], dydx[1][0], dydx[1][1]);
        printf(" dydp = { ");

        for (int i=0; i<N_x; i++) {
            for (int j=0; j<N_p; j++) {

                printf(" %f ", dydp[i][j]);

            }
        }

        printf(" }\n");
        //printf("dydp = { %f , %f , %f , %f }\n", dydp[0][0], dydp[0][1], dydp[1][0], dydp[1][1]);

    }
    */
    MPI_Finalize();
    
    return 0;

}
