#include <iostream>
#include <cmath>
#include <array>
#include "dco.hpp"
#include <math.h>
#include <chrono>

const double heat_diff_constant = 0.0005;
const double heat_source_term = 0;
const double delta_t = 0.001;
const double delta_x = 0.005;
const int x_span = 1;
const int N_x = (int)(x_span/delta_x)+1;
const int N_p = 2;
const int number_of_time_steps = 1;

// #define N_active 4;

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


// Here starts the function definition

template<typename T, typename PT>
vec_x<T> f(vec_x<T> x, const vec_p<PT> p) {
    
    vec_x<T> y;
    
    for (int i=1; i<N_x-1; i++) {

        y[i] = x[i] +p[0] * (delta_t/pow(delta_x,2)) * (x[i+1] - 2*x[i] + x[i-1]) + p[1] * delta_t;
    }
    
    y[0] = x[0];
    y[N_x-1] = x[N_x - 1];
    
    return y;

}


// passive evolution

template<typename T, typename PT>
vec_x<T> passive_evolution(int n, vec_x<T> x0, const vec_p<PT> p) {
    vec_x<T> x = x0;

    for (int i=0; i<n; i++) {

        x = f(x,p);

    }
    return x;
}


template<typename T, typename PT> 
std::pair<vec_x<T>, jacobian_x_p<T>> adjoint_evolution (int n, vec_x<T> x0, const vec_p<PT> p0) {

    // Augment forward
    
    using DCO_M=typename dco::ga1s<T>; // Adjoint mode
    using DCO_T=typename DCO_M::type; // Adjoint type
    using DCO_TT=typename DCO_M::tape_t; // Tape type
    
    // activate

    using vec_x_dco = vec_x<DCO_T>;
    using vec_p_dco = vec_p<DCO_T>;
    vec_x_dco x_in; vec_p_dco p_in;


    // assign values

    std::copy(x0.begin(), x0.end(), x_in.begin());
    std::copy(p0.begin(), p0.end(), p_in.begin());
    
    // record the augment forward time
    // start recording
    //auto t_forward_start = std::chrono::high_resolution_clock::now();

    // create tape
    DCO_M::global_tape=DCO_TT::create();
    
    jacobian_x_p<T> jac;
    vec_x<T> y;
    
    DCO_M::global_tape->register_variable(x_in.begin(), x_in.end());
    DCO_M::global_tape->register_variable(p_in.begin(), p_in.end());
    
    auto t_passive_start = std::chrono::high_resolution_clock::now();
    vec_x<T> output_test = passive_evolution(n, x0, p0);
    auto t_passive_end = std::chrono::high_resolution_clock::now();
    std::cout<< output_test[0]<< std::endl;
    auto t_forward_start = std::chrono::high_resolution_clock::now();
    vec_x<DCO_T> output = passive_evolution(n, x_in, p_in);
    // end recording
    auto t_forward_end = std::chrono::high_resolution_clock::now();

    // print the recorded time
    std::chrono::duration<double, std::milli> t_duration_augment_forward = t_forward_end-t_forward_start;
    std::chrono::duration<double, std::milli> t_passive = t_passive_end - t_passive_start;
    std::cout << "alpha: ratio value is:" << t_duration_augment_forward/t_passive << std::endl;

    printf(" Augment forward time is: %f ms\n", t_duration_augment_forward);
    printf(" Passive evolution time is: %f ms\n", t_passive);

    // All the values that are NOT given values are default to zero.
    long l_bytes;

    for (int i=0; i<N_x ; i++) {

	// get the values of y
	y[i] = dco::value(output[i]);

	auto t_start_bp = std::chrono::high_resolution_clock::now();
        // Reverse mode
        dco::derivative(output[i]) = 1.0;

	// Inrerpret adjoint
	
	//l_bytes = dco::size_of(dco::ga1s<T>::global_tape);

	DCO_M::global_tape->interpret_adjoint();
	auto t_end_bp = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> t_bp = t_end_bp - t_start_bp;
	printf("The required time for bp is: %f\n", t_bp);
	
	for (int k1=0; k1<N_x;k1++) {
	  (jac.first)[i][k1] = dco::derivative(x_in[k1]);
	}
	for (int k2=0; k2<N_p;k2++) {
	  (jac.second)[i][k2] = dco::derivative(p_in[k2]);
	}

        DCO_M::global_tape->zero_adjoints();
        
        dco::derivative(output[i]) = 0.0;

        

    }
    
    //printf("Tape size in MB: %f\n", l_bytes / 1024.0 / 1024.0);
    DCO_TT::remove(DCO_M::global_tape);
   
    return std::make_pair(y,jac);



}


int main(int argc, char** argv) {
    
    vec_x<double> zeros;
    zeros.fill(0.0);
    vec_x<double> x = zeros;
    x[0] = 2.0; x[N_x-1] = 0.0;
    vec_p<double> p = {heat_diff_constant, heat_source_term};
    
    auto t_start = std::chrono::high_resolution_clock::now();
    auto res = adjoint_evolution(number_of_time_steps,x,p);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t_end - t_start;

    printf(" Time duration is: %f ms\n", ms_double);

    /*
    mat_x_res<double> dydx = (res.second).first;
    mat_p_res<double> dydp = (res.second).second;
    vec_x<double> y = res.first;
    
    // print out the value of y
    printf("y = { ");
  
    for (int i=0; i< N_x; i++) {

         printf("%f ", y[i]);

    }
    printf("} \n ");
    
    //std::cout<< "y = " << '{' << y[0] <<" , "<<y[1]<<'}'<<std::endl;

    // print out the jacobian for x and for p
    
    printf(" dydx = { "); 
    for (int i=0; i< N_x; i++) {
      for (int j=0; j<N_x; j++) {
	
	printf("%f ", dydx[i][j]);
      }
    }
    
    printf(" }\n");
    //std::cout<<" dy0/d {x, p}" << (jac.first)[0][0] <<' '<< (jac.first)[0][1] <<' '<<
    //   (jac.second)[0][0] <<' '<< (jac.second)[0][1] << std::endl;
    
    printf("dydp = { ");
    for (int i=0; i<N_x; i++) {
      for (int j=0; j<N_p; j++) {
	printf(" %f ",dydp[i][j]);
      }
    }
    
    printf("}\n");

    //std::cout<<" dy1/d {x, p}" << (jac.first)[1][0] <<' '<< (jac.first)[1][1] <<' '<<
    //   (jac.second)[1][0] <<' '<< (jac.second)[1][1] << std::endl;
    */

    return 0;
}
