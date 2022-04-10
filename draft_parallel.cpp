#include <iostream>
#include <cmath>
#include <tuple>
#include <mpi.h>
#include "dco.hpp"
#include <stdio.h>

int I, N;

template<typename T, typename PT>
T f(T x, const PT p) { return p*sin(x); }

template<typename T, typename PT>
T passive_evolution(int n, T x0, const PT p) {
    T x=x0;
    for (int i=0; i< n; i++) {
        x=f(x,p);
    }
    return x;
}

template<typename T, typename PT>
std::tuple<T, T, T> adjoint_evolution(int n, T x0, const PT p) {

    //double dx = 1;  //////////
    //double global_dp;
    T x=x0;
    x=passive_evolution(I, x, p);      T y=f(x, p);

    // Augment forward
    using DCO_M=typename dco::ga1s<T>;  // Adjoint mode
    using DCO_T=typename DCO_M::type;  // adjoint type
    using DCO_TT=typename DCO_M::tape_t;  // Tape type
    DCO_T x_in;  dco::value(x_in)=x; // activate
    DCO_T p_in;  dco::value(p_in)=p;

    DCO_M::global_tape=DCO_TT::create();  // create tape
    DCO_M::global_tape->register_variable(x_in);
    DCO_M::global_tape->register_variable(p_in);

    DCO_T output = f(x_in, p_in);

    //T y=output::value(output);
    T dx, transfer_buffer;
    PT dp, global_dp; 
    
    //dco::derivative(output) = 1;  // This part is changed 
 
    if (I==N-1) {

   	    //dco::derivative(output) = 1;
        transfer_buffer = 1.0;

    }
    else {
   
        MPI_Recv(&transfer_buffer, 1, MPI_DOUBLE, I+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
	
    dco::derivative(output) = transfer_buffer;

    DCO_M::global_tape->interpret_adjoint();

    dx = dco::derivative(x_in);
    dp = dco::derivative(p_in);
    
    /*
    dp+=dco::derivative(p_in)*dx;
    dx=dco::derivative(x_in)*dx;
    */

    std::cout<<"Processor: "<< I <<"  dp = " << dp << '\t' << " dx = " << dx <<std::endl;

    if (I>0) {

        MPI_Send(&dx, 1, MPI_DOUBLE, I-1, 0, MPI_COMM_WORLD);

    }
    

    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&dp, &global_dp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    DCO_TT::remove(DCO_M::global_tape);

    return std::tuple<T, T, T>(y, dx, global_dp);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &I); // I is the processor index
    MPI_Comm_size(MPI_COMM_WORLD, &N); // N is the size of the communicator

    double x=1.0, p=1.0;
    
    std::tuple<double, double, double> res= adjoint_evolution(N, x, p);
    double y= std::get<0>(res);
    double dydx = std::get<1>(res);
    double dydp = std::get<2>(res);
    //auto [y, dydx, dydp]=adjoint_evolution(N, x, p);

    if (I==N-1) {
	printf("y = %f\n", y);
        //std::cout<< "y = "<< y << std::endl;
    }

    if (I==0) {
        printf("dydx = %f, \t dydp = %f\n", dydx, dydp);
	/*
	std::cout<< "dydx = "<<dydx<< std::endl;
        std::cout<< "dydp = "<<dydp<< std::endl;
	*/
    }

    MPI_Finalize();
}
