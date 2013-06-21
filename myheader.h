#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#include <cmath>
#include <string>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include "time.h"


using namespace std;

// Function prototypes
double genrand_real3(void); // Real number in (0,1)


// Simple convert to string function
template <typename T>
string to_string(T const& value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}


// Generate a pseudo-normal random number via Box-Muller transformation
void gen_norm(double &n1, double &n2)
{
	double rand1 = genrand_real3();
	double rand2 = genrand_real3();
	
	n1 = (sqrt(-2*log(rand1))*cos(2*M_PI*rand2));
	n2 = (sqrt(-2*log(rand1))*sin(2*M_PI*rand2));
}

int imod(int x, int m) {
    int r = x%m;
    return r<0 ? r+m : r;
}

bool* poisson(int N, int firing_rate, double dt);
double** new_it_array_2step(int N,double IC);
int bool_sum(bool* B,int size);
void IF_plastic (int N_ex,              // Number of excitiory inputs
                                                int N_inh,             // Number of inhibitory inputs
                                                int N_out,             // Number of output cells
                                                double E_ex,			// Leak reversal potential (mV) Excitory Neurons
                                                double E_inh,          // Leak reversal potential (mV) Inhibitory Neurons
                                                double theta,          // Threshold for spiking neurons (mV)
                                                double V_reset,		// Reset voltage after spiking  (mV)
                                                double V_rest,         // Resting voltage of cell (mV)
                                                double tau_m,
                                                double tau_ex,         // Time constant for excitiory synapses (msec)
                                                double tau_inh,        // Time constant for inhibitory synapses (msec)
                                                double g_bar_max,      // Max synaptic conductance (dimless)
                                                double g_bar_inh,      // Step in inhibory conductance (dimless)
                                                double tau_plus,		// Absolute refractory period (ms)
                                                double tau_minus,      // Total simulation time (ms)
                                                double A_plus,         // increment of weight change
                                                double A_minus,         // derement of weight change
                                                double dt,              // timestep
                                                double t_max,          // maximum time
                                                double fr_inh,          // inhibitory firing rate (Hz)
                                                double fr_ex0          // initial excitory firing rate (Hz)
                 );


#endif
