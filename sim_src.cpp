/*
 *  sim_src.cpp
 *  EIFSims/Users/yonibrowning/Downloads/sim_src.cpp
 *
 *  Created by Yoni Browning on 5/19/13. (modified last on 5/21/13)
 *  Copyright 2013 University of Washington. All rights reserved.
 *
 */

#include "myheader.h"
#include "mt19937ar.h"

using namespace std;
int main(int argc, char *argv[]) {
    int N_ex= 1000;
    int N_inh= 200;
    int N_out= 1;
    double E_ex= 0;
    double E_inh =  -70;
    double theta= -54;
    double V_reset= -60;
    double V_rest= -70;
    double tau_m = 20;
    double tau_ex= 5;
    double tau_inh= 5;
    double g_bar_max= .015;
    double g_bar_inh= .05;
    double tau_plus= 20;
    double tau_minus= 20;
    double A_plus= .005;
    double A_minus= (A_plus*1.05);
    double dt=  .01;
    double t_max= (1000*500);
    double fr_inh= 10;
    double fr_ex0 = 40;
	
    /*for (int t = 0;t<(int)t_max/dt;t++){
        bool*B = poisson(N_ex,100,dt);
        for(int i = 0;i<N_ex;i++){
            if (B[i]){
                cout<<t;cout<<'\n';
            }
        }
    }
    for(int i =0;i<100;i++){
    cout<<(double)random()/RAND_MAX;cout<<'\n';
     }*/
    IF_plastic(N_ex,N_inh,N_out,E_ex,E_inh,theta,V_reset,V_rest,tau_m,tau_ex,tau_inh,g_bar_max,g_bar_inh,tau_plus,tau_minus,A_plus,A_minus,dt,t_max,fr_inh,fr_ex0);
    return 0;
}

//Produce a binary vector of N poisson neurons with a specified firing rate
bool* poisson(int N, int firing_rate, double dt){
    bool* spikes = new bool[N];
    for(int i = 0;i<N;++i){
        spikes[i] = (genrand_real3()<(firing_rate/1000.00*dt));
    }
    
    return spikes;
}

double** new_it_array_2step(int N,double IC){
    double**array = new double*[N];
    for(int i = 0;i<N;i++){
        array[i] = new double[2];
        array[i][0] = 0; array[i][1] = 0;
    }
    return array;
}

int bool_sum(bool* B,int size){
    int count = 0;
    for(int i = 0; i<size;i++){
        if (B[i]){
            count++;
        }
    }
    return count;
}

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
                 )
{   
    
    //Define time
    long long int time_steps = (long long int)ceil(t_max/dt);	// Number of time steps in the simulation.
    cout<<"Hey, this works!"<<endl;
    /*
     //Setup Raster
     double binsize = 5; //time (ms)
     long long int tsteps_per_bin = (long long int)ceil(binsize/dt); //Number of timesteps in a given bin;
     long long int num_bins = (long long int)ceil(time_steps/tsteps_per_bin);
     */
    
    //V stores the "present" and "next" value of the membrane potential of the output cell
    double* V  = new double[2];
    V[0] = V_rest;V[1] = V_rest;
    
    
    double dV, g_inh, g_ex, M;
    //Similar to the structure of V, stores:
    //1) g_inh, the inhibitory conductance of the neuron
    g_inh = 0;
    //2) g_ex, the excititory conductance of the neuron
    g_ex= 0;
    //3) M
    M = 0;
    //4) P
    double* P = new double[N_ex];
    for (int i = 0;i<1;++i){
        P[i] = 0;
    }
    //5) g_bar (cannot use method, as it requires random initialization);
    
    //Read a prexisting text file into gbar array
    double* g_bar = new double[N_ex];
    ifstream steadyState ("gbar_SteadyState30.dat");
    double val;
    if (steadyState.is_open()){
        for (int i = 0;i<N_ex;i++)
        {
            //Save each new line as a new value
            steadyState>>val;
            g_bar[i] = val;
        }
    }else{
        cout<<"Steady state file not available...Starting Over"<<endl;
        for (int i =0; i<N_ex;i++){
            g_bar[i] = genrand_real3()*g_bar_max;
        }
    }

    
    //Determine the activity of the input neurons. This is done according to a poisson process with a given firing rate
    
    int sum_inh;
    bool* spikes_ex;
    bool* spikes_inh;
    bool fires;
    for (int t = 1; t<time_steps; t++){
        spikes_ex = poisson(N_ex,fr_ex0,dt);
        spikes_inh = poisson(N_inh,fr_inh,dt);
        //Modify voltage. If statement handles spiking
        if (V[1]>theta){
            fires =true;
            V[1] = 10;//Arbitrary height of spikes. Strictly cosmetic, for the sake of plotting.
            V[2] = V_reset;
        }else{
            fires = false;
            dV = (1/tau_m)*((V_rest-V[1])+g_ex*(E_ex-V[1])+g_inh*(E_inh-V[1]));
            V[2] = V[1]+dV*dt;
        }
        
        //Adjust conductances
        
        //Inhibitory
        sum_inh = bool_sum(spikes_inh,N_inh);
        if (sum_inh>0){
            g_inh = g_inh+(sum_inh*g_bar_inh);
        }else{
            g_inh = g_inh+(-g_inh/tau_inh)*dt;
        }
        
        //excititory
        if (bool_sum(spikes_ex,N_ex)>0){
            for(int i = 0;i<N_ex;i++){
                if (spikes_ex[i]){
                    g_ex += g_bar[i];
                }
            }
            
        }else{
            g_ex = g_ex+(-g_ex/tau_ex)*dt;
        }
        
        
        //Depression
        if (fires){
            M = M-A_minus;
            for(int i =0;i<N_ex;i++){
                g_bar[i] = min(g_bar_max,(g_bar[i]+P[i]*g_bar_max));
            }
        }else{
            M = M+(-M/tau_minus)*dt;
        }

        
        //Potentiation
        for(int k = 0;k<N_ex;k++){
            if (spikes_ex[k]){
                P[k] = P[k]+A_plus;
                g_bar[k] = max((double)0, g_bar[k]+M*g_bar_max);
            }else{
                P[k]= P[k]+(-P[k]/tau_plus)*dt;
            }
            
        }
        //THIS CODE IS FOR DEBUGGING ONLY
        /*if (bool_sum(spikes_ex,N_ex)>0){
         printf("thats a spike!!");
         }*/
        V[1] = V[2];
        if (t % (int)(1000/dt)==0){
            cout<<(t*(dt/1000));cout<<" second"<<endl;
        }

    }

    /*
     freopen("gbar.txt","w",stdout);
     for(int i =0;i<N_ex;i++){
     cout<<g_bar[i];
     }
     fclose(stdout);
     */
   ofstream outfile;
    outfile.open(("gbar_SteadyState30.dat"));
    
    for(int i = 0; i < N_ex; ++i) {
        outfile << g_bar[i] << endl;
        cout<<g_bar[i]<<endl;
    }
    
    outfile.close();
    
}

