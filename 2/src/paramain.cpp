#include <iostream>
#include <time.h>
#include <mpi.h>
#include <slatervmc.h>


using namespace std;

int main(int nargs, char* args[])
{
    clock_t start,stop;

    start = clock();

    //Set up paralell threads
    int rank,nthreads;
    MPI_Init(&nargs,&args);
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //cout<<"hello from thread"<<rank<<endl;

    //do the thing

    slatervmc sim(12, //particles
                  0.01, //stepsize
                  0.01, //omega
                  0.150,//alpha
                  0.1, //beta
                  rank+1 //seed
                  );
    sim.useJastrow = 0;

    sim.run(1e7, //cycles
            500000 //blocksize
            );

    double  localresults[6];
    double  globalresults[6];

    localresults[0] = sim.kinetic_energy;
    localresults[1] = sim.kineticMeanVar;
    localresults[2] = sim.potential_energy;
    localresults[3] = sim.potentialMeanVar;
    localresults[4] = sim.energy;
    localresults[5] = sim.energyMeanVar;
    //localresults[6] = sim.mean_distance;

    MPI_Reduce(&localresults,&globalresults,6,MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

    if(rank == 0)
    {
        double kinetic = globalresults[0]/nthreads;
        double Kerror = sqrt(globalresults[1])/nthreads;
        double potential = globalresults[2]/nthreads;
        double Verror = sqrt(globalresults[3])/nthreads;
        double energy = globalresults[4]/nthreads;
        double Eerror = sqrt(globalresults[5])/nthreads;
        //double distance = globalresults[6]/nthreads;

        cout<<"K \t sK \t V \t sV \t E \t sE"<<endl;
        cout<<kinetic<<"\t"<<Kerror<<"\t"
            <<potential<<"\t"<<Verror<<"\t"
            <<energy<<"\t"<<Eerror<<endl;

        stop = clock();

        double time = (float)(stop-start)/CLOCKS_PER_SEC;
        cout<<"time [s]: "<<time<<endl;

    }


    MPI_Finalize();


    return 0;
}

