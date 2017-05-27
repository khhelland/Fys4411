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

    slatervmc sim(2, //particles
                  0.01, //stepsize
                  1, //omega
                  1,//alpha
                  0.4, //beta
                  rank+1 //seed
                  );
    sim.run(1e7, //cycles
            2e5 //blocksize
            );

    double  localresults[2];
    double  globalresults[2];

    localresults[0] = sim.energy;
    localresults[1] = sim.energyError;

    cout<<"From thread "<<rank<<". Energy: "<<sim.energy<<", error: "<<sim.energyError<<endl;
    MPI_Reduce(&localresults,&globalresults,2,MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

    if(rank == 0)
    {
        double energy = globalresults[0]/nthreads;
        double error = globalresults[1]/pow(nthreads,1.5);
        cout<<"Energy \t Error"<<endl;
        cout<<energy<<"\t"<<error<<endl;
        stop = clock();

        double time = (float)(stop-start)/CLOCKS_PER_SEC;
        cout<<"time [s]: "<<time<<endl;

    }


    MPI_Finalize();


    return 0;
}

