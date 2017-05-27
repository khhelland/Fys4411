#include <iostream>
#include <mpi.h>


using namespace std;

int main(int nargs, char* args[])
{
    int rank,nthreads;
    MPI_Init(&nargs,&args);
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    //cout << "Hello from thread " <<rank<<" of "<<nthreads<<endl;
    int ans[2];
    int arr[2];
    arr[0] = rank;
    arr[1] = rank*rank;
    cout<<arr[0]<<arr[1]<<endl;
    MPI_Reduce(&arr,&ans,2,MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);
    if(rank == 0)
    {
        cout<<"sum is "<<ans[0]<<endl;
        cout<<"square sum is "<<ans[1]<<endl;
    }
    MPI_Finalize();
    return 0;
}

