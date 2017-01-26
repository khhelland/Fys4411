#include "harmonicoscillator2d.h"
/*Things the class must be able to do:
    -translate between single index and 3 nx,ny,sz
    -compute Hermite polynomials
    -compute single particle wave functions
*/
HarmonicOscillator2D::HarmonicOscillator2D()
{

}

int HarmonicOscillator2D::getindex(int nx, int ny, int spin)
{
    int i = 2*((nx+ny)*(nx+ny+1)+ny) + (1-spin)/2;
    return i;
}

void HarmonicOscillator2D::getquantumnumbers(int i)
{
    spin = 2(i % 2)-1;
    int j = i/2;
    for(int n = 1, s = 0 ;n < 11; n++)
    {
        s+=n;
        if (s>=j)
        {
            int shell = n;
            ny = j-s;//figure this out

        }
    }

}
