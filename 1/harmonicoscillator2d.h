#ifndef HARMONICOSCILLATOR2D_H
#define HARMONICOSCILLATOR2D_H


class HarmonicOscillator2D
{
public:
    HarmonicOscillator2D();
    int getindex(int nx, int ny, int spin);
    void getquantumnumbers(int i);


private:
    int levels = 10;
    int spin = 0;
    int nx = 0;
    int ny = 0;

};

#endif // HARMONICOSCILLATOR2D_H
