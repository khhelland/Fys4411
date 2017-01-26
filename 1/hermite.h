#ifndef HERMITE_H
#define HERMITE_H


class Hermite
{
public:
    Hermite(int);
    double evaluate(double x);
private:
    int degree;
};

#endif // HERMITE_H
