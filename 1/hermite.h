#ifndef HERMITE_H
#define HERMITE_H


class Hermite
{
public:
    Hermite(int);
    Hermite();
    void set_degree(int);
    double evaluate(double x);
private:
    int degree;
};

#endif // HERMITE_H
