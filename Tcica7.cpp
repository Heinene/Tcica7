#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <random>
#include <stdexcept>

struct zn
{
    double h;
    double S;
    std::vector<double> aa;
    double om;
    double d;
};

double n(const double& P, const double& E, const double& X, const double& Xm)
{
    return (log(1 - P) / log(1 - E / (X - Xm)));
}

std::vector<double> AA(const int& r)
{
    std::vector<double> aa(r);
    std::mt19937 gen(time(0));
    std::uniform_real_distribution<double> urd(0, 1);
    double A = urd(gen);
    aa[(r - 1) / 2] = A;

    for (auto i = (r - 1) / 2 - 1; i > 0; i--)
    {
        double sum = 0;
        for (auto j = i; j < r - i - 1; j++)
            sum += aa[j];
        std::uniform_real_distribution<double> urd1(0, 1 - sum);
        double B = urd1(gen);
        aa[i] = 0.5 * B;
    }
    for (auto i = (r - 1) / 2 + 1; i < r - 1; i++)
        aa[i] = aa[r - i - 1];
    double E = 0;
    for (auto i = 1; i < r - 1; i++)
        E += aa[i];
    aa[0] = 0.5 * (1 - E);
    aa[r - 1] = aa[0];
    return aa;
}

std::vector<double> fil(const std::vector<double>& f, std::vector<double>& aa, const int& r)
{
    std::vector<double> Res(101);
    int M = (r - 1) / 2;
    for (auto i = 0; i <= 100; i++)
    {
        double kva = 0;
        for (auto j = i - M; j <= i + M; j++)
        {
            try
            {
                f.at(j);
                kva += f[j] * f[j] * aa[j + M - i];
            }
            catch (const std::out_of_range& s)
            {
                kva += 0;
            }
        }
        Res[i] = sqrt(kva);
    }
    return Res;
}

double om(const std::vector<double>& F)
{
    double max = 0;
    for (auto i = 1; i <= 100; i++)
        if (max < abs((F[i] - F[i - 1]))) {
            max = abs(F[i] - F[i - 1]);
        };
    return max;
}

double d(const std::vector<double>& F, const std::vector<double>& F_noise)
{
    double max = 0;
    for (auto i = 0; i <= 100; i++)
        if (max < abs((F[i] - F_noise[i]))) {
            max = abs(F[i] - F_noise[i]);
        };
    return max;
}

double S(const double& om, const double& d)
{
    return std::max(om,d);
}

double J(double h, double om, double d)
{
    return h * om + (1 - h) * d;
}

int main()
{
    int N3 = 3, N5 = 5;
    double Xm = 0, XM = M_PI;
    int K = 100;
    double a1 = -0.25, a2 = 0.25;
    int L = 10;
    double H = 0;
    double P = 0.95;
    double E = 0.01;
    double Om, D, Omm, Dm;
    std::vector<zn> zna;
    std::vector<double> aa;
    std::vector<double> AAm;
    double Dict, Sm = 100000;
    double N = n(P, E, XM, Xm);

    std::vector<double> X(101);
    for (auto i = 0; i < K; i++)
        X[i] = Xm + i * (XM - Xm) / K;

    std::vector<double> F(101);
    for (auto i = 0; i <= K; i++)
        F[i] = sin(X[i]) + 0.5;

    std::vector<double> FShu(101);
    for (auto i = 0; i <= K; i++)
    {
        std::mt19937 gen(time(0));
        std::uniform_real_distribution<double> urd(a1, a2);
        double A = urd(gen);
        FShu[i] = F[i] + A;
    }

    for (auto i = 0; i <= K; i++)
        std::cout << "f = " << F[i] << " F s shumom = " << FShu[i] << "\n";


    std::vector<double> FN3;
    std::cout << " H       Distanc             alpha              Omega        delta" << "\n";
    for (auto l = 0; l <= L; l++)
    {
        for (auto i = 0; i < K / 2; i++)
        {
            aa = AA(N3);
            FN3 = fil(FShu, aa, N3);
            Om = om(FN3);
            D = d(FN3, FShu);
            Dict = S(Om, D);

            if (Dict < Sm)
            {
                Sm = Dict;
                Omm = Om;
                Dm = D;
                AAm = aa;
            }
        }
        H = (double)l / L;
        zna.push_back({ H, Sm, AAm, Omm, Dm });
        std::cout << std::fixed << std::setprecision(1) << H << "  " << std::setprecision(4) << Sm << "   [ ";
        for (auto i : AAm)
            std::cout << i << " ";
        std::cout << "]  " << Omm << "   " << Dm << "\n";
        Sm = 100000;
    }

    std::sort(zna.begin(), zna.end(), [](auto x, auto y) {return x.S < y.S; });
    std::cout << "\n" << " h*     J        om        d" << "\n";
    std::cout << std::setprecision(1) << zna[0].h << "   " << std::setprecision(4) << J(zna[0].h, zna[0].om, zna[0].d) << "   " << zna[0].om << "   " << zna[0].d << "\n" << "\n";
    FN3 = fil(FShu, AAm, N3);
    for (auto i : FN3)
        std::cout << i << "\n";
    std::cout << "\n";


    zna.clear();
    AAm.clear();
    Dm = 0;
    Omm = 0;
    std::vector<double> FN5;
    std::cout << " h    dictanc                 alpha                  omega            d" << "\n";
    for (auto l = 0; l <= L; l++)
    {
        for (auto i = 0; i < N; i++)
        {
            aa = AA(N5);
            FN5 = fil(FShu, aa, N5);
            Om = om(FN5);
            D = d(FN5, FShu);
            Dict = S(Om, D);

            if (Dict < Sm)
            {
                Sm = Dict;
                Omm = Om;
                Dm = D;
                AAm = aa;
            }
        }
        H = (double)l / L;
        zna.push_back({ H, Sm, AAm, Omm, Dm });
        std::cout << std::fixed << std::setprecision(1) << H << "  " << std::setprecision(4) << Sm << "   [ ";
        for (auto i : AAm)
            std::cout << i << " ";
        std::cout << "]  " << Omm << "   " << Dm << "\n";
        Sm = 100000;
    }

    std::sort(zna.begin(), zna.end(), [](auto x, auto y) {return x.S < y.S; });
    std::cout << "\n" << " h*     J        w        d" << "\n";
    std::cout << std::setprecision(1) << zna[0].h << "   " << std::setprecision(4) << J(zna[0].h, zna[0].om, zna[0].d) << "   " << zna[0].om << "   " << zna[0].d << "\n" << "\n";
    FN5 = fil(FShu, AAm, N5);
    for (auto i : FN5)
        std::cout << i << "\n";
    std::cout << "\n";

    return 0;
}