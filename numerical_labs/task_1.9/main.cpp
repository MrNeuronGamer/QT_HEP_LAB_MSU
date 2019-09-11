#include <iostream>
#include <cmath>

using namespace std;
typedef unsigned int ui;


const ui R = UINT_MAX;
ui A = 999;
ui B = 9999;
size_t N = 33000000;

void generate_ab(ui& a, ui& b)
{

    a = a * 1664525 + 1013904223;
    b = b * 1664525 + 1013904223;
}

void generate_x(double& x_1, double& x_2)
{

    double s = 2;
    double u, v;
    
    while (s > 1)
    {
        generate_ab(A, B);

        u = 2 * double(A) / R - 1;
        v = 2 * double(B) / R - 1;
        s = u*u + v*v;

    }

    x_1 = u*sqrt(-2.0 / s*log(s));
    x_2 = v*sqrt(-2.0 / s*log(s));

    
}

void fill_x(double* x)
{

    for (size_t i = 0; i < N; i += 2)
    {

        generate_x(x[i], x[i + 1]);

    }

}

void heapsort(double* arr)
{


}


void Neuman_sort(double* x)
{




}

int main()
{
      

    double* x = NULL;

    x = new double [N];

    fill_x(x);



    


   
    delete[] x;


    return 0;
}