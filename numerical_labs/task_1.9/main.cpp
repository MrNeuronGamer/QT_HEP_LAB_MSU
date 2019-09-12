#include <iostream>
#include <cmath>

using namespace std;
typedef unsigned int ui;


const ui R = UINT_MAX;
ui A = 999;
ui B = 9999;
const size_t N = 33000000;

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

void direct_ins_sort(double* x, size_t Size)
{
    

    size_t min_ind = 0;
    for (size_t i = 0; i < Size; i++)
    {
        min_ind = i;

        for (size_t j = i ; j < Size; j++)
            if (x[j] < x[min_ind]) min_ind = j;

        swap(x[i], x[min_ind]);


    }

    

}

void heapsort(double* arr, size_t Size)
{



}

void merge(double* a, size_t as, double* b, size_t bs)
{
    size_t rs = as + bs;
    double* res = new double[rs];
    size_t i = 0;
    size_t j = 0;

    while ((i < as) && (j < bs))
    {

        if (a[i] < b[j])
        {
            res[i + j] = a[i];
            i++;

        }
        else
        {
            res[i + j] = b[j];
            j++;
        }

      

   }


    if (i == as)
    {
        for (; j < bs; j++)
            res[i + j] = b[j];
    }
    else
    {
        for (; i < as; i++)
            res[i + j] = a[i];
    }

    for (i = 0; i < rs; i++)
    {
        a[i] = res[i];
    }


    delete[] res;
    


}

void Neuman_sort(double* x, size_t Size)
{


    if (Size > 8)
    {

        Neuman_sort(x, Size / 2);
        Neuman_sort(x + Size / 2, Size - Size / 2);

        merge(x, Size / 2, x + Size / 2, Size - Size / 2);

    }
    else direct_ins_sort(x, Size);



}

int main()
{

    

    return 0;
}