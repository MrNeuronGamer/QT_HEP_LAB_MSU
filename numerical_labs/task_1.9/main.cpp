#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;
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

        for (size_t j = i; j < Size; j++)
            if (x[j] < x[min_ind]) min_ind = j;

        swap(x[i], x[min_ind]);


    }



}

double find_max(double* arr, size_t Size)
{
    double res = arr[0];

    for (size_t i = 1; i < Size; i++)
    {

        if (arr[i] > res) res = arr[i];

    }

    return res;
}

void pop_heap(double* heap,size_t Size, double* dest,double max)
{
    double stopper = max;
    *(dest) = heap[0];

    size_t i = 0;
    size_t target = 0;
    i = (target << 1) + 1;
    while (i < Size)
    {

        if (i + 1 == Size)
        {

            heap[target] = heap[i];
            target = i;
            

        }

        else
        {

            if (heap[i] < heap[i + 1])
            {
                heap[target] = heap[i];
                target = i;
            }

            else
            {
                heap[target] = heap[i+1];
                target = i+1;
            }
        }

        i = (target << 1) + 1;

    }
    heap[target] = stopper;


}

void build_heap(double* arr, size_t Size)
{

    int target;
    int i = Size - 1;

    while (i > 0)
    {

        

        target = ((i - 1) >> 1);
        if (((i-1)>>1) !=  ((i-2)>>1) )
        {
            if(arr[target] > arr[i])
                swap(arr[target] , arr[i]);
            
        }
        else
        {
            if (arr[i] < arr[i -1])
            {
                if (arr[target] > arr[i])
                    swap(arr[target], arr[i]);
            }
            else
            {
                if (arr[target] > arr[i-1])
                    swap(arr[target], arr[i -1]);
            }
                


        }

        i = ((target - 1) << 1) + 2;
    }



}

void heap_to_arr(double* arr, size_t Size)
{

    double* buff = new double[Size];
    double max = find_max(arr, Size);

    for (size_t i = 0; i < Size; i++)
    {
        pop_heap(arr, Size , buff + i, max);

    }


    for (size_t i = 0; i < Size; i++)
    {
        arr[i] = buff[i];

    }



    delete[] buff;

}

void heapsort(double* arr, size_t Size)
{

    build_heap(arr, Size);
    heap_to_arr(arr, Size);

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
    
    double *x, *y;
    x = new double[N];
    y = new double[N];

    fill_x(x);
    fill_x(y);

    auto start = high_resolution_clock::now();
    Neuman_sort(x, N);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Neuman sort has executed for :  " << duration.count() << " mcs" << endl;

    delete[] x;

    start = high_resolution_clock::now();
    heapsort(y, N);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Heap sort has executed for :  " << duration.count() << " mcs" << endl;

    for (int i = 0; i < 7; i++)
        cout << y[i] << "  ";

    delete[] y;
   


   


    int a;
    cin >> a;

    return 0;
}