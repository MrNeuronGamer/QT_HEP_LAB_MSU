#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;
static int Epsi_global = 10000;

class LE_System
{

public:
    double*     lhs_matrix;
    double*     rhs;
    size_t      Size;
    double*     solution;

            LE_System       (void(*fill_f), size_t size);               // creates a system of linear equations of size "size" using "fill_f" generating function; lhs an rhs get created; field Size == size
            LE_System       ();
            ~LE_System      ();                                         // frees all the fields' memory and destroys object
    void    operator=       (const LE_System&);
    void    operator++      (int);
    void    FindSolution    ();
    void    operator*       (double*);


    friend ostream& operator<<(ostream&, LE_System&);


    // debug section:
    void printTables() const;           // prints lhs and rhs 
    LE_System(double* Lines, size_t size);

private:
    void    SwapColumn      (size_t start_pos, size_t target_pos);      // swaps columns indexed by start_pos and target_pos
    void    SwapLine        (size_t start_pos, size_t target_pos);      // swaps lines indexed by start_pos and target_pos
    double  GetNorm         (size_t line) const;                        // determines a L1 norm for a line and returns it
    void    MultiplyLine    (size_t line, double coeff);                // multiplyes a line with a coeff
    bool    CheckResult     (double error = 1) const;                   // checks whether a solution found recreates a column b with a residual less then "error" lines-wise / default error = 1 / returns true if solutions is ok
    void    SubtractLines   (size_t from, size_t what);                 // subtracts line "what" from line "from"
    void    SettleLine      (size_t line);                              // makes shure that any line under "line" gets subtracted with "line" and gets rid of a variable indexed as "line"
    void    ReduceBigNorms  ();                                         // makes sure that no line's norm is 10^2 times bigger comparing with average norm 
    
};


static LE_System Test;

int main()
{
    LE_System A;
    
    ofstream OUT;
    OUT.open("results_LE.dat", ios_base::ate);

    for (size_t i = 0 ;i < 10; i++)
    {

        A.FindSolution();
        
        OUT << A;
        A++;
        cout << " iteration :: " << i << "\n";

    }

   

    OUT.close();



    return 0;
}


ostream& operator<<(ostream& out, LE_System& A)
{

    double E_1 = 0;
    double E_2 = 0;
    double x = 0;


    for (size_t i = 0; i < A.Size; i++)
    {
        E_1 += (A.solution[i] - x)*(A.solution[i] - x);
        
        x += 0.001;
    }   
    E_1 = sqrt(E_1);


    Test*A.solution;

    for (size_t i = 0; i < A.Size; i++)
    {
        E_2 += (A.solution[i] - A.rhs[i])*(A.solution[i] - A.rhs[i]);
       
    }

    E_2 = sqrt(E_2);

    out << E_1 << "     :::     " << E_2 << "\n\n";
    
    return out;
}

void LE_System::SwapColumn(size_t start_pos, size_t target_pos)
{   
    
    for (size_t ci = 0; ci < Size; ci++)
        swap(lhs_matrix[ci*Size + start_pos], lhs_matrix[ci*Size + target_pos]);

}

void LE_System::SwapLine(size_t start_pos, size_t target_pos)
{
    for (size_t ci = 0; ci < Size; ci++)
        swap(lhs_matrix[start_pos*Size + ci], lhs_matrix[target_pos*Size + ci]);

    swap(rhs[start_pos], rhs[target_pos]);

}

double LE_System::GetNorm(size_t line) const 
{
    double norm = 0;
    for (size_t i = 0; i< Size; i++)
    {
        norm += abs(lhs_matrix[line*Size + i]);
    }

    return norm;

}

void LE_System::MultiplyLine(size_t line, double coeff)
{

    for (size_t i = 0; i< Size; i++)
    {
        lhs_matrix[line*Size + i]*=coeff;
    }

    rhs[line] *= coeff;


}

void LE_System::SubtractLines(size_t from, size_t what)
{

    for (size_t le = 0; le < Size; le++)
    {
        lhs_matrix[le + from*Size] -= lhs_matrix[le + what*Size];
    }

    rhs[from] -= rhs[what];

}

void LE_System::SettleLine(size_t line)
{
    MultiplyLine(line, 1 / lhs_matrix[line + Size*line]);
    for (size_t l = line + 1; l < Size; l++)
    {
        MultiplyLine(l, 1/lhs_matrix[line+Size*l]);        
        SubtractLines(l, line);
    }

}

void LE_System::operator= (const LE_System& A)
{

    Size = A.Size;  
    for (size_t i = 0; i < Size; i++)
    {
        for (size_t j = 0; j < Size; j++)
            lhs_matrix[i*Size + j] = A.lhs_matrix[i*Size + j];

        rhs[i] = A.rhs[i];
        solution[i] = A.solution[i];
    }

}

void LE_System::ReduceBigNorms()
{
    bool flag = true;
    while (flag)
    {
        flag = false;
        double AvNorm = 0;
        vector<double> Norms;
        Norms.resize(Size);

        for (size_t i = 0; i < Size; i++)
        {
            Norms[i] = GetNorm(i);
            AvNorm += Norms[i];
        }

        AvNorm = AvNorm / Size;

        for (size_t i = 0; i < Size; i++)
        {
            if (Norms[i] > AvNorm * 100)
            {
                MultiplyLine(i, Norms[i] / AvNorm);
                flag = true;
            }
        }


    }

}



void LE_System::FindSolution()
{

    ReduceBigNorms();           // get rid of big norms    
    for (size_t i = 0; i < Size; i++)
    {
        SettleLine(i);           // triangulate
    }
   

    for (size_t i = Size - 1; Size - i != Size+1; i--)
    {
        solution[i] = rhs[i];
        for (size_t j = i + 1; j < Size; j++)
        {
            solution[i] -= solution[j] * lhs_matrix[i*Size+j];
        }
    }
   
    cout <<"Solution is correct :  " <<  CheckResult(1e-8) << "    ";
   
}

bool LE_System::CheckResult(double error ) const
{

    for (size_t i = 0; i < Size; i++)
    {
        double LHS = 0;
        for (size_t j = 0; j < Size; j++)
        {
            LHS += solution[j] * lhs_matrix[i*Size + j];
        }

        if (abs(LHS - rhs[i]) > error) return false;
    }

    return true;
}

void LE_System::operator++(int)
{

    for (size_t i = 0; i < Size; i++)
        lhs_matrix[i*Size + i] += Epsi_global;
        
    Epsi_global = Epsi_global / 2;

}

LE_System::~LE_System()
{
    delete[] lhs_matrix;
    delete[] rhs;
    delete[] solution;

}

LE_System::LE_System()
{

    Size = 1001;
    lhs_matrix = new double[Size*Size];
    rhs = new double[Size];
    solution = new double[Size];


    double dx = 0.001;
    double a = 0.765;
    double x[1001];

    for (size_t i = 0; i < Size; i++)
        x[i] = i*dx;



    auto c = [](size_t j)
    {   
        if (j == 0) return 1. / 3;
        if (j == 1000) return 1. / 3;
        if (j % 2 == 0) return 4. / 3;
        return 2. / 3;
    
    };

    auto g = [&a](double x)
    {
        return 0.5*log(1 + (1 - 2 * x) / (a*a + x*x)) + (x* (atan((1 - x) / a) + atan(x / a))) / (a);

    };


    for (size_t i = 0; i < Size; i++)
    {
        for (size_t j = 0; j < Size; j++)
            lhs_matrix[i*Size + j] = dx*c(j) / ((x[i] - x[j])*(x[i] - x[j]) + a*a);

        rhs[i] = g(x[i]);
    }




}




LE_System::LE_System(double* arr, size_t size)
{

    lhs_matrix = new double[size*size];
    solution = new double[size];
    rhs = new double[size];
    Size = size;

    for (size_t i = 0; i < Size; i++)
    {

        for (size_t j = 0; j < Size; j++)
        {
            lhs_matrix[i*Size + j] = arr[i*(Size + 2) + j];

        }
        solution[i] = arr[i*(Size + 2) + Size];
        rhs[i] = arr[i*(Size + 2) + Size+1];

    }


}

void LE_System::operator*(double* b)
{

    double* buffer = new double[Size];


    for (size_t i = 0; i < Size; i++)
    {
        buffer[i] = b[i];
        b[i] = 0;
    }

    for (size_t i = 0; i < Size; i++)
        for (size_t j = 0; j < Size; j++)
        {
            b[i] += buffer[j] * lhs_matrix[i*Size + j];            
        }



    delete[] buffer;


}

void LE_System::printTables() const
{

    for (size_t i = 0; i < Size; i++)
    {
        for (size_t j = 0; j < Size; j++)
            cout << lhs_matrix[i*Size + j] << " ";

        cout << "  ::  " << solution[i] << "  ::  " << rhs[i] << "\n";
    }

    int c = 0;
    cin >> c;
    cout << "\n\n\n";
}