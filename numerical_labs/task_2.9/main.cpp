#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class LE_System
{

public:
    double*     lhs_matrix;
    double*     rhs;
    size_t      Size;
    double*     solution;

            LE_System       (void(*fill_f), size_t size);               // creates a system of linear equations of size "size" using "fill_f" generating function; lhs an rhs get created; field Size == size
            ~LE_System      ();                                         // frees all the fields' memory and destroys object
    void    operator=       (const LE_System&);
    void    operator++      (int);
    void    FindSolution    ();


    friend ostream& operator<<(ostream&, const LE_System&);


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



int main()
{
    double arr[] = { 12,0,0,3,8,7,0,1 };

    LE_System A(arr, 2);
    A.printTables();
    A.FindSolution();



    return 0;
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

    printTables();              //check matrix view

    for (size_t i = Size - 1; Size - i != Size+1; i--)
    {
        solution[i] = rhs[i];
        for (size_t j = i + 1; j < Size; j++)
        {
            solution[i] -= solution[j] * lhs_matrix[i*Size+j];
        }
    }

   

    if (!CheckResult())
        cout << "\nSmth gone wrong\n\n\n\n";
    else cout << "\n\n All's fine Dude! ;) \n";

    printTables();



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


LE_System::~LE_System()
{
    delete[] lhs_matrix;
    delete[] rhs;
    delete[] solution;

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