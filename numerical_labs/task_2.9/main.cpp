#include <iostream>
#include <fstream>

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

private:
    void    SwapColumn      (size_t start_pos, size_t target_pos);      // swaps columns indexed by start_pos and target_pos
    void    SwapLine        (size_t start_pos, size_t target_pos);      // swaps lines indexed by start_pos and target_pos
    double  GetNorm         (size_t line);                              // determines a L1 norm for a line and returns it
    void    MultiplyLine    (size_t line, double coeff);                // multiplyes a line with a coeff
    bool    CheckResult     (double error = 1);                         // checks whether a solution found recreates a column b with a residual less then "error" lines-wise / default error = 1
    void    SubtractLines   (size_t from, size_t what);                 // subtracts line "what" from line "from"
    void    SettleLine      (size_t line);                              // makes shure that any line under "line" gets subtracted with "line" and gets rid of a variable indexed as "line"

};



int main()
{



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

double LE_System::GetNorm(size_t line)
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




LE_System::~LE_System()
{
    delete[] lhs_matrix;
    delete[] rhs;
    delete[] solution;

}
