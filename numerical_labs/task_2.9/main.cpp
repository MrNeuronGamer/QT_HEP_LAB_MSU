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
    void    MoveColumn      (size_t start_pos, size_t target_pos);      // swaps columns indexed by start_pos and target_pos
    void    MoveLine        (size_t start_pos, size_t target_pos);      // swaps lines indexed by start_pos and target_pos
    double  GetNorm         (size_t line_begin);                        // determines a L1 norm for a line and returns it
    void    MultiplyLine    (double coeff);                             // multiplyes a line with a coeff
    bool    CheckResult     (double error);                             // checks whether a solution found recreates a column b with a residual less then "error" lines-wise


};



int main()
{



    return 0;
}




void LE_System::MoveColumn(size_t start_pos, size_t target_pos)
{


}