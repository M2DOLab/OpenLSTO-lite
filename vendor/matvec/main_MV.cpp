#include "./MatrixM2DO.cpp" // Lazy compiliing
// #include "./VectorM2DO.h"

int main(){

    // testing static matrix
    Matrix<double, 3, 2> matrixd_3_2; 
    Matrix<double, 3, 3> matrixd_3;
    Matrix<double, -1, -1> matrixd_3i(3,3);

    matrixd_3i.data = {{1,2,3},{4,5,6},{7,8,9}};


    matrixd_3_2.fill(2.4);
    matrixd_3_2(0,1) = 7.0;

    matrixd_3_2.print();

    matrixd_3.fill(8.0);
    matrixd_3(1,2) = 0.0;
    matrixd_3(0,1) = 0.0;

    Matrix<double, 3, 2> matrixdot = matrixd_3.dot(matrixd_3_2);
    matrixdot.print();

    Matrix<double, 3, 3> matrixd_inv = matrixd_3.inverse();
    Matrix<double, 2, 3> matrixd_T = matrixd_3_2.transpose();
    double det = matrixd_3.determinant();

    matrixd_inv.print();
    matrixd_T.print();
    std::cout << det << std::endl;




    std::cout<<"\n\nfor debug1 ----\n";
    //Matrix<int,2,2> A1,A2,A3;
    Matrix<int,-1,-1> A1(2,2),A2(2,2),A3(2,2);
    A1.data = {{1,2},{3,4}};
    A2.data = {{1,2},{3,4}};
    // A3 = A1.dot(A2);
    std::cout << A3.size() << endl;
    A3.print();
    std::cout<<"\nfor debug1 ----\n\n";




    std::cout << "\n\nfor debug2 ----\n";
    Matrix<int,-1,-1> mat_A(3,4);
    mat_A.data = {{1,-2,0,4},{-7,3,1,-3},{4,5,8,0}};
    std::cout << "\n[A] =\n";
    mat_A.print();

    Matrix<int,-1,-1> mat_B(4,4);
    mat_B.data = {{4,2,0,-5},{2,6,-2,2},{4,3,4,3},{-1,8,-1,0}};
    std::cout << "\n[B] =\n";
    mat_B.print();

    Matrix<int,-1,-1> mat_C(3,4);
    mat_C.data = {{-3,2,3,0},{1,0,-1,1},{4,0,-8,5}};
    std::cout << "\n[C] =\n";
    mat_C.print();
    
    std::cout << "\n ******* TEST *******\n";
    mat_A = mat_A.dot(mat_B);
    std::cout << "\n[A] = [A]*[B] =\n";
    mat_A.print(); // js985: Works fine

    Matrix<int,-1,-1> mat_D; // size unknown a-priori
    mat_D = mat_A.dot(mat_B.dot(mat_C.transpose()));
    std::cout << "\n[D] = [A]*[B]*[C]^T =\n";
    mat_D.print();
    std::cout << "\nfor debug2 ----\n\n\n";



    Vector<double,3> v3;
    Vector<double,-1> v3i(3);

    v3 *= 2;
    v3.data = {1,2,3};
    v3.print();
    v3i.data = {5,7,9};
    v3i += 1;
    v3i.print();
    std::cout << "\nDEBUG" << v3(0) << std::endl;

    Vector<double,3> v3_dot = matrixd_3.dot(v3);
    v3_dot.print();

    Vector<double,-1> v3_doti = matrixd_3i.dot(v3i);
    v3_doti.print();

    return 0;
}