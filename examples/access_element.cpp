#include <iostream>

#include "matrix.hpp"

int main()
{
    // playing with a ekumath::Matrix instance
    ekumath::Matrix A(3,3,1.);

    std::cout << "A.rows() = " << A.rows() << std::endl;
    std::cout << "A.cols() = " << A.cols() << std::endl;
    std::cout << "A(2, 2) = " << A(2,2) << std::endl;
    std::cout << "Modifying A(2, 2) to 2.0" << std::endl;
    A(2,2) = 2.;
    std::cout << "A(2, 2) = " << A(2,2) << std::endl;
    std::cout << std::endl;

    // playing with a ekumath::Matrix const instance
    const ekumath::Matrix B(3,3,3.1);

    std::cout << "B.rows() = " << B.rows() << std::endl;
    std::cout << "B.cols() = " << B.cols() << std::endl;
    std::cout << "B(2, 2) = " << B(2,2) << std::endl;
    std::cout << "B can't be modified" << std::endl;
    // B(2,2) = 2.;

    ekumath::Matrix C{{1.0, 2.0, 3.0}, {3.0, 1.0, 2.0}, {2.0, 3.0, 1.0}};
    std::cout << C << std::endl;

    ekumath::Matrix D{1.0, 2.0, 3.0};
    std::cout << D << std::endl;

    ekumath::Matrix E(3, 1);
    E(0,0) = 1.0;
    E(1,0) = 2.0;
    E(2,0) = 3.0;

    std::cout << (D == E) << std::endl;
    return 0;
}