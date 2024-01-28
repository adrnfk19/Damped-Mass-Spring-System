// Header file to define class Matrix for matrix multiplication
#ifndef CLASS_MATRIX
#define CLASS_MATRIX

#include <vector>

class matrix{
private:
    std::vector<std::vector<double>> Mat; // The matrix itself declared as a 2D vector

public:
    matrix(std::vector<std::vector<double>> vec);  // Constructor
    ~matrix(); // Destructor
    
    // Operator Overloading for Matrix Operations
    matrix operator+(matrix& nextMat); // Matrix addition
    matrix operator*(matrix& nextMat); // Matrix multiplication
    matrix operator*(double scalar); // Scalar multiplication
    
    // Function declaration
    std::vector<std::vector<double>> getVector(){
        return Mat; // returns vector instead of matrix datatype
    }
    
};

#endif