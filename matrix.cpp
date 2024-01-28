// cpp file to define class member functions
#include "matrix.h"

// cpp Libraries
#include <vector>
#include <iterator>
#include <iostream>

// Constructor definition
matrix::matrix(std::vector<std::vector<double>> vec):Mat(vec) {
    }

// Destructor definition
matrix::~matrix(){
    // Do nothing
}


// Operator Overloading for matrix operations
// Matrix Addition
matrix matrix::operator+(matrix& nextMat){
    std::vector<std::vector<double>> tempMat; // Temporary matrix to store results
    
    // Ensure both matrices are of the same size
    if (Mat.size() != nextMat.Mat.size() || Mat[0].size() != nextMat.Mat[0].size()){
        throw std::logic_error("Vectors are not the same size!");
        return matrix(tempMat); // Return empty tempMat to break out of operation
    }
    
    for (unsigned int i = 0; i < Mat.size(); i++) { // Note to self: .size() gives number of rows; [0].size() gives col
        std::vector<double> tempRow; // Temporary 1D vector to hold results
        for (unsigned int j = 0; j < Mat[0].size(); j++){
            tempRow.push_back(Mat[i][j] + nextMat.Mat[i][j]); // Push back sum to temporary row vector
        }
        tempMat.push_back(tempRow); // Push back temporary row vector to temporary matrix
    }
    
    return matrix(tempMat);
}

// Matrix multiplication
matrix matrix::operator*(matrix& nextMat){

    // Initialize the resulting matrix with zeros
    std::vector<std::vector<double>> tempMat(Mat.size(), std::vector<double>(nextMat.Mat[0].size(), 0.0));
    
    // Check if the matrices can be multiplied correctly
    if (Mat[0].size() != nextMat.Mat.size()) {
        throw std::logic_error("Matrices cannot be multiplied: Incorrect dimensions");
    }

    // Perform matrix multiplication (very cursed O(N^3))
    for (unsigned int i = 0; i < Mat.size(); i++) {
        for (unsigned int j = 0; j < nextMat.Mat[0].size(); j++) {
            for (unsigned int k = 0; k < Mat[0].size();  k++) {
                tempMat[i][j] += Mat[i][k] * nextMat.Mat[k][j]; // Store product into temporary matrix
            }
        }
    }
    
    return matrix(tempMat);
}

// Scalar Multiplication
matrix matrix::operator*(double scalar){
    std::vector<std::vector<double>> tempMat = Mat; // Temporary matrix to store matrix to be scalar multiplied
    
    // Perform scalar multiplication
    for (auto& row : tempMat) {
        for (double& element : row) {
            element *= scalar;
        }
    }
    
    return matrix(tempMat);
}





