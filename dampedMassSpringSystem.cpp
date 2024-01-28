// cpp file to define class member functions
#include "dampedMassSpringSystem.h"
#include "matrix.h"

// cpp Libraries
#include <vector>
#include <iterator>
#include <iostream>

// Constructor definition
dampedMassSpringSystem::dampedMassSpringSystem(double simTime, double k, double d, double dt, 
                                    std::vector<double>& m, std::vector<double>& x, std::vector<double>& v)
                :cSimulationTime(simTime), cSpringConstant(k), cDampingCoeff(d), cTimeStep(dt),massVector(m), 
                 displacementVector(x), velocityVector(v){
                     }

// Destructor defintion
dampedMassSpringSystem::~dampedMassSpringSystem(){
// Do nothing
}

// Function definitions
std::vector<std::vector<double>> dampedMassSpringSystem::createMatrix(double k, double d, const std::vector<double>& m){
    unsigned int n = m.size(); // Size of mass vector determines size of system matrix (2nx2n)
    
    // 2D Vector Structure: DATATYPE IDENTIFIER(ROW_COUNT, DATATYPE(COLUMN_COUNT))
    // Initialising vectors needed to make system matrix for uncoupling ODE -> derived by hand
    std::vector<std::vector<double>> systemVector(2*n,std::vector<double>(2*n, 0.0)); // System vector used to solve ODE
    
    // The system matrix (if 3 <= n) should look like the following: 
    // Top left - nxn zeroes matrix
    // Top right - nxn identity matrx
    // Bottom left - nxn banded {-1,2,-1} matrix with each column multiplied by k/m[i]
    // Bottom right - nxn identity matrix, scaled by -d/m[i]
    
    // Note a different matrix would be needed for n == 1 and n == 2, in particular the banded matrix, but should be trivial
    
    // Top left, no work needs to be done as systemVector has already been initialised to be zeros
    
    // Top right, identity matrix
    for (unsigned int i = 0; i < n; i++){
        systemVector[i][i+n] = 1.0;
    }   
    
    // Bottom left, spring terms
    double vScaleFactor; // Variable that computes the scale factor needed for each column (i.e. k/m[i])
    std::vector<double> band = {1.0, -2.0, 1.0}; // Banded vector
    
    if (n == 1){
        systemVector[1][0] = band[1] * k / m[0];
    }
    
    else if (n == 2){
        systemVector[2][0] = band[1] * k / m[0];
        systemVector[3][0] = band[2] * k / m[0];
        systemVector[2][1] = band[0] * k / m[1];
        systemVector[3][1] = band[1] * k / m[1];
    }
    
    else {
        for (unsigned int i = n; i < 2*n; i++){
            systemVector[i][i-n] = band[1]; // Diagonal elements
            if(i < 2 * n - 1) {
                systemVector[i][i-n+1] = band[2]; // Upper off-diagonal elements
            }
            if(i > n) {
                systemVector[i][i-n-1] = band[0]; // Lower off-diagonal elements
            }
        }
        
        // Multiplying each column by scale factor, could be changed later when scalar multiply is implemented???
        for (unsigned int c = 0; c < n; c++){ // For each column
            vScaleFactor = k/m[c];
            for (unsigned int r = 0; r < n; r++){ // For each row
                systemVector[r+n][c] *= vScaleFactor;
            } 
        }
    }

    // Bottom right - damping terms
    if (!(d == 0)){
        for (unsigned int i = n; i < 2*n; i++){
            systemVector[i][i] = -d/m[i-n];
        }     
    }
   
    return systemVector;
}


std::vector<std::vector<double>> dampedMassSpringSystem::createStateVector(const std::vector<double>& x, const std::vector<double>& v){
    std::vector<std::vector<double>> stateVector((x.size()*2), std::vector<double>(1,0.0)); // Initialise state vector
    unsigned int n = x.size();
    for (unsigned int i = 0; i < n; i++){
        stateVector[i][0] = x[i]; // Input initial displacement values
        stateVector[i+n][0] = v[i]; // Input initial velocity values
    }
    
    return stateVector;
}



std::vector<std::vector<double>> dampedMassSpringSystem::explicitEuler(double cSimulationTime, double cTimeStep){
    
    int nSteps = cSimulationTime/cTimeStep; 
    
    std::vector<std::vector<double>> systemMatrix = createMatrix(cSpringConstant, cDampingCoeff, massVector); // Create system matrix
    std::vector<std::vector<double>> stateVector = createStateVector(displacementVector, velocityVector); // Create state vector
    std::vector<std::vector<double>> timeVector(nSteps, std::vector<double>(1, 0.0));
    std::vector<std::vector<double>> result; // Initialising result vector to store results
    
    matrix sm = matrix(systemMatrix);
    matrix sv = matrix(stateVector);
    
    // Push back initial state
    std::vector<double> stateCopy;
    stateCopy.insert(stateCopy.end(),0.0);
    for (const auto& row : sv.getVector()) {
        stateCopy.insert(stateCopy.end(), row.begin(), row.end()); // Flatten the 2D vector to a 1D vector
    }
    result.push_back(stateCopy);
    
    for (int i = 0; i < nSteps; i++){
        sv = (sm * sv) * cTimeStep + sv; // Explicit Euler Solver
        std::vector<double> stateCopy;
        stateCopy.insert(stateCopy.end(), (i+1)*cTimeStep);
        for (const auto& row : sv.getVector()) {
            stateCopy.insert(stateCopy.end(), row.begin(), row.end()); // Flatten the 2D vector to a 1D vector
        }
        result.push_back(stateCopy); // Store the flattened state vector
    }
    
    return result;
}

std::vector<std::vector<double>> dampedMassSpringSystem::rungeKuttaFour(double cSimulationTime, double cTimeStep){
    
    int nSteps = cSimulationTime/cTimeStep; 
    
    std::vector<std::vector<double>> systemMatrix = createMatrix(cSpringConstant, cDampingCoeff, massVector); // Create system matrix
    std::vector<std::vector<double>> stateVector = createStateVector(displacementVector, velocityVector); // Create state vector
    std::vector<std::vector<double>> timeVector(nSteps, std::vector<double>(1, 0.0));
    std::vector<std::vector<double>> result; // Initialising result vector to store results
    
    matrix sm = matrix(systemMatrix);
    matrix sv = matrix(stateVector);
    
    // Push back initial state
    std::vector<double> stateCopy;
    stateCopy.insert(stateCopy.end(),0.0);
    for (const auto& row : sv.getVector()) {
        stateCopy.insert(stateCopy.end(), row.begin(), row.end()); // Flatten the 2D vector to a 1D vector
    }
    result.push_back(stateCopy);
    
    for (int i = 0; i < nSteps; i++){
        
        // Initialise temporary 2d vectors to store k coefficients
        std::vector<std::vector<double>> a;
        std::vector<std::vector<double>> b;
        std::vector<std::vector<double>> c;
        std::vector<std::vector<double>> d;
        std::vector<std::vector<double>> e;
        
        matrix k1 = matrix(a);
        matrix k2 = matrix(b);
        matrix k3 = matrix(c);
        matrix k4 = matrix(d);
        matrix ph1 = matrix(e); // Placeholder matrix
        
        // RK4 Begins Here
        
        k1 = sm * sv; // Find k1
        
        k2 = k1 * (cTimeStep/2.0) + sv; // Intermediate step, cannot do parantheses in operator overload T_T
        k2 = sm * k2; // Find k2
        
        k3 = k2 * (cTimeStep/2.0) + sv; // Intermediate step
        k3 = sm * k3; // Find k3
        
        k4 = k3 * cTimeStep + sv;  // Intermediate step
        k4 = sm * k4; // Find k4
        
        ph1 = k2 * 2.0;
        ph1 = k3 * 2.0 + ph1; // Placeholder for final calculation
        
        sv = (k1 + k4 + ph1) * (cTimeStep/6.0) + sv; // Update stateVector
        
        // RK4 Ends Here
        
        std::vector<double> stateCopy;
        stateCopy.insert(stateCopy.end(), (i+1)*cTimeStep);
        for (const auto& row : sv.getVector()) {
            stateCopy.insert(stateCopy.end(), row.begin(), row.end()); // Flatten the 2D vector to a 1D vector
        }
        result.push_back(stateCopy); // Store the flattened state vector
    }
    
    return result;
}

