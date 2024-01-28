// Header file to define class dampedMassSpringSystem
#ifndef CLASS_DAMPEDMASSSPRINGSYSTEM
#define CLASS_DAMPEDMASSSPRINGSYSTEM

#include <vector>

class dampedMassSpringSystem{
    
private:
    double cSimulationTime;
    double cSpringConstant;
    double cDampingCoeff;
    double cTimeStep;
    
    std::vector<double> massVector;
    std::vector<double> displacementVector;
    std::vector<double> velocityVector;
        
public:
    dampedMassSpringSystem(double simTime, double k, double d, double dt, std::vector<double>& m, std::vector<double>& x, std::vector<double>& v); // Constructor
    ~dampedMassSpringSystem(); // Destructor
            
    // Declare member functions, definitions of member functions defined in "dampedMassSpringSystem.cpp"
    std::vector<std::vector<double>> createMatrix(double cSpringConstant, double cDampingCoeff, const std::vector<double>& massVector); // Create system vector for system of ODES 
    std::vector<std::vector<double>> createStateVector(const std::vector<double>& displacementVector, const std::vector<double>& velocityVector); // Create state vector here 

    // Declare functions for numerical solvers
    std::vector<std::vector<double>> explicitEuler(double cSimulationTime, double cTimeStep); // Explicit Euler
    std::vector<std::vector<double>> rungeKuttaFour(double cSimulationTime, double cTimeStep); // Runge-Kutta 4
    
};

#endif 