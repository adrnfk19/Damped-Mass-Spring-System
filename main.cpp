// Import classes
#include "dampedMassSpringSystem.h"
#include "matrix.h"

// Import libraries
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

// Function declaration
static void readFile(int&, double&, double&, double&, double&, std::vector<double>&, std::vector<double>&, std::vector<double>&); // Function that reads from parameters.txt and overwrites k,d,dt and T and pushes back into mass vectors
static void writeFile(std::vector<std::vector<double>>&); // Function that writes results into output file


int main(int argc, char **argv)
{
    // Initialise constants
    int solverChoice = 0;
    double T, dt, k ,d;
    
    // Initialise vectors
    std::vector<double> vMass, vDisplacement, vVelocity;
    std::vector<std::vector<double>> systemMatrix, stateVector, result;
    
    // Main script
    try{
        readFile(solverChoice, T, dt, k ,d, vMass, vDisplacement, vVelocity); // Import data from file and store in glocal scoped variables
    }
    catch (const std::domain_error& e) {
        std::cout << "Error occured: " << e.what() << std::endl;
    }
    catch (const std::logic_error& e) {
        std::cout << "Error occured: " << e.what() << std::endl;
    }
    
    dampedMassSpringSystem sys(T, k, d, dt, vMass, vDisplacement, vVelocity); // Construct the mass spring system

    // Switch case for solverChoice
   switch (solverChoice){
        case 0:
            result = sys.explicitEuler(T, dt); // Explicit Euler Solver if solverChoice == 0
            break;
        case 1:
            result = sys.rungeKuttaFour(T, dt); // RK4 Solver if solverChoice == 1
            break;
        default:
            return 0;
    }

    // Destructor called automatically
    
    try{
        writeFile(result); // Write to output file
    }
    catch (const std::logic_error& e){
        std::cout << "Error occured: " << e.what() << std::endl;
        }
    
	return 0;
}



// Function definition
void readFile(int& solverChoice, double& T, double& dt, double& k, double& d,
            std::vector<double>& massVector, std::vector<double>& displacementVector, std::vector<double>& velocityVector){
    
    // Initialise variables
    int lineCount = 1;
    double mass, displacement, velocity; 
                
    std::ifstream vMyFile("parameters.txt"); // Find input file
    std::string firstLine; // Variable that reads first line of text file

    if (vMyFile.good()){ // If file is in good state, run the following block of code
    
        if (std::getline(vMyFile,firstLine)){
            std::istringstream firstString(firstLine); // Use of stringstream to store first line of txt file
            firstString >> solverChoice >> T >> dt >> k >> d; // Store values from first line to respective variables
        }
        
        while (vMyFile >> mass >> displacement >> velocity){
 
            lineCount += 1;

            if (mass < 0) {
                throw std::domain_error("Mass cannot be negative!"); // Throws exception if mass is negative
            }

            // Push back to vectors to store respective data
            massVector.push_back(mass);
            displacementVector.push_back(displacement);
            velocityVector.push_back(velocity);
                
            if (vMyFile.eof()){
                break;
            }
        }
  
        if (lineCount == 1) {
            throw std::logic_error("No data was input!"); // Check if file has no data
        }

        vMyFile.close(); // Close file
    }
    
    else{
        throw std::logic_error ("Unable to read file"); // Throw error to let user know that file cannot be read
    }
}

void writeFile(std::vector<std::vector<double>>& result){
    
    std::ofstream outputFile("output.txt"); // Find output file
    
    if (outputFile.is_open()){ // If file is open, run the following block of code
    
        // Construct header
        std::string header = "t";
        for (unsigned int i = 1; i <= (result[0].size()-1)/2; i++) {
            header += " x_" + std::to_string(i) + "(t)" + " v_" + std::to_string(i) + "(t)";
        }
        outputFile << header << "\n";
    
        for (const auto& row : result) {

            outputFile << row[0]; // Write time, which was added in solver function

            // Write dynamic motion of each mass
            for (unsigned int i = 1; i <= ((result[0].size() - 1) / 2); i++) {
                outputFile << " " << std::setprecision(5) << row[i]; // Write displacement values
                outputFile << " " << std::setprecision(5) << row[i+(result[0].size()-1)/2]; // Write velocity values
            }

            outputFile << "\n"; // Move to next line for next state
        }

        outputFile.close(); // Close file
    }
    
    else {
        throw std::logic_error ("Unable to write on file!"); // Throw error to let user know that file cannot be written on
    }
}

