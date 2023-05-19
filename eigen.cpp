#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <string>
#include <Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;

vector<vector<double>> readData(){
    // code to read data from file
    std::vector<std::vector<double>> data; 
    std::ifstream file("./data/AAPL.txt");
    if (!file.is_open()) {
        std::cout << "Khong the mo file." << std::endl;
        return data;
    }
    
    std::string line;
    while (std::getline(file, line)) { 
        std::vector<double> row; 
        std::stringstream ss(line); 
        std::string cell;
        while (std::getline(ss, cell, ',')) { 
            double value = std::stold(cell); 
            row.push_back(value); 
        }
        data.push_back(row); 
    }
        
    file.close();
    return data;
}

void makeData(vector<vector<double>>& A, vector<vector<double>>& y){
    // code to modify data
    int r = A.size();
    int c = A[0].size();
    
    for (int i = 0; i < r;i++){
        y[i][0] = A[i][c - 1];
        A[i][c - 1] = 1;
    }
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0 ,9);
    
    for (int i = 0; i < r;i++){
        for(int j = 0; j < c - 1; j++){
            if(A[i][j] == A[i][j + 1]){
                A[i][j] = A[i][j] + dist(gen) / 1000000000000;
            }
        }
    }
}

int main() {
    // read data
    vector<vector<double>> data = readData();
	int rows = data.size();
    int cols = data[0].size();
	
	vector<vector<double>> y_data(rows, vector<double>(1));
    makeData(data, y_data);
	
    // create Eigen matrix
   
    MatrixXd A(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A(i, j) = data[i][j];
        }
    }
    
    // create y matrix
    MatrixXd y(rows, 1);
    
    for (int i = 0; i < rows; i++) {
        y(i, 0) = y_data[i][0];
    }
    
    // multiply matrices
    MatrixXd result = (A.transpose() * A).inverse() * A.transpose() * y;

    return 0;
}
