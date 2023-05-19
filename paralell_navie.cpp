 #include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <random>
#include <omp.h>
#include <chrono>
using namespace std;

std::vector<std::vector<double>> readData(){
	std::vector<std::vector<double>> data; 
	std::ifstream file("H:/Courses/Parallel Commputing/Project/data/AAPL1.txt");
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

void makeData(vector<vector<double>>& A, vector<vector<double>>& y,int num){
	int r = A.size();
	int c = A[0].size();
	omp_set_num_threads(num);
	for (int i = 0; i < r;i++){
		y[i][0] = A[i][c - 1];
		A[i][c - 1] = 1;
	}
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dist(0 ,9);
	
	# pragma omp parallel for collapse (2)
	for (int i = 0; i < r;i++){
		for(int j = 0; j < c; j++){
			if(A[i][j] == A[i][j + 1]){
				A[i][j] = A[i][j] + dist(gen) / 1000000000000;
			}
		}
	}
		
}

vector<vector<double>> generate_identity(int size)
{
    vector<vector<double>> I;

    for (int i = 0; i < size; i++)
    {
        vector<double> row;
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                row.push_back(1);
                continue;
            }
            row.push_back(0);
        }
        I.push_back(row);
    }
    return I;
}


vector<vector<double>> Mul(const vector<vector<double>> A, const vector<vector<double>> B,int num) {
    int m = A.size(), n = A[0].size(), p = B[0].size();
    if (n != B.size()) {
       	printf("Khong nhan duoc\n");
        return {};
    }
		omp_set_num_threads(num);
    vector<vector<double>> C(m, vector<double>(p));
   # pragma omp parallel for collapse (2)
	for (int i = 0; i < m; i++) {
    for (int j = 0; j < p; j++) {
      double sum = 0;
      for (int k = 0; k < n; k++) {
        sum += A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
	return C;
}

vector<vector<double>> Inv(vector<vector<double>> &input_matrix,int num){
 signed int size = input_matrix.size();
    signed int i = 0;
    	omp_set_num_threads(num);
    vector<vector<double>> I = generate_identity(size);
    #pragma omp parallel for
	for (i = 0; i < size; i++)
    {
        if (input_matrix[i][i] == 0)
        {
            for (int j = i + 1; j < size; j++)
            {
                if (input_matrix[j][i] != 0.0)
                {
                    swap(input_matrix[i], input_matrix[j]);
                    break;
                }
                if (j == size - 1)
                {
                    cout << "Inverse does not exist for this matrix";
                    exit(0);
                }
            }
        }
        double scale = input_matrix[i][i];
        for (int col = 0; col < size; col++)
        {
            input_matrix[i][col] = input_matrix[i][col] / scale;
            I[i][col] = I[i][col] / scale;
        }
        if (i < size - 1)
        {
            for (int row = i + 1; row < size; row++)
            {
                double factor = input_matrix[row][i];
                for (int col = 0; col < size; col++)
                {
                    input_matrix[row][col] -= factor * input_matrix[i][col];
                    I[row][col] -= factor * I[i][col];
                }
            }
        }
    }
    # pragma omp parallel for collapse (2)
    for (int zeroing_col = size - 1; zeroing_col >= 1; zeroing_col--)
    {
        for (int row = zeroing_col - 1; row >= 0; row--)
        {
            double factor = input_matrix[row][zeroing_col];
            for (int col = 0; col < size; col++)
            {
                input_matrix[row][col] -= factor * input_matrix[zeroing_col][col];
                I[row][col] -= factor * I[zeroing_col][col];
            }
        }
    }
    return I;
}

vector<vector<double>> Trans(const vector<vector<double>>& A,int num){
	int m = A.size();
  	int n = A[0].size();
	omp_set_num_threads(num);
  	vector<vector<double>> C(n, vector<double>(m));
	# pragma omp parallel for collapse (2)
  	for (int i = 0; i < m; i++) {
    	for (int j = 0; j < n; j++) {
      	C[j][i] = A[i][j];
    	}
  	}

  	return C;
}

void print_Matrix(vector<vector<double>> data,int r,int c){
	for (int i = 0; i < r; i++) {
    	for (int j = 0; j < c; j++) {
      		cout<< data[i][j] << " ";
    	}
    	cout<< endl;
	}
}

void Run_in_Threat(vector<vector<double>> data, vector<vector<double>> y, int num){

	vector<vector<double>> data_trans = Trans(data,num);
        
	vector<vector<double>> temp1 = Mul(data_trans,data,num);
    
	vector<vector<double>> temp2 = Inv(temp1,num);

	vector<vector<double>> temp3 = Mul(data_trans, y,num);

	vector<vector<double>> temp4 = Mul(temp1, temp3,num);
}


int main() {
	vector<vector<double>> data = readData();

    vector<vector<double>> y(data.size(), vector<double>(1));
    makeData(data,y,10);

	for(int i = 10; i < 3000;i+=50){
		auto start_serial = chrono::high_resolution_clock::now();
		Run_in_Threat(data,y,i);
		auto end_serial = chrono::high_resolution_clock::now();
		auto duration_serial = chrono::duration_cast<chrono::milliseconds>(end_serial - start_serial);
		cout<< duration_serial.count() <<endl;
	}

	return 0;
}






