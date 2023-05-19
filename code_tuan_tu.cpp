 #include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <random>
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

void makeData(vector<vector<double>>& A, vector<vector<double>>& y){
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


vector<vector<double>> Mul(const vector<vector<double>> A, const vector<vector<double>> B) {
    int m = A.size(), n = A[0].size(), p = B[0].size();
    if (n != B.size()) {
       	printf("Khong nhan duoc\n");
        return {};
    }

    vector<vector<double>> C(m, vector<double>(p));
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

vector<vector<double>> Inv(vector<vector<double>> &input_matrix){
 signed int size = input_matrix.size();
    signed int i = 0;
    vector<vector<double>> I = generate_identity(size);
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

vector<vector<double>> Trans(const vector<vector<double>>& A){
	int m = A.size();
  	int n = A[0].size();

  	vector<vector<double>> C(n, vector<double>(m));

  	for (int i = 0; i < m; i++) {
    	for (int j = 0; j < n; j++) {
      	C[j][i] = A[i][j];
    	}
  	}

  	return C;
}

vector<double> matrixVectorMult(vector<vector<double>> A, vector<double> x) {
	int m = A.size();
	int n = A[0].size();

	if (n != x.size()){
		printf("Khong nhan duoc \n");
		return {};
	}
	vector<double> y(m);

	for (int i = 0; i < m; i++) {
		double sum = 0;
		for (int j = 0; j < n; j++) {
			sum += x[j] * A[i][j];
		}
		y[i] = sum;
	}
	
	return y;
}

void print_Matrix(vector<vector<double>> data,int r,int c){
	for (int i = 0; i < r; i++) {
    	for (int j = 0; j < c; j++) {
      		cout<< data[i][j] << " ";
    	}
    	cout<< endl;
	}
}


int main() {
	vector<vector<double>> data = readData();
    vector<vector<double>> y(data.size(), vector<double>(1));
    makeData(data,y);
    
	auto start_serial = chrono::high_resolution_clock::now();
    vector<vector<double>> data_trans = Trans(data);
    
    
	vector<vector<double>> temp1 = Mul(data_trans,data);
   
    
    
	vector<vector<double>> temp2 = Inv(temp1);
	
	
	vector<vector<double>> temp3 = Mul(data_trans, y);


	vector<vector<double>> temp4 = Mul(temp2, temp3);
	
	
	
	auto end_serial = chrono::high_resolution_clock::now();
	auto duration_serial = chrono::duration_cast<chrono::milliseconds>(end_serial - start_serial);
	cout<< duration_serial.count() <<endl;
	print_Matrix(temp4,temp4.size(),temp4[0].size());
	return 0;
}






