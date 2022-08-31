//Function Definitions for the Matrix Class

#include "matrix.h"

//Constructors
Matrix::Matrix(){ //initialized to 1x1, shouldn't be used too much, though
    rows = 1;
    cols = 1;
    mat[0][0] = 0.0;
}

Matrix::Matrix(int r, int c){ //this should be the go-to constructor
    rows = r;
    cols = c;
    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            mat[i][j] = 0.0;
        }
    }
}

Matrix::Matrix(const Matrix& old){ //copy constructor
    rows = old.rows;
    cols = old.cols;
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            mat[i][j] = old.mat[i][j];
        }
    }
}

//getters
int Matrix::getRows(){ return rows; }
int Matrix::getCols(){ return cols; }

//setters
void Matrix::setRows(int r){
    int oldRows = rows;
    rows = r;
    if (oldRows < rows){
        for (int i = oldRows; i < rows; i++){
            for (int j = 0; j < cols; j++){
                mat[i][j] = 0.0;
            }
        }
    }
    if (oldRows > rows){
        for (int i = oldRows; i > rows; i--){
            for (int j = 0; j < cols; j++){
                mat[i][j] = 0.0;
            }
        }
    }
}

void Matrix::setCols(int c){
    int oldCols = cols;
    cols = c;
    if (oldCols < cols){
        for (int i = 0; i < rows; i++){
            for (int j = oldCols; j < cols; j++){
                mat[i][j] = 0.0;
            }
        }
    }
    if (oldCols > cols){
        for (int i = 0; i < rows; i++){
            for (int j = oldCols; j > cols; j--){
                mat[i][j] = 0.0;
            }
        }
    }
}

//general matrix stuff
void Matrix::setElement(int r, int c, double val){
    mat[r][c] = val;
}

double Matrix::readElement(int r, int c){
    return mat[r][c];
}

void Matrix::printMatrix(){ //do output formatting
    for (int i = 0; i < rows; i++){
        std::cout << "[";
        for (int j = 0; j < cols; j++){
            std::cout.width(10); std::cout << std::left << mat[i][j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "\n";
}

void Matrix::swapRows(int r1, int r2){
    double tempRow[cols];
    for (int i = 0; i < cols; i++){
        tempRow[i] = mat[r2][i];
    }
    for (int j = 0; j < cols; j++){
        mat[r2][j] = mat[r1][j];
    }
    for (int k = 0; k < cols; k++){
        mat[r1][k] = tempRow[k];
    }
}

Matrix Matrix::reducedEchelon(){
    Matrix result = ((*this).rowEchelon()).first;
    double ratio;
    
    for (int i = rows - 1; i > -1; i--){ //start at last row
        for (int j = 0; j < cols; j++){ //start at left side
            if (result.mat[i][j] != 0){ //if the current value is nonzero; a PIVOT VALUE
                ratio = (double) 1/(result.mat[i][j]); //calculate reciprocal
                result.scaleRow(i, ratio); //scale the row down such that the leading term is now 1
                result.mat[i][j] = 1; //set explicitly equal for precision
                for (int k = 0; k < i; k++){ //for every row above the current one
                    result.rowAddition(i, k, -(result.mat[k][j])); //add such that every value above the pivot is 0
                    result.mat[k][j] = 0; //explicitly set equal to 0 for precision. 
                }
                break;
            }
        }
    }
    result.roundDown();

    return result;
}

std::pair<Matrix, bool> Matrix::rowEchelon(){ 
    Matrix result = *this;
    bool swapFlag;
    int r = 0; //coords of current pivot
    int c = 0;
    double ratio;

    while (r < rows && c < cols){ //while within the matrix
        if (result.mat[r][c] == 0){ //if the current pivot is 0...
            bool colIsZero = true; //make bool that stores if the column below the pivot is filled with 0
            for (int j = r+1; j < rows; j++){ //search through remaining rows for nonzero leading term
                if (result.mat[j][c] != 0){ //if found, swap the rows and break out of the loop. 
                    result.swapRows(r, j);
                    swapFlag = true;
                    colIsZero = false;
                    break;
                } 
            }
            if (colIsZero){ //if the column is filled with zeroes below the pivot
                c++; //move to the next column. 
            }
        } 
        
        for (int i = r+1; i < rows; i++){ //for every row after the pivot
            ratio = (result.mat[i][c]) / (result.mat[r][c]); //calculate ratio between first terms
            result.rowAddition(r, i, -ratio); //subtract with ratio, eliminating first term of row below
            result.mat[i][c] = 0; //set column directly below to exactly 0
        }
        
        r++; c++; //move to the next pivot coordinates        
    }
    result.roundDown();
    return std::make_pair(result, swapFlag);
}

void Matrix::scaleRow(int r, double s){
    for (int i = 0; i < cols; i++){
        mat[r][i] *= s;
    }
}

void Matrix::scaleCol(int c, double s){
    for (int i = 0; i < rows; i++){
        mat[i][c] *= s;
    }
}

void Matrix::rowAddition(int start, int dest, double s){
    for (int i = 0; i < cols; i++){
        mat[dest][i] += mat[start][i] * s;
    }
}

Matrix Matrix::operator+(const Matrix& m){
    Matrix result(rows, cols);
    if (rows == m.rows && cols == m.cols){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                double tempSum = m.mat[i][j] + mat[i][j];
                result.setElement(i, j, tempSum);
            }
        }
    } else {
        throw std::invalid_argument("Attempted to Add Incompatible Matrix Sizes.\n");
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& m){ 
    if (cols != m.rows){ //check for mxn * nxp matrix multiplication
        throw std::invalid_argument("Invalid Matrix Dimensions for Multiplication.\n");
    }
    Matrix result(rows, m.cols); //initialize result as mxp matrix

    int size = cols;
    double temp[size];
    double temp2[size];

    for (int i = 0; i < rows; i++){ //current row
        for (int j = 0; j < m.cols; j++){ //current column
            for (int k = 0; k < cols; k++){ //copy current row of mat
                temp[k] = mat[i][k];
            }
            for (int l = 0; l < m.rows; l++){ //copy current column of m
                temp2[l] = m.mat[l][j];
            }
            result.setElement(i, j, dotProduct(size, temp, temp2));
        }
    }
    result.roundDown();

    return result;
}

Matrix Matrix::operator*(const double s){
    Matrix result = *this;
    for (int i = 0; i < rows; i++){
        result.scaleRow(i, s);
    }
    return result;
}

std::pair<int, int>* Matrix::pivotIndices(){ //Must Test
    Matrix result = (this->rowEchelon()).first;
    int count = 0;
    std::pair<int, int> myCoords[MAX_DIM];
    
    for (int i = rows - 1; i > -1; i--){ //start at last row
        for (int j = 0; j < cols; j++){ //start at left side

            //std::cout << '\n' << '(' << i << ',' << j << ')'; 

            if (result.mat[i][j] != 0){ //if the current value is nonzero; a PIVOT VALUE
                myCoords[count] = std::make_pair(i, j);
                //std::cout << "PIVOT";
                count++;
                break;
            }
        }
    }

    std::pair<int, int>* pivotList = new std::pair<int, int>[count + 1];
    pivotList[0] = std::make_pair(count, 0); //jank way of storing the rank
    for (int k = 0; k < count; k++){
        pivotList[k + 1] = myCoords[count - k - 1];
    }


    return pivotList;
}

double** Matrix::columnSpace(){ //pivot columns in ROWS of result
    std::pair<int, int>* pivots = this->pivotIndices();
    int rank = pivots[0].first;

    double** result = new double*[rank];

    for (int i = 1; i <= rank; i++){
        result[i - 1] = new double[rows];
        int index = pivots[i].second;
        for (int j = 0; j < rows; j++){
            result[i - 1][j] = mat[j][index];
        }
    }

    delete[] pivots; //cleanup

    return result;
}

//hasnt been tested
double** Matrix::rowSpace(){ //pivot rows stored in rows.
    std::pair<int, int>* pivots = this->pivotIndices();
    int rank = pivots[0].first;

    double** result = new double*[rank];

    for (int i = 1; i <= rank; i++){
        result[i-1] = new double[cols];
        int index = pivots[i].first;
        for (int j = 0; j < cols; j++){
            result[i-1][j] = abs(mat[j][index]); //test

        }
    }
    
    delete[] pivots;

    return result;
}

int Matrix::rank(){
    std::pair<int, int>* pivots = this->pivotIndices();
    int rank = pivots[0].first;
    delete[] pivots;

    return rank;
}

//Square Matrix Exclusive
double Matrix::findDet(){
    double det = 0.0;
    if (rows != cols){ //check if valid matrix dimensions
        throw std::invalid_argument("Non-Square Matrix: Cannot Find Determinant.\n");
    }
    if (rows == 1){ //base case
        det = mat[0][0];
    }
    for (int i = 0; i < cols; i++){
        Matrix temp(rows - 1, cols - 1);
        for (int j = 1; j < rows; j++){
            for (int k = 0; k < cols; k++){
                if (k < i){
                    temp.mat[j-1][k] = mat[j][k];
                } else if (k > i){
                    temp.mat[j-1][k-1] = mat[j][k];
                }
            }
        }
        det += mat[0][i] * pow(-1, i) * temp.findDet();
    }
    return det;
}

Matrix Matrix::invert(){
    if (rows != cols){
        throw std::invalid_argument("Cannot Invert Non-Square Matrix.\n");
    }
    Matrix temp = (*this);
    if (temp.findDet() == 0){
        throw std::runtime_error("Cannot Invert Singular Matrix.\n");
    }
    temp.setCols(2 * cols); //DOUBLE the array size
    for (int i = 0; i < rows; i++){
        temp.setElement(i, i + cols, 1);
    }

    Matrix temp2 = temp.reducedEchelon();

    Matrix result(rows, cols);
    double element;

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < rows; j++){
            element = temp2.mat[i][j + cols];
            result.setElement(i, j, element);
        }
    }
    
    return result;
}

double Matrix::trace(){
    if (rows != cols){
        throw std::invalid_argument("Cannot trace non-square matrix.\n");
    }
    double sum = 0;
    for (int i = 0; i < rows; i++){
        sum += mat[i][i];
    }
    return sum;
}

double Matrix::dotProduct(int n, double v1[], double v2[]){
    double result = 0; //implement size failsafe
    for (int i = 0; i < n; i++){
        result += v1[i] * v2[i];
    }
    return result;
}

void Matrix::roundDown(){
    for (int i = 0; i < rows; i++){ //round all small numbers down to 0. 
        for (int j = 0; j < cols; j++){
            if ((mat[i][j] < PRECISION) && (mat[i][j] > -PRECISION)){
                mat[i][j] = 0;
            }
        }
    }
}

void Matrix::getInput(){ //refine this later
    double val;
    int m, n;
    std::cout << "Set Rows: "; std::cin >> m;
    std::cout << "Set Columns: "; std::cin >> n;

    setRows(m); setCols(n);

    for (int i = 0; i < rows; i++){
        std::cout << "Enter row " << i << ": ";
        for (int j = 0; j < cols; j++){
            std::cin >> val;
            mat[i][j] = val;
        }
    }
}


