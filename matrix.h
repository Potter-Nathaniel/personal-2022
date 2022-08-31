//Defintion for Class Matrix for Linear Algebra
//Nathaniel Potter, Started 6/29/2022

#ifndef MATRIX_H
#define MATRIX_H

//include libraries
#include <iostream>
#include <cmath> 
#include <utility>
#include <string>
#include <fstream>

#define MAX_DIM 64
#define PRECISION 0.000000000001 //1e-12

class Matrix{  
    int rows, cols;
    double mat[MAX_DIM][MAX_DIM];

    //internal methods
    double dotProduct(int, double[], double[]);

    public:
        //constructors
        Matrix();
        Matrix(int, int);
        Matrix(const Matrix&);

        //getters
        int getRows();
        int getCols();

        //setters, Matrix auto resizes
        void setRows(int); 
        void setCols(int);

        //general matrix stuff
        void setElement(int, int, double);
        double readElement(int, int);
        void printMatrix();
        void swapRows(int, int);
        std::pair<Matrix, bool> rowEchelon();
        Matrix reducedEchelon(); 
        void scaleRow(int, double);
        void scaleCol(int, double);
        void rowAddition(int, int, double);
        Matrix operator+(const Matrix&);
        Matrix operator*(const Matrix&); //Matrix multiplication
        Matrix operator*(const double);  //Scalar Multiplication

        //std::pair<Matrix, Matrix> factorizeLU(); probably won't implement
        //Add vector space tools, Nul A, Col A, Row A, rank A, Basis stuff, Eigenvalues

        std::pair<int, int>* pivotIndices();
        Matrix transform();
        double** columnSpace();
        double** rowSpace();
        int rank();

        double** nullSpace();
        double* eigenvalues(); //see power method.


        //exclusive to SQUARES
        double findDet();
        Matrix invert(); 
        double trace(); 
        
        //misc methods
        void roundDown();
        void getInput();
};  

#endif