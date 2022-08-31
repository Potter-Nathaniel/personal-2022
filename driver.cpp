//Driver Code for Linear Algebra Tools
//Nathaniel Potter, Version 1.0
#include "matrix.h"

#define MAX_STORE 25

int loadMatrices(std::string, int, Matrix[]);
void unloadMatrices(std::string, int, Matrix[]);
int getMenuChoice();
int getMatrixChoice(int, Matrix[], bool);
bool saveAs(int, Matrix[], Matrix);
void printAsBasis(double**&, int, int, bool);
void refinedChoice(int&, int, Matrix[], bool, Matrix&);
void pressContinue();

int main(int argc, char* argv[]){
    /* To Test Inverse
    Matrix box(5, 5); //random stuff

    box.getInput();

    Matrix box2 = box.invert();
    Matrix box3 = box2 * (box.findDet());
    
    std::cout << "\nMatrix A:\n";
    box.printMatrix();
    std::cout << "\nMatrix A(-1):\n";
    box2.printMatrix();
    std::cout << "\nMatrix A(-1) * det(A):\ndet(A) = " << box.findDet() << "\n";
    box3.printMatrix();
    std::cout << "\nVerify AA(-1): \n";
    (box * box2).printMatrix();
    std::cout << "\nVerify that A(-1)(-1) = A:\n";
    (box2.invert()).printMatrix();
    */
    
    /* Testing Matrix Multiplication
    Matrix bong(3, 8);
    std::cout << "\nMatrix A, 8 elements per row:\n";
    bong.getInput();

    Matrix burp(8, 2);
    std::cout << "\n\nMatrix B, 5 elements per row:\n";
    burp.getInput();

    std::cout << "\nMatrix A:\n";
    bong.printMatrix();
    std::cout << "\nMatrix B:\n";
    burp.printMatrix();
    std::cout << "\nMatrix B*A:\n";
    try{
        (burp * bong).printMatrix();
    } catch (std::invalid_argument e){
        std::cout << "\n" << e.what() << "\n";
    }
    std::cout << "\nMatrix A*B\n";
    (bong * burp).printMatrix();
    */

    /* Testing REF
    int m, n;
    std::cout << "Dimension 1 (m): ";
    std::cin >> m;
    std::cout << "Dimension 2 (n): ";
    std::cin >> n;

    Matrix box(m, n);

    box.getInput();
    std::cout << "\nMatrix in REF:\n";
    box.rowEchelon().first.printMatrix();
    std::cout << "\nMatrix in RREF:\n";
    box.reducedEchelon().printMatrix();
    
    return 0;
    
    */
    
    if (argc != 2){
        std::cout << "Correct Usage: ./linAlgTools fileName";
        return 0;
    }

    std::cout << "Welcome to the Linear Algebra Tools!\n";

    Matrix entries[MAX_STORE];
    int size = loadMatrices(argv[1], MAX_STORE, entries);
    int choice; //choosing the option
    int selection; //within option
    double tempVal;
    Matrix temp, copy; //two temporary variables

    double** basis ; //2D array to store the basis
    int dim, len;
    //dimensions for the basis array. vectors within the basis will be stored in the ROWS
    //note that dim is the dimension of the space spanned by basis; dim V.
    
    bool isRunning = true;
    while (isRunning){
        choice = getMenuChoice();
        switch (choice){
            case 0: //Simple Exit
                isRunning = false;
                exit(0); //directly terminate
                break;
            case 1: //Edit Matrix
                selection = getMatrixChoice(size, entries, true);
                temp.getInput();
                if (selection == 0){
                    entries[size] = temp;
                    size++;
                } else {
                    entries[selection - 1] = temp;
                }
                break;
            case 2: //Echelon Form
                selection = getMatrixChoice(size, entries, true);
                if (selection == 0){
                    copy.getInput();
                    temp = copy.rowEchelon().first;
                } else {
                    temp = entries[selection - 1].rowEchelon().first;
                }
                temp.printMatrix();
                if(saveAs(size, entries, temp)){size++;}
                break;
            case 3: //RREF
                selection = getMatrixChoice(size, entries, true);
                if (selection == 0){
                    copy.getInput();
                    temp = copy.reducedEchelon();
                } else {
                    temp = entries[selection - 1].reducedEchelon();
                }
                
                temp.printMatrix();
                if(saveAs(size, entries, temp)){size++;}
                break;
            case 4: //Matrix Multiplication
                selection = getMatrixChoice(size, entries, true); //first matrix
                if (selection == 0){
                    std::cout << "Enter Matrix 1:\n";
                    temp.getInput();
                } else {
                    temp = entries[selection - 1];
                }

                selection = getMatrixChoice(size, entries, true); //second matrix
                if (selection == 0){
                    std::cout << "Enter Matrix 2:\n";
                    copy.getInput();
                } else {
                    copy = entries[selection - 1];
                }

                (temp * copy).printMatrix();
                if(saveAs(size, entries, temp * copy)){size++;}
                break;
            case 5:
                selection = getMatrixChoice(size, entries, false);
                entries[selection - 1].printMatrix();
                break;
            case 6:
                std::cout << "Matrix MUST be SQUARE.\n";

                refinedChoice(selection, size, entries, true, temp);

                try {
                    copy = temp.invert();
                    std::cout << "Inverted Matrix:\n";
                    copy.printMatrix();
                    if(saveAs(size, entries, copy)){size++;}

                    std::cout << "[A]^(-1) * det([A])\n";
                    (copy * (temp.findDet())).printMatrix();
                    if(saveAs(size, entries, copy * (temp.findDet()))){size++;}
                } catch (std::exception& e){
                    std::cout << e.what();
                }
                break;
            case 7:
                std::cout << "Matrix MUST be SQUARE.\n";
                refinedChoice(selection, size, entries, true, temp);

                try {
                    tempVal = temp.findDet();
                    std::cout << "The Determinant is: " << tempVal << '\n';
                } catch (std::exception& e){
                    std::cout << e.what();
                }
                break;
            case 8:
                std::cout << "Coming Soon!\n";
                break;
            case 9:
                std::cout << "Null Space Coming Soon!\n";
                break;
            case 10:
                refinedChoice(selection, size, entries, true, temp);

                basis = temp.columnSpace();
                dim = temp.rank();
                len = temp.getRows();

                std::cout << "\nColumn Space of Matrix:\n";
                std::cout << "Given as a Column-Vector Basis:\n";
                
                printAsBasis(basis, dim, len, true);
                break;
            case 11: //Row Space not confirmed yet.
                refinedChoice(selection, size, entries, true, temp);

                basis = temp.columnSpace();
                dim = temp.rank();
                len = temp.getCols();

                std::cout << "\nRow Space of Matrix:\n";
                std::cout << "Given as a Row-Vector Basis:\n";

                printAsBasis(basis, dim, len, false);
                break;
            case 12:
                refinedChoice(selection, size, entries, true, temp);

                std::cout << "The Rank of the matrix is: ";
                std::cout << temp.rank() << '\n';
                break;
            case 13:
                size = 0;
                std::cout << "\n===+ALL MATRICES WIPED+===\n";
                break;
            case 14:
                unloadMatrices(argv[1], size, entries);
                std::cout << "\n===+PROGRESS SAVED+===\n";
                break;
        }
        pressContinue();
    }

    unloadMatrices(argv[1], size, entries);
    return 0;
} 

/*
1 <- # of matrices contained
3 5 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15(\n)
(another matrix)

First two numbers determine dimension, rest is content.
[1 2 3 4 5]
[6 7 8 9 10]
[11 12 13 14 15]

*/

int loadMatrices(std::string fileName, int size, Matrix entries[]){ //failsafe later
    std::ifstream myFile(fileName);
    int count = 0;
    if (!myFile.is_open()){
        std::cout << "Sorry, could not open given file for reading.";
        return 0;
    }
    int head = 0;
    myFile >> head;
    while (count < head){
        int m, n;
        double val;
        myFile >> m;
        myFile >> n;
        Matrix temp(m, n);
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                myFile >> val;
                temp.setElement(i, j, val);
            }
        }
        entries[count] = temp;
        count++;
    }
    myFile.close();

    return count;
}

void unloadMatrices(std::string fileName, int size, Matrix entries[]){
    std::ofstream myFile(fileName); //this will REWRITE the file FROM SCRATCH
    myFile << size << '\n';
    int count = 0;
    for (int i = 0; i < size; i++){
        int m = entries[i].getRows();
        int n = entries[i].getCols();
        myFile << m << ' ';
        myFile << n << ' ';
        for (int j = 0; j < m; j++){
            for (int k = 0; k < n; k++){
                myFile << entries[i].readElement(j, k) << ' ';
            }
        }
        if (!(i == size - 1)){ //cannot have any spaces or \n's at the end
            myFile << '\n';
        }
        count++;
    }
    myFile.close();
}

int getMenuChoice(){
    std::cout << "\nWhat would you like to do?\n";
    std::cout << "0. Exit\n";
    std::cout << "1. Edit Matrix\n";
    std::cout << "2. Echelon Form of Matrix\n";
    std::cout << "3. RREF of a Matrix\n";
    std::cout << "4. Multiply Two Matrices\n";
    std::cout << "5. Print Matrix\n";

    std::cout << "\nSquare Matrix Options\n";
    std::cout << "6. Invert Matrix\n";
    std::cout << "7. Find Determinant\n";
    std::cout << "8. Determine Eigenvalues\n";

    std::cout << "\nVector Space Options\n";
    std::cout << "9. Find Null Space\n";
    std::cout << "10. Find Column Space\n";
    std::cout << "11. Find Row Space\n";
    std::cout << "12. Find Rank of Matrix\n";

    std::cout << "\nMemory Options\n";
    std::cout << "13. Clear ALL Memory\n";
    std::cout << "14. Save Progress\n";
    //std::cout << "15. Remove Matrix\n";

    int choice;
    std::cin >> choice; //add failsafes later
    return choice;  
}

int getMatrixChoice(int size, Matrix options[], bool enable){
    int choice;
    if (size == 0){
        choice = 0;
    } else {
        std::cout << "Select Matrix:\n";
        for (int i = 0; i < size; i++){
            std::cout << i+1 << ". ";
            std::cout << '(' << options[i].getRows() << 'x' << options[i].getCols() << ")\n";
        }
        if (enable){
            std::cout << "Enter 0 for new Matrix.\n";
        }

        std::cin >> choice;
    }
    return choice;
}

bool saveAs(int size, Matrix options[], Matrix toSave){
    bool save;
    std::cout << "Would you like to save this Matrix?\n";
    std::cout << "Enter 1 to save, 0 to discard: ";
    std::cin >> save;
    if (save){
        std::cout << "Saving to Matrix...\n";
        int spot = getMatrixChoice(size, options, true);
        if (spot == 0){
            options[size] = toSave;
        } else {
            options[spot - 1] = toSave;
            save = false;
        }
    }

    return save;
}

void printAsBasis(double** &basis, int dim, int length, bool mode){ //must fix
    if (mode) { //let column = 1 and row = 0
        for (int i = 0; i < length; i++){
            for (int j = 0; j < dim; j++){
                std::cout << '[';
                std::cout.width(8);
                std::cout << std::left << basis[j][i];
                std::cout << "] ";
            }
            std::cout << '\n';
        }
    }
    if (!mode){
        for (int i = 0; i < dim; i++){
            std::cout << "[";
            for (int j = 0; j < length; j++){
                std::cout.width(10);
                std::cout << std::left << basis[i][j];
            }
            std::cout << "]\n\n";
        }
    }

    for (int i = 0; i < dim; i++){ //2D Cleanup
        delete[] basis[i];
    }

    delete[] basis;
}

void refinedChoice(int& s, int si, Matrix e[], bool wr, Matrix& t){
    s = getMatrixChoice(si, e, wr);
    if (s == 0){
        t.getInput();
    } else {
        t = e[s - 1];
    }
}

void pressContinue(){
    std::cout << "\nEnter 0 to continue...\n";
    int pass = 1;
    std::cin >> pass;
}

//Notes
/*
Use std::cin.fail() to detect bad input or menu controls
File storage system?
*/
