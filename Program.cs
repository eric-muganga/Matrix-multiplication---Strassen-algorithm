
double[,] A =
{
    {1,2,3},
    {4,5,6}
};

double[,] B =
{
    {7,8},
    {9,10},
    {11,12}
};

try
{
    { 1, 2, 3 }, 
    for (int i = 0; i < C.GetLength(0); i++)
    {
        for (int j = 0; j < C.GetLength(1); j++) Console.Write(C[i, j] + " ");
        Console.WriteLine();
    }

double[,] MatrixMultiply(double[,] A, double[,] B)
{
    //check if the dimension match

    int rowsA = A.GetLength(0);
    int colsB = B.GetLength(1);
    if(rowsA != colsB)
    {
        throw new Exception("Error/Exception (Dimension mismatch)");
    }
   
    int pA = GetNearestPowerOfTwo(A);
    int pB = GetNearestPowerOfTwo(B);

    int n = Math.Max(pA, pB);

    A = MatrixPadding(A, n);
    B = MatrixPadding(B, n);

    double[,] product = new double[n, n];

    if (n == 1)
    {
        product[0, 0] = A[0, 0] * B[0, 0];
    }
    else
    {
        // Dividing matrices into four submatrices
        double[,] a11 = SplitMatrix(A, 0, 0);
        double[,] a12 = SplitMatrix(A, 0, n / 2);
        double[,] a21 = SplitMatrix(A, n / 2, 0);
        double[,] a22 = SplitMatrix(A, n / 2, n / 2);

        double[,] b11 = SplitMatrix(B, 0, 0);
        double[,] b12 = SplitMatrix(B, 0, n / 2);
        double[,] b21 = SplitMatrix(B, n / 2, 0);
        double[,] b22 = SplitMatrix(B, n / 2, n / 2);

        /** 
        M1 = (A11 + A22)(B11 + B22)
        M2 = (A21 + A22) B11
        M3 = A11 (B12 - B22)
        M4 = A22 (B21 - B11)
        M5 = (A11 + A12) B22
        M6 = (A21 - A11) (B11 + B12)
        M7 = (A12 - A22) (B21 + B22)
        **/
        double[,] m1 = MatrixMultiply(AddMatrices(a11, a22), AddMatrices(b11, b22));
        double[,] m2 = MatrixMultiply(AddMatrices(a21, a22), b11);
        double[,] m3 = MatrixMultiply(a11, SubtractMatrices(b12, b22));
        double[,] m4 = MatrixMultiply(a22, SubtractMatrices(b21, b11));
        double[,] m5 = MatrixMultiply(AddMatrices(a11, a12), b22);
        double[,] m6 = MatrixMultiply(SubtractMatrices(a21, a11), AddMatrices(b11, b12));
        double[,] m7 = MatrixMultiply(SubtractMatrices(a12, a22), AddMatrices(b21, b22));


        /**
        C11 = M1 + M4 - M5 + M7
        C12 = M3 + M5
        C21 = M2 + M4
        C22 = M1 - M2 + M3 + M6
        **/
        double[,] c11 = AddMatrices(SubtractMatrices(AddMatrices(m1, m4), m5), m7);
        double[,] c12 = AddMatrices(m3, m5);
        double[,] c21 = AddMatrices(m2, m4);
        double[,] c22 = AddMatrices(AddMatrices(SubtractMatrices(m1, m2), m3), m6);


        double[,] merged = MergeMatrices(c11, c12, c21, c22);
        product = RemovePadding(merged);
    }
    

    // spliting the matrix method
    double[,] SplitMatrix(double[,] matrix, int startRow, int startCol)
    {
        int rows = matrix.GetLength(0) / 2;
        double[,] newMatrix = new double[rows, rows];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < rows; j++)
            {
                newMatrix[i, j] = matrix[startRow + i, startCol + j];
            }
        }
        return newMatrix;
    }

    //Adding Matrices method
    double[,] AddMatrices(double[,] A, double[,] B)
    {
        int n = A.GetLength(0);
        double[,] addition = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                addition[i, j] = A[i, j] + B[i, j];
            }
        }

        return addition;
    }

    //Subtracting Matrices Method
    double[,] SubtractMatrices(double[,] A, double[,] B)
    {
        int n = A.GetLength(0);
        double[,] result = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result[i, j] = A[i, j] - B[i, j];
            }
        }

        return result;
    }

    //Merging the Matrices
    double[,] MergeMatrices(double[,] C11, double[,] C12, double[,] C21, double[,] C22)
    {
        int n = C11.GetLength(0) * 2;
        double[,] result = new double[n, n];

        for (int i = 0; i < n / 2; i++)
        {
            for (int j = 0; j < n / 2; j++)
            {
                result[i, j] = C11[i, j];
                result[i, j + n / 2] = C12[i, j];
                result[i + n / 2, j] = C21[i, j];
                result[i + n / 2, j + n / 2] = C22[i, j];
            }
        }
        return result;
    }

    double[,] MatrixPadding(double[,] matrix, int n)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        double[,] paddedMatrix = new double[n, n];

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                paddedMatrix[i, j] = matrix[i, j];
            }
        }

        return paddedMatrix;
    }

    int GetNearestPowerOfTwo(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        int max = Math.Max(rows, cols);

        int p = 1;
        while (p < max)
        {
            p *= 2;
        }
        return p;
    }


     double[,] RemovePadding(double[,] originalArray)
    {
        double[,] subarray = new double[rowsA, colsB];

        for (int i = 0; i < rowsA; i++)
        {
            for (int j = 0; j < colsB; j++)
            {
                subarray[i, j] = originalArray[i, j];
            }
        }

        return subarray;
    }


    return product;
}







