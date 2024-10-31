package Assignment2;

import java.util.Arrays;
import java.util.Scanner;

/**
 * Main class that handles user input and applies the Interior-Point algorithm for solving
 * linear programming problems.
 */
public class InteriorPointAlgorithm {

    /**
     * The main method to read user input and start the Interior-Point algorithm.
     *
     * @param args command-line arguments (not used)
     */
    public static void main(String[] args) {
        // Input
        Scanner scanner = new Scanner(System.in);
        // Number of variables that are involved in the equation
        int numVars = scanner.nextInt();
        // Coefficients of the objective function: Maximize 9x1 + 10x2 + 16x3 + 0x4 + 0x5 + 0x6
        double[] c = new double[numVars];
        for (int i = 0; i < numVars; i++) {
            c[i] = scanner.nextDouble(); // The coefficients for additional variables are zero
        }

        // Number of constraints
        int constraints = scanner.nextInt();
        // Constraint coefficients (including additional variables)
        //{ {18, 15, 12, 1, 0, 0},
        //  {6, 4, 8, 0, 1, 0},
        //  {5, 3, 3, 0, 0, 1}
        // }
        double[][] A = new double[constraints][numVars];
        for (int i = 0; i < constraints; i++) {
            for (int j = 0; j < numVars; j++) {
                A[i][j] = scanner.nextDouble();
            }
        }

        // Right-hand side of the constraints
        double[] b = new double[constraints];
        for (int i = 0; i < constraints; i++) {
            b[i] = scanner.nextDouble();
        }

        // Initial solution x0 = (x1, x2, x3, x4, x5, x6) = (1, 1, 1, 315, 174, 169)
        double[] x = new double[numVars];
        for (int i = 0; i < numVars; i++) {
            x[i] = scanner.nextDouble();
        }

        // Result of the algorithm for alpha = 0.5
        algorithm(0.5, x, A, c, b);
        System.out.println();
        // Result of the algorithm for alpha = 0.9
        algorithm(0.9, x, A, c, b);

    }

    /**
     * Method that represents all interior point algorithm
     *
     * @param alpha measures the fraction used of the distance that could be loved before the feasible region is left
     * @param x     solutions for equation on each iteration
     * @param A     the matrix representing the constraints' coefficients
     * @param c     the vector representing the objective function's coefficients
     * @param b     the vector representing the right-hand side of the constraints
     */
    public static void algorithm(double alpha, double[] x, double[][] A, double[] c, double[] b) {
        int iteration = 1;

        System.out.println("For alpha = " + alpha);

        while (true) {
            double[] v = Arrays.copyOf(x, x.length); // Save a copy x
            Matrix D = Matrix.diag(x);               // Create a diagonal matrix from x
            Matrix AA = new Matrix(A).dot(D);        // A * D
            Vector cc = D.dot(new Vector(c));        // D * c
            Matrix I = Matrix.eye(c.length);         // Identity matrix
            Matrix F = AA.dot(AA.transpose());       // AA * AA^T
            Matrix FI = F.inverse();                 // Inverse matrix F
            Matrix H = AA.transpose().dot(FI);       // AA^T * FI
            Matrix P = I.subtract(H.dot(AA));        // P = I - H * AA
            Vector cp = P.dot(cc);                   // cp = P * cc

            double nu = Math.abs(cp.min());          // Compute nu
            if (nu == 0) {
                System.out.println("Warning: nu is zero, cp: " + cp);
                break; // Обработка случая с нулевым nu
            }

            // y = 1 + (alpha / nu) * cp
            Vector y = Vector.ones(cp.size()).add(cp.multiply(alpha / nu));
            Vector yy = D.dot(y);                      // yy = D * y
            x = yy.toArray();                          // Update x

            if (iteration <= 4) {
                System.out.println("In iteration " + iteration + " we have x = " + Arrays.toString(x) + "\n");
            }
            iteration++;

            // Условие завершения
            if (yy.subtract(new Vector(v)).norm(2) < 0.0001) {
                break;
            }
        }

        // Objective function value
        double obj = 0.0;
        for (int i = 0; i < c.length; i++) {
            obj += c[i] * x[i];
        }

        // Results
        System.out.println("In the last iteration we have x = " + Arrays.toString(x) + "\n");
        System.out.println("Objective function value: " + obj);
    }
}

/**
 * Matrix class representing a mathematical matrix.
 */
class Matrix {
    private int rows;
    private int cols;
    private double[][] values;

    /**
     * Constructor for creating a matrix with specified dimensions.
     *
     * @param rows number of rows
     * @param cols number of columns
     */
    public Matrix(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.values = new double[rows][cols];
    }

    /**
     * Constructor for creating a matrix with specified dimensions and values.
     *
     * @param values 2D array containing matrix elements
     */
    public Matrix(double[][] values) {
        this.rows = values.length;
        this.cols = values[0].length;
        this.values = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            this.values[i] = Arrays.copyOf(values[i], cols);
        }
    }

    /**
     * Returns the value of the matrix at a specific position.
     *
     * @param i row index
     * @param j column index
     * @return the value at the specified position
     */
    public double getValue(int i, int j) {
        return values[i][j];
    }

    /**
     * Returns the value of the matrix at a specific position.
     *
     * @param i row index
     * @param j column index
     * @param value new value
     */
    public void setValue(int i, int j, double value) {
        values[i][j] = value;
    }

    /**
     * Returns the number of rows in the matrix.
     *
     * @return the number of rows
     */
    public int getNumRows() {
        return rows;
    }

    /**
     * Returns the number of columns in the matrix.
     *
     * @return the number of columns
     */
    public int getNumCols() {
        return cols;
    }

    /**
     * Returns the transpose matrix.
     *
     * @return the transpose matrix
     */
    public Matrix transpose() {
        double[][] transposed = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = values[i][j];
            }
        }
        return new Matrix(transposed);
    }

    /**
     * Returns the inverse matrix.
     *
     * @return the inverse matrix
     */
    public Matrix inverse() {
        if (rows != cols) {
            throw new IllegalArgumentException("Matrix must be square.");
        }

        double[][] augmentedMatrix = new double[rows][2 * rows];

        // Creation [A | I]
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                augmentedMatrix[i][j] = values[i][j];
            }
            augmentedMatrix[i][i + rows] = 1;
        }

        for (int i = 0; i < rows; i++) {
            double maxEl = Math.abs(augmentedMatrix[i][i]);
            int maxRow = i;
            for (int k = i + 1; k < rows; k++) {
                if (Math.abs(augmentedMatrix[k][i]) > maxEl) {
                    maxEl = Math.abs(augmentedMatrix[k][i]);
                    maxRow = k;
                }
            }

            double[] temp = augmentedMatrix[maxRow];
            augmentedMatrix[maxRow] = augmentedMatrix[i];
            augmentedMatrix[i] = temp;

            for (int k = i + 1; k < rows; k++) {
                double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
                for (int j = 0; j < 2 * rows; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }

        for (int i = rows - 1; i >= 0; i--) {
            double leadingElement = augmentedMatrix[i][i];
            for (int j = 0; j < 2 * rows; j++) {
                augmentedMatrix[i][j] /= leadingElement;
            }
            for (int k = i - 1; k >= 0; k--) {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * rows; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }

        double[][] inverseMatrix = new double[rows][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                inverseMatrix[i][j] = augmentedMatrix[i][j + rows];
            }
        }

        return new Matrix(inverseMatrix);
    }

    /**
     * Returns the multiplication AB, where A is this matrix.
     *
     * @param  B the second matrix in multiplication
     * @return the multiplication of two matrices
     */
    public Matrix dot(Matrix B) {
        if (cols != B.getNumRows()) {
            throw new IllegalArgumentException("Matrix dimensions do not match for multiplication.");
        }

        double[][] result = new double[rows][B.getNumCols()];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < B.getNumCols(); j++) {
                for (int k = 0; k < cols; k++) {
                    result[i][j] += values[i][k] * B.getValue(k, j);
                }
            }
        }
        return new Matrix(result);
    }

    /**
     * Returns the multiplication matrix * vector.
     *
     * @param vector for multiplication
     * @return       the multiplication of matrix and vector
     */
    public Vector dot(Vector vector) {
        if (cols != vector.size()) {
            throw new IllegalArgumentException("Matrix and vector dimensions do not match.");
        }

        double[] result = new double[rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i] += values[i][j] * vector.get(j);
            }
        }
        return new Vector(result);
    }

    /**
     * Returns the diagonal matrix from vector b.
     *
     * @param values values of vector b (right-hand side constraints)
     * @return       the diagonal matrix
     */
    public static Matrix diag(double[] values) {
        int n = values.length;
        double[][] diagMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            diagMatrix[i][i] = values[i];
        }
        return new Matrix(diagMatrix);
    }

    /**
     * Returns the identity matrix.
     *
     * @param size size of identity matrix
     * @return       the identity matrix
     */
    public static Matrix eye(int size) {
        double[][] identity = new double[size][size];
        for (int i = 0; i < size; i++) {
            identity[i][i] = 1.0;
        }
        return new Matrix(identity);
    }

    /**
     * Returns the subtraction of 2 matrices.
     *
     * @param other the second matrix
     * @return      subtraction of 2 matrices
     */
    public Matrix subtract(Matrix other) {
        if (this.rows != other.getNumRows() || this.cols != other.getNumCols()) {
            throw new IllegalArgumentException("Matrix dimensions do not match for subtraction.");
        }

        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = values[i][j] - other.getValue(i, j);
            }
        }
        return new Matrix(result);
    }
}

/**
 * Vector class representing a mathematical vector.
 */
class Vector {
    private double[] values;

    /**
     * Constructor for creating a vector with specified size and values.
     *
     * @param values array containing vector elements
     */
    public Vector(double[] values) {
        this.values = Arrays.copyOf(values, values.length);
    }

    /**
     * Returns the value of the vector at a specific index.
     *
     * @param index index of the element to be returned
     * @return the value at the specified index
     */
    public double get(int index) {
        return values[index];
    }

    /**
     * Returns the size of the vector.
     *
     * @return the size of the vector
     */
    public int size() {
        return values.length;
    }

    /**
     * Returns the copy of vector.
     *
     * @return the copy of the vector
     */
    public Vector copy() {
        return new Vector(Arrays.copyOf(values, values.length));
    }

    /**
     * Returns the multiplication by scalar.
     *
     * @param scalar for multiplication
     * @return the multiplication vector by scalar
     */
    public Vector multiply(double scalar) {
        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] * scalar;
        }
        return new Vector(result);
    }

    /**
     * Returns the summation of two vectors.
     *
     * @param other vector for addition
     * @return the summation of two vectors
     */
    public Vector add(Vector other) {
        if (values.length != other.size()) {
            throw new IllegalArgumentException("Vector dimensions do not match for addition.");
        }

        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] + other.get(i);
        }
        return new Vector(result);
    }

    /**
     * Returns the summation of vector and scalar.
     *
     * @param scalar that is added to each value of vector
     * @return the summation of vector and scalar
     */
    public Vector add(double scalar) {
        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] + scalar;
        }
        return new Vector(result);
    }

    /**
     * Returns the multiplication of two vectors.
     *
     * @param other vector for multiplication
     * @return the multiplication of vectors
     */
    public Vector dot(Vector other) {
        if (values.length != other.size()) {
            throw new IllegalArgumentException("Vector dimensions do not match for dot product.");
        }

        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] * other.get(i);
        }
        return new Vector(result);
    }

    /**
     * Returns the normalization of vector.
     *
     * @param p vector for multiplication
     * @return the multiplication of vectors
     */
    public double norm(int p) {
        double sum = 0.0;
        for (double v : values) {
            sum += Math.pow(Math.abs(v), p);
        }
        return Math.pow(sum, 1.0 / p);
    }

    /**
     * Returns the minimum value of vector.
     *
     * @return the min value of vector
     */
    public double min() {
        return Arrays.stream(values).min().orElse(Double.NaN);
    }

    /**
     * Returns the identity vector (all values are 1).
     *
     * @return the identity vector
     */
    public static Vector ones(int size) {
        double[] onesArray = new double[size];
        Arrays.fill(onesArray, 1.0);
        return new Vector(onesArray);
    }

    /**
     * Returns the subtraction of two vectors.
     *
     * @param other vector that is subtracted to each value of vector
     * @return the subtraction of two vectors
     */
    public Vector subtract(Vector other) {
        if (values.length != other.size()) {
            throw new IllegalArgumentException("Vector dimensions do not match for subtraction.");
        }

        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] - other.get(i);
        }
        return new Vector(result);
    }

    /**
     * Returns the copy of array.
     *
     * @return the copy of array
     */
    public double[] toArray() {
        return Arrays.copyOf(values, values.length);
    }

    /**
     * Returns the string from values from array.
     *
     * @return the string from values from array
     */
    @Override
    public String toString() {
        return Arrays.toString(values);
    }
}