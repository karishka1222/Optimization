package Assignment2;

import java.util.Arrays;
import java.util.Scanner;

public class InteriorPointAlgorithm {

    public static void main(String[] args) {
        // Коэффициенты целевой функции: Максимизировать 9x1 + 10x2 + 16x3
        double[] c = {9, 10, 16, 0, 0, 0}; // Коэффициенты при дополнительных переменных равны нулю

        // Коэффициенты ограничений (включая дополнительные переменные)
        double[][] A = {
                {18, 15, 12, 1, 0, 0}, // Первое ограничение
                {6, 4, 8, 0, 1, 0},    // Второе ограничение
                {5, 3, 3, 0, 0, 1}     // Третье ограничение
        };


        // Правая часть ограничений
        double[] b = {360, 192, 180};

        // Начальное решение: (x1, x2, x3, x4, x5, x6) = (1, 1, 1, 315, 174, 169)
        double[] x = {1, 1, 1, 315, 174, 169};

//        Scanner inputData = new Scanner(System.in);
//        String action = inputData.nextLine();
//
//        if (action.equals("min")){
//            //коэффициенты функции нашей делаем противоположными по знаку, дальше всё также
//        }

        algorithm(0.5, x, A, c, b);
        System.out.println();
        algorithm(0.9, x, A, c, b);

    }

    public static void algorithm(double alpha, double[] x, double[][] A, double[] c, double[] b) {
        int iteration = 1;

        System.out.println("For alpha = " + alpha);

        while (true) {
            double[] v = Arrays.copyOf(x, x.length); // Сохранить копию x
            Matrix D = Matrix.diag(x);               // Создать диагональную матрицу из x
            Matrix AA = new Matrix(A).dot(D); // A * D
            Vector cc = D.dot(new Vector(c));         // D * c
            Matrix I = Matrix.eye(c.length);          // Единичная матрица
            Matrix F = AA.dot(AA.transpose());        // AA * AA^T
            Matrix FI = F.inverse();                  // Обратная матрица F
            Matrix H = AA.transpose().dot(FI);        // AA^T * FI
            Matrix P = I.subtract(H.dot(AA));         // P = I - H * AA
            Vector cp = P.dot(cc);                    // cp = P * cc

            double nu = Math.abs(cp.min());           // Вычислить nu
            if (nu == 0) {
                System.out.println("Warning: nu is zero, cp: " + cp);
                break; // Обработка случая с нулевым nu
            }

            // y = 1 + (alpha / nu) * cp
            Vector y = Vector.ones(cp.size()).add(cp.multiply(alpha / nu));
            Vector yy = D.dot(y);                      // yy = D * y
            x = yy.toArray();                          // Обновить x

            if (iteration <= 4) {
                System.out.println("In iteration " + iteration + " we have x = " + Arrays.toString(x) + "\n");
            }
            iteration++;

            // Условие завершения
            if (yy.subtract(new Vector(v)).norm(2) < 0.0001) {
                break;
            }
        }

        // results
        System.out.println("In the last iteration we have x = " + Arrays.toString(x) + "\n");
        System.out.println("Objective function value: " + (9 * x[0] + 10 * x[1] + 16 * x[2]));
    }
}

// class Matrix
class Matrix {
    private int rows;
    private int cols;
    private double[][] values;

    public Matrix(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.values = new double[rows][cols];
    }

    public Matrix(double[][] values) {
        this.rows = values.length;
        this.cols = values[0].length;
        this.values = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            this.values[i] = Arrays.copyOf(values[i], cols);
        }
    }

    public double getValue(int i, int j) {
        return values[i][j];
    }

    public void setValue(int i, int j, double value) {
        values[i][j] = value;
    }

    // number of row
    public int getNumRows() {
        return rows;
    }

    // number of column
    public int getNumCols() {
        return cols;
    }

    // T for Matrix
    public Matrix transpose() {
        double[][] transposed = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = values[i][j];
            }
        }
        return new Matrix(transposed);
    }

    // Inverse
    public Matrix inverse() {
        if (rows != cols) {
            throw new IllegalArgumentException("Matrix must be square.");
        }

        double[][] augmentedMatrix = new double[rows][2 * rows];

        // creation [A | I]
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

    // MUX for Matrix
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

    // MUX Matrix*Vector
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

    // D Matrix
    public static Matrix diag(double[] values) {
        int n = values.length;
        double[][] diagMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            diagMatrix[i][i] = values[i];
        }
        return new Matrix(diagMatrix);
    }

    // Identity
    public static Matrix eye(int size) {
        double[][] identity = new double[size][size];
        for (int i = 0; i < size; i++) {
            identity[i][i] = 1.0;
        }
        return new Matrix(identity);
    }

    // Subtraction
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

class Vector {
    private double[] values;

    public Vector(double[] values) {
        this.values = Arrays.copyOf(values, values.length);
    }

    public double get(int index) {
        return values[index];
    }

    public int size() {
        return values.length;
    }

    public Vector copy() {
        return new Vector(Arrays.copyOf(values, values.length));
    }

    // MUX
    public Vector multiply(double scalar) {
        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] * scalar;
        }
        return new Vector(result);
    }

    // Add
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

    public Vector add(double scalar) {
        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = values[i] + scalar;
        }
        return new Vector(result);
    }

    // MUV Vector*D
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

    public double norm(int p) {
        double sum = 0.0;
        for (double v : values) {
            sum += Math.pow(Math.abs(v), p);
        }
        return Math.pow(sum, 1.0 / p);
    }

    public double min() {
        return Arrays.stream(values).min().orElse(Double.NaN);
    }

    // Identity Vector
    public static Vector ones(int size) {
        double[] onesArray = new double[size];
        Arrays.fill(onesArray, 1.0);
        return new Vector(onesArray);
    }

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

    public double[] toArray() {
        return Arrays.copyOf(values, values.length);
    }

    @Override
    public String toString() {
        return Arrays.toString(values);
    }
}