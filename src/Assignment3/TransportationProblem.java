package Assignment3;
import java.util.Arrays;
import java.util.Scanner;

public class TransportationProblem {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);

        //example of input:
        /*
        3
        4 3 2
        3
        6 3 5
        4 20 7
        3 4 3
        4 3 2
        */
        // input of a vector of coefficients of supply - S.
        int number_Of_S_Coefficients = sc.nextInt();// по-идее у нас всегда 3
        double[] S_Coefficients = new double[number_Of_S_Coefficients];
        for (int i = 0; i < number_Of_S_Coefficients; i++) {
            S_Coefficients[i] = sc.nextDouble();
        }


        // input of a matrix of coefficients of costs - C
        int number_Of_D_Coefficients = sc.nextInt();// по-идее у нас всегда 4
        double[][] C_Coefficients = new double[number_Of_S_Coefficients][number_Of_D_Coefficients];
        for (int i = 0; i < number_Of_S_Coefficients; i++) {
            for (int j = 0; j < number_Of_D_Coefficients; j++) {
                C_Coefficients[i][j] = sc.nextDouble();
            }
        }


        // input of a vector of coefficients of demand - D
        double[] D_Coefficients = new double[number_Of_D_Coefficients];
        for (int i = 0; i < number_Of_D_Coefficients; i++) {
            D_Coefficients[i] = sc.nextDouble();
        }


        // creation of Vectors with Supply coefficients and Matrix with Coefficients_Of_Cost
        Vector Supply = new Vector(S_Coefficients);
        Vector Destination = new Vector(D_Coefficients);
        Matrix Coefficients_Of_Costs = new Matrix(C_Coefficients);

        North_West_Corner_Method(Supply, Destination, Coefficients_Of_Costs);

        Vogel_s_Approximation_Method(Supply, Destination, Coefficients_Of_Costs);

        Russell_s_Approximation_Method(Supply, Destination, Coefficients_Of_Costs);


        //New methods examples:

        // возвращает разность двух наименьших в строке
        double a = Matrix.sub_two_min_Elements_From_Row(1, Coefficients_Of_Costs);
        // возвращает max в строке
        double b = Matrix.max_Elements_From_Row(1, Coefficients_Of_Costs);


        // возвращает разность двух наименьших в столбце
        double c = Matrix.sub_two_min_Elements_From_Column(0, Coefficients_Of_Costs);
        // возвращает max в столбце
        double d = Matrix.max_Elements_From_Column(2, Coefficients_Of_Costs);

        System.out.println("a = " + a + ", b = " + b + ", c = " + c + ", d = " + d);
    }



    public static void North_West_Corner_Method(Vector Supply, Vector Destination, Matrix Coefficients_Of_Costs) {

    }

    public static void Vogel_s_Approximation_Method(Vector Supply, Vector Destination, Matrix Coefficients_Of_Costs){

    }

    public static void Russell_s_Approximation_Method(Vector Supply, Vector Destination, Matrix Coefficients_Of_Costs){

    }

}


// class Matrix
class Matrix {
    private static int rows;
    private static int cols;
    private static double[][] values;

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

    public static double[] getRow(int i) {return Arrays.copyOf(values[i], rows);}

    public static double[] getColumn(int i) {
        double[] column = new double[rows];
        for (int row = 0; row < rows; row++) {
            column[row] = values[row][i];
        }
        return column;
    }

    public void setValue(int i, int j, double value) {
        values[i][j] = value;
    }

    // number of row
    public static int getNumRows() {
        return rows;
    }
    // number of column
    public static int getNumCols() {
        return cols;
    }

    public static double sub_two_min_Elements_From_Row(int row, Matrix matrix) {
        double[] Row = getRow(row);
        Arrays.sort(Row);
        return Math.abs(Row[1]-Row[0]);
    }

    public static double max_Elements_From_Row(int row, Matrix matrix) {
        int n = getNumRows();
        double[] Row = getRow(row);
        Arrays.sort(Row);
        return Row[n-1];
    }


    public static double sub_two_min_Elements_From_Column(int column, Matrix matrix) {
        double[] Column = getColumn(column);
        Arrays.sort(Column);
        return Math.abs(Column[1]-Column[0]);
    }

    public static double max_Elements_From_Column(int column, Matrix matrix) {
        int n = getNumCols();
        double[] Column = getColumn(column);
        Arrays.sort(Column);
        return Column[n-1];
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