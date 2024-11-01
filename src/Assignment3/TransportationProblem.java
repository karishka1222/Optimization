package Assignment3;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Scanner;

public class TransportationProblem {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);

        //example of input:
        /*
4 3 2
6 3 5
4 20 7
3 4 3
4 3 2
        */
        // input all data
        ArrayList<String[]> inputList = new ArrayList<>();
        // input supply vector
        inputList.add(sc.nextLine().split(" "));
        // input costs
        for (int i = 0; i < inputList.getFirst().length; i++) {
            inputList.add(sc.nextLine().split(" "));
        }
        // input demands
        inputList.add(sc.nextLine().split(" "));

        // process input
        // 1. Supply Vector S
        double[] supplyArr = new double[inputList.getFirst().length];
        for (int i = 0; i < supplyArr.length; i++) {
            supplyArr[i] = Double.parseDouble(inputList.getFirst()[i]);
        }
        Vector Supply = new Vector(supplyArr);

        //2. demand & destination Vector D
        double[] destArr = new double[inputList.getLast().length];
        for (int i = 0; i < destArr.length; i++) {
            destArr[i] = Double.parseDouble(inputList.getLast()[i]);
        }
        Vector Demand = new Vector(destArr);

        //3. Cost Matrix C
        double[][] costArr = new double[inputList.getFirst().length][inputList.size()-2];
        String[][] costs = new String[costArr.length][costArr[0].length];
        for (int row = 0; row < inputList.size()-2; row++) {
            String[] arr = inputList.get(row+1);
            for (int col = 0; col < arr.length; col++) {
                costArr[row][col] = Double.parseDouble(arr[col]);
                costs[row][col] = arr[col];
            }
        }
        Matrix Coefficients_Of_Costs = new Matrix(costArr);

        // demonstrate the input:
        System.out.println("Initial table:");
        printTable(inputList.getFirst(), costs, inputList.getLast());

        System.out.println();
        System.out.println("Result of North-West Corner Method:");
        North_West_Corner_Method(Supply, Demand, Coefficients_Of_Costs);

        System.out.println();
        System.out.println("Result of Vogel's Approximation Method:");
        Vogel_s_Approximation_Method(Supply, Demand, Coefficients_Of_Costs);

        System.out.println();
        System.out.println("Result of Russell's Approximation Method:");
        Russell_s_Approximation_Method(Supply, Demand, Coefficients_Of_Costs);

        //New methods examples:

        // возвращает разность двух наименьших в строке
//        double a = Matrix.sub_two_min_Elements_From_Row(1, Coefficients_Of_Costs);
        // возвращает max в строке
//        double b = Matrix.max_Elements_From_Row(1, Coefficients_Of_Costs);


        // возвращает разность двух наименьших в столбце
//        double c = Matrix.sub_two_min_Elements_From_Column(0, Coefficients_Of_Costs);
        // возвращает max в столбце
//        double d = Matrix.max_Elements_From_Column(2, Coefficients_Of_Costs);

//        System.out.println("a = " + a + ", b = " + b + ", c = " + c + ", d = " + d);
    }

    private static void printTable(String[] supply, String[][] costs, String[] demand) {
        int S = supply.length;
        int D = costs.length;

        // Create the header
        System.out.printf("%30s%n", "Destination");
        System.out.printf("%8s", " "); // Empty top left cell
        for (int i = 0; i < D; i++) {
            System.out.printf("%8d", i + 1);
        }
        System.out.printf("%8s%n", "Supply");

        // Print the cost rows
        for (int i = 0; i < S; i++) {
            System.out.printf("Source %d", i + 1);
            for (int j = 0; j < D; j++) {
                System.out.printf("%8s", costs[i][j]);
            }
            System.out.printf("%8s%n", supply[i]);
        }

        // Print the demand row
        System.out.printf("%8s", "Demand");
        for (String d : demand) {
            System.out.printf("%8s", d);
        }
        System.out.println();
    }

    public static void North_West_Corner_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs) {
        int n = Coefficients_Of_Costs.getNumRows();
        int m = Coefficients_Of_Costs.getNumCols();
        // If method is not applicable
        if (Supply.size() != n || Demand.size() != m) {
            System.out.println("The method is not applicable!");
            return;
        }
        // If method is not balanced
        double sumSupply = Supply.getSumValues();
        double sumDemand = Demand.getSumValues();
        if (sumSupply != sumDemand) {
            System.out.println("The problem is not balanced!");
            return;
        }

        Matrix resultPath = new Matrix(n, m);
        double result = 0;
        int i = 0;
        int j = 0;
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];

        while (i < n && j < m) {
            double currentValue = Math.min(Supply.get(i), Demand.get(j));
            resultPath.setValue(i, j, currentValue);
            result += currentValue * Coefficients_Of_Costs.getValue(i, j);
            Supply.set(i, Supply.get(i) - currentValue);
            Demand.set(j, Demand.get(j) - currentValue);
            if (Supply.get(i) == 0) {
                i++;
            }
            if (Demand.get(j) == 0) {
                j++;
            }
        }
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < m; l++) {
                resultsCost[k][l] = DecimalFormat.getInstance().format(resultPath.getValue(k, l));
                resultsDemand[l] = DecimalFormat.getInstance().format(Demand.get(l));
            }
            resultsSupply[k] = DecimalFormat.getInstance().format(Supply.get(k));
        }
        printTable(resultsSupply, resultsCost, resultsDemand);
        System.out.println("Optimal path is " + result);

    }

    public static void Vogel_s_Approximation_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs){

    }

    public static void Russell_s_Approximation_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs){

    }

}


// class Matrix
class Matrix {
    private int rows;
    private int cols;
    private static double[][] values;

    public Matrix(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.values = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.values[i][j] = 0.0;
            }
        }
    }

    public Matrix(double[][] values) {
        this.rows = values.length;
        this.cols = values[0].length;
        this.values = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            this.values[i] = Arrays.copyOf(values[i], cols);
        }
    }

    public double[][] getValues() { return values; }

    public double getValue(int i, int j) {
        return values[i][j];
    }

    public double[] getRow(int i) {return Arrays.copyOf(values[i], rows);}

    public double[] getColumn(int i) {
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
    public int getNumRows() {
        return rows;
    }
    // number of column
    public int getNumCols() {
        return cols;
    }

    public double sub_two_min_Elements_From_Row(int row, Matrix matrix) {
        double[] Row = getRow(row);
        Arrays.sort(Row);
        return Math.abs(Row[1]-Row[0]);
    }

    public double max_Elements_From_Row(int row, Matrix matrix) {
        int n = getNumRows();
        double[] Row = getRow(row);
        Arrays.sort(Row);
        return Row[n-1];
    }

    public double sub_two_min_Elements_From_Column(int column, Matrix matrix) {
        double[] Column = getColumn(column);
        Arrays.sort(Column);
        return Math.abs(Column[1]-Column[0]);
    }

    public double max_Elements_From_Column(int column, Matrix matrix) {
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

    public void set(int index, double value) { values[index] = value; }

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

    public double getSumValues() {
        double sum = 0.0;
        for (double v : values) {
            sum += v;
        }
        return sum;
    }

    public double[] toArray() {
        return Arrays.copyOf(values, values.length);
    }

    @Override
    public String toString() {
        return Arrays.toString(values);
    }
}
