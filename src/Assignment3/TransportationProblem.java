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

        if (!isApplicable(Supply, Demand, Coefficients_Of_Costs)) {
            System.out.println("The method is not applicable!");
            return;
        }

        if (!isBalanced(Supply, Demand)) {
            System.out.println("The problem is not balanced!");
            return;
        }

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

    public static boolean isApplicable(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs) {
        if (Supply.size() != Coefficients_Of_Costs.getNumRows() || Demand.size() != Coefficients_Of_Costs.getNumCols()) {
            return false;
        }
        return true;
    }

    public static boolean isBalanced(Vector Supply, Vector Demand) {
        // If method is not balanced
        double sumSupply = Supply.getSumValues();
        double sumDemand = Demand.getSumValues();
        if (sumSupply != sumDemand) {
            return false;
        }
        return true;
    }

    public static void North_West_Corner_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs) {
        int n = Coefficients_Of_Costs.getNumRows();
        int m = Coefficients_Of_Costs.getNumCols();

        // Create copies for Supply, Demand and Coefficients
        Vector supplyCopy = new Vector(Supply.size());
        supplyCopy.copy(Supply);
        Vector demandCopy = new Vector(Demand.size());
        demandCopy.copy(Demand);
        Matrix costsCopy = new Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                costsCopy.setValue(i, j, Coefficients_Of_Costs.getValue(i, j));
            }
        }

        Matrix resultPath = new Matrix(n, m);
        double result = 0;
        int i = 0;
        int j = 0;
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];

        while (i < n && j < m) {
            double currentValue = Math.min(supplyCopy.get(i), supplyCopy.get(j));
            resultPath.setValue(i, j, currentValue);
            result += currentValue * costsCopy.getValue(i, j);
            supplyCopy.set(i, supplyCopy.get(i) - currentValue);
            demandCopy.set(j, demandCopy.get(j) - currentValue);
            if (supplyCopy.get(i) == 0) {
                i++;
            }
            if (demandCopy.get(j) == 0) {
                j++;
            }
        }
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < m; l++) {
                if (resultPath.getValue(k, l) == -1) {
                    resultPath.setValue(k, l, 0);
                }
                resultsCost[k][l] = DecimalFormat.getInstance().format(resultPath.getValue(k, l));
                resultsDemand[l] = DecimalFormat.getInstance().format(demandCopy.get(l));
            }
            resultsSupply[k] = DecimalFormat.getInstance().format(supplyCopy.get(k));
        }
        printTable(resultsSupply, resultsCost, resultsDemand);
        System.out.println("Optimal path is " + DecimalFormat.getInstance().format(result));

    }

    public static void Vogel_s_Approximation_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs){
        int n = Coefficients_Of_Costs.getNumRows();
        int m = Coefficients_Of_Costs.getNumCols();

        // Create copies for Supply, Demand and Coefficients
        Vector supplyCopy = new Vector(Supply.size());
        supplyCopy.copy(Supply);
        Vector demandCopy = new Vector(Demand.size());
        demandCopy.copy(Demand);
        Matrix costsCopy = new Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                costsCopy.setValue(i, j, Coefficients_Of_Costs.getValue(i, j));
            }
        }
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < m; j++) {
//                System.out.print(costsCopy.getValue(i, j) + " ");
//            }
//            System.out.println();
//        }

        // To fix result table and optimal value
        Matrix resultPath = new Matrix(n, m);
        double result = 0;
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];

//        int step= 0;
        while (true) {
//            step++;
//            System.out.println("Iteration " + step);
            // Searching the difference between two minimal elements in rows and columns
            double[] differenceRows = new double[supplyCopy.size()];
            double[] differenceCols = new double[demandCopy.size()];
            boolean notEnd = false;
            for (int i = 0; i < n; i++) {
                differenceRows[i] = costsCopy.sub_two_min_Elements_From_Row(i);
                if (differenceRows[i] != -1) {
                    notEnd = true;
                }
            }
            for (int i = 0; i < m; i++) {
                differenceCols[i] = costsCopy.sub_two_min_Elements_From_Column(i);
                if (differenceCols[i] != -1) {
                    notEnd = true;
                }
            }
//            System.out.println("Difference row:");
//            for (int i = 0; i < differenceRows.length; i++) {
//                System.out.print(differenceRows[i] + " ");
//            }
//            System.out.println();
//            System.out.println("Difference column:");
//            for (int i = 0; i < differenceCols.length; i++) {
//                System.out.print(differenceCols[i] + " ");
//            }
//            System.out.println();
//            System.out.println();
            if (!notEnd) {
                break;
            }
            // Searching maximum value in differences
            double maxElementInRow = -1;
            int maxIndexInRow = -1;
            double maxElementInCol = -1;
            int maxIndexInCol = -1;
            for (int i = 0; i < differenceRows.length; i++) {
                if (differenceRows[i] > maxElementInRow) {
                    maxElementInRow = differenceRows[i];
                    maxIndexInRow = i;
                }
            }
            for (int i = 0; i < differenceCols.length; i++) {
                if (differenceCols[i] > maxElementInCol) {
                    maxElementInCol = differenceCols[i];
                    maxIndexInCol = i;
                }
            }

            // Determine what we choose, row or column
            double maxElement;
            boolean isCol = false;
            if (maxElementInCol > maxElementInRow) {
                maxElement = maxElementInCol;
                isCol = true;
            } else {
                maxElement = maxElementInRow;
            }

            // Take values from this line
            double[] currentLine;
            if (isCol) {
                currentLine = costsCopy.getColumn(maxIndexInCol);
            } else {
                currentLine = costsCopy.getRow(maxIndexInRow);
            }

            // Find minimum value from chosen line
            double minElementInLine = Integer.MAX_VALUE;
            int miniIndexLine = -1;
            for (int i = 0; i < currentLine.length; i++) {
                if (minElementInLine > currentLine[i]) {
                    miniIndexLine = i;
                    minElementInLine = currentLine[i];
                }
            }
            if (isCol) {
                double currentValue = Math.min(supplyCopy.get(miniIndexLine), demandCopy.get(maxIndexInCol));
                resultPath.setValue(miniIndexLine, maxIndexInCol, currentValue);
                result += currentValue * Coefficients_Of_Costs.getValue(miniIndexLine, maxIndexInCol);
                supplyCopy.set(miniIndexLine, supplyCopy.get(miniIndexLine) - currentValue);
                demandCopy.set(maxIndexInCol, demandCopy.get(maxIndexInCol) - currentValue);
                if (supplyCopy.get(miniIndexLine) == 0) {
                    for (int i = 0; i < m; i++) {
                        costsCopy.setValue(miniIndexLine, i, Integer.MAX_VALUE);
                    }
                }
                if (demandCopy.get(maxIndexInCol) == 0) {
                    for (int i = 0; i < m; i++) {
                        costsCopy.setValue(i, maxIndexInCol, Integer.MAX_VALUE);
                    }
                }
            } else {
                double currentValue = Math.min(supplyCopy.get(maxIndexInRow), demandCopy.get(miniIndexLine));
                resultPath.setValue(maxIndexInRow, miniIndexLine, currentValue);
                result += currentValue * costsCopy.getValue(maxIndexInRow, miniIndexLine);
                supplyCopy.set(maxIndexInRow, supplyCopy.get(maxIndexInRow) - currentValue);
                demandCopy.set(miniIndexLine, demandCopy.get(miniIndexLine) - currentValue);
                if (supplyCopy.get(maxIndexInRow) == 0) {
                    for (int i = 0; i < m; i++) {
                        costsCopy.setValue(maxIndexInRow, i, Integer.MAX_VALUE);

                    }
                }
                if (demandCopy.get(miniIndexLine) == 0) {
                    for (int i = 0; i < m; i++) {
                        costsCopy.setValue(i, miniIndexLine, Integer.MAX_VALUE);
                    }
                }
            }
//            System.out.println("Coeff main:");
//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < m; j++) {
//                    System.out.print(costsCopy.getValue(i, j) + " ");
//                }
//                System.out.println();
//            }
//            System.out.println();
//            System.out.println("Result costs:");
//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < m; j++) {
//                    System.out.print(resultPath.getValue(i, j) + " ");
//                }
//                System.out.println();
//            }
//            System.out.println();
//            System.out.println("Supply costs:");
//            for (int i = 0; i < n; i++) {
//                System.out.print(supplyCopy.get(i) + " ");
//            }
//            System.out.println();
//            System.out.println("Demand costs:");
//            for (int i = 0; i < m; i++) {
//                System.out.print(demandCopy.get(i) + " ");
//            }
//            System.out.println();
//            System.out.println();

        }
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < m; l++) {
                if (resultPath.getValue(k, l) == -1) {
                    resultPath.setValue(k, l, 0);
                }
                resultsCost[k][l] = DecimalFormat.getInstance().format(resultPath.getValue(k, l));
                resultsDemand[l] = DecimalFormat.getInstance().format(demandCopy.get(l));
            }
            resultsSupply[k] = DecimalFormat.getInstance().format(supplyCopy.get(k));
        }
        printTable(resultsSupply, resultsCost, resultsDemand);
        System.out.println("Optimal path is " + DecimalFormat.getInstance().format(result));
    }

    public static void Russell_s_Approximation_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs){

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
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.values[i][j] = -1;
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

    public double[] getRow(int i) {
        return values[i];
    }

    public double[] getColumn(int i) {
//        System.out.println(i);
//        for (int j = 0; j < rows; j++) {
//            for (int k = 0; k < cols; k++) {
//                System.out.print(values[j][k] + " ");
//            }
//            System.out.println();
//        }
        double[] column = new double[rows];
        for (int row = 0; row < rows; row++) {
            column[row] = values[row][i];
//            System.out.println(column[row]);
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

    public double sub_two_min_Elements_From_Row(int row) {
        double[] Row = getRow(row);
        double[] rowLine = new double[Row.length];
        for (int i = 0; i < Row.length; i++) {
            rowLine[i] = Row[i];
        }
        Arrays.sort(rowLine);
        if (rowLine[0] == Integer.MAX_VALUE && rowLine[1] == Integer.MAX_VALUE) {
            return -1;
        }
        if (rowLine[1] == Integer.MAX_VALUE) {
            return rowLine[0];
        }
        if (rowLine[0] == Integer.MAX_VALUE) {
            return rowLine[1];
        }
        return Math.abs(rowLine[1]-rowLine[0]);
    }

    public double max_Elements_From_Row(int row, Matrix matrix) {
        int n = getNumRows();
        double[] Row = getRow(row);
        Arrays.sort(Row);
        return Row[n-1];
    }

    public double sub_two_min_Elements_From_Column(int column) {
        double[] Column = getColumn(column);
        double[] columnLine = new double[Column.length];
        for (int i = 0; i < Column.length; i++) {
            columnLine[i] = Column[i];
        }
        Arrays.sort(columnLine);
//        for (int i = 0; i < Column.length; i++) {
//            System.out.print(columnLine[i] + " ");
//        }
//        System.out.println();
        if (columnLine[0] == Integer.MAX_VALUE && columnLine[1] == Integer.MAX_VALUE) {
            return -1;
        }
        if (columnLine[1] == Integer.MAX_VALUE) {
            return columnLine[0];
        }
        if (columnLine[0] == Integer.MAX_VALUE) {
            return columnLine[1];
        }

//        System.out.println(columnLine[0] + " " + columnLine[1]);
        return Math.abs(columnLine[1]-columnLine[0]);
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

    public Vector(int n) { this.values = new double[n]; }

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

    public void copy(Vector vector) {
        for (int i = 0; i < vector.size(); i++) {
            values[i] = vector.get(i);
        }
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
