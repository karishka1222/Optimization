package Assignment3;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Scanner;

/**
 * Main class that handles user input and applies the Transportation problem for solving
 * linear programming problems.
 */
public class TransportationProblem {
    /**
     * The main method to read user input and start the Interior-Point algorithm.
     *
     * @param args command-line arguments (not used)
     */
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        // Input all data
        ArrayList<String[]> inputList = new ArrayList<>();
        // Input supply vector
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
        double[][] costArr = new double[inputList.getFirst().length][Demand.size()];
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

    }

    private static void printTable(String[] supply, String[][] costs, String[] demand) {
        int S = supply.length;
        int D = costs[0].length;

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
            double currentValue = Math.min(supplyCopy.get(i), demandCopy.get(j));
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

        // To fix result table and optimal value
        Matrix resultPath = new Matrix(n, m);
        double result = 0;
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];

        while (true) {
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
            boolean isCol = maxElementInCol > maxElementInRow;

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
            int curRow, curColumn;
            if (isCol) {
                curRow = miniIndexLine;
                curColumn = maxIndexInCol;
            } else {
                curRow = maxIndexInRow;
                curColumn = miniIndexLine;
            }
            double currentValue = Math.min(supplyCopy.get(curRow), demandCopy.get(curColumn));
            resultPath.setValue(curRow, curColumn, currentValue);
            result += currentValue * Coefficients_Of_Costs.getValue(curRow, curColumn);
            supplyCopy.set(curRow, supplyCopy.get(curRow) - currentValue);
            demandCopy.set(curColumn, demandCopy.get(curColumn) - currentValue);
            if (supplyCopy.get(curRow) == 0) {
                for (int i = 0; i < m; i++) {
                    costsCopy.setValue(curRow, i, Integer.MAX_VALUE);
                }
            }
            if (demandCopy.get(curColumn) == 0) {
                for (int i = 0; i < n; i++) {
                    costsCopy.setValue(i, curColumn, Integer.MAX_VALUE);
                }
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

    public static void Russell_s_Approximation_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs){
        Vector supplyCopy = new Vector(Supply.size());
        supplyCopy.copy(Supply);
        Vector demandCopy = new Vector(Demand.size());
        demandCopy.copy(Demand);
        int n = Coefficients_Of_Costs.getNumRows();
        int m = Coefficients_Of_Costs.getNumCols();
        Matrix resultPath = new Matrix(n, m);
        Matrix resultPathCopy = new Matrix(n, m);
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];
        double result = 0;
        Matrix costsCopy = new Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                costsCopy.setValue(i, j, Coefficients_Of_Costs.getValue(i, j));
            }
        }
        while (true){
            int x = 0, y = 0;
            double max = Double.MIN_VALUE;
            for (int i = 0; i < n; i++) {
                double row_max = costsCopy.max_Elements_From_Row(i, costsCopy);
                for (int j = 0; j < m; j++) {
                    double col_max = costsCopy.max_Elements_From_Column(j, costsCopy);
                    if (costsCopy.getValue(i,j) != 0 && max < Math.abs(costsCopy.getValue(i, j) - col_max - row_max)) {
                        x = i;
                        y = j;
                        max = Math.abs(costsCopy.getValue(i, j) - col_max - row_max);

                    }
                }
            }
            if (demandCopy.getSumValues() == 0){
                break;
            }
            if (resultPath.getValue(x, y) == -1) {
                resultPath.setValue(x, y, 0);
            }
            if ((demandCopy.get(y) - supplyCopy.get(x)) == 0) {
                resultPathCopy.setValue(x, y, resultPath.getValue(x, y) + costsCopy.getValue(x, y) * demandCopy.get(y));
                resultPath.setValue(x, y, Math.min(demandCopy.get(y),supplyCopy.get(x)));
                supplyCopy.set(x, 0);
                demandCopy.set(y, 0);
                for (int i = 0; i < n; i++) {
                    costsCopy.setValue(i, y, 0);
                }
            } else if ((demandCopy.get(y) - supplyCopy.get(x)) > 0) {
                resultPathCopy.setValue(x, y, resultPath.getValue(x, y) + costsCopy.getValue(x, y) * supplyCopy.get(x));
                resultPath.setValue(x, y, Math.min(demandCopy.get(y),supplyCopy.get(x)));
                demandCopy.set(y, demandCopy.get(y) - supplyCopy.get(x));
                supplyCopy.set(x, 0);
                costsCopy.setValue(x, y, 0);
            } else if ((demandCopy.get(y) - supplyCopy.get(x)) < 0) {
                resultPathCopy.setValue(x, y, resultPath.getValue(x, y) + costsCopy.getValue(x, y) * demandCopy.get(y));
                resultPath.setValue(x, y, Math.min(demandCopy.get(y),supplyCopy.get(x)));
                supplyCopy.set(x, supplyCopy.get(x) - demandCopy.get(y));
                demandCopy.set(y, 0);
                for (int i = 0; i < n; i++) {
                    costsCopy.setValue(i, y, 0);
                }
            }
            costsCopy.setValue(x, y, 0);
        }
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < m; l++) {
                if (resultPathCopy.getValue(k, l) == -1) {
                    resultPathCopy.setValue(k, l, 0);
                } else {
                    result += resultPathCopy.getValue(k, l);
                }
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

    public Matrix(double[][] values)  {
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

    public double[] getRow(int i) {
        return values[i];
    }

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
        int n = getNumCols();
        double[] Row = getRow(row);
        double[] newRow = new double[Row.length];
        for (int i = 0; i < Row.length; i++) {
            newRow[i] = Row[i];
        }
        Arrays.sort(newRow);
        return newRow[n-1];
    }

    public double sub_two_min_Elements_From_Column(int column) {
        double[] Column = getColumn(column);
        double[] columnLine = new double[Column.length];
        for (int i = 0; i < Column.length; i++) {
            columnLine[i] = Column[i];
        }
        Arrays.sort(columnLine);
        if (columnLine[0] == Integer.MAX_VALUE && columnLine[1] == Integer.MAX_VALUE) {
            return -1;
        }
        if (columnLine[1] == Integer.MAX_VALUE) {
            return columnLine[0];
        }
        if (columnLine[0] == Integer.MAX_VALUE) {
            return columnLine[1];
        }

        return Math.abs(columnLine[1]-columnLine[0]);
    }

    public double max_Elements_From_Column(int column, Matrix matrix) {
        int n = getNumRows();
        double[] Column = getColumn(column);
        double[] newColumn = new double[Column.length];
        for (int i = 0; i < Column.length; i++) {
            newColumn[i] = Column[i];
        }
        Arrays.sort(newColumn);
        return newColumn[n-1];
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

    public double getSumValues() {
        double sum = 0.0;
        for (double v : values) {
            sum += v;
        }
        return sum;
    }

    @Override
    public String toString() {
        return Arrays.toString(values);
    }
}
