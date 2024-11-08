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
        // Input costs
        for (int i = 0; i < inputList.getFirst().length; i++) {
            inputList.add(sc.nextLine().split(" "));
        }
        // Input demands
        inputList.add(sc.nextLine().split(" "));

        // Creating Supply and Demand vectors and Coefficient matrix
        // Supply Vector S
        double[] supplyArr = new double[inputList.getFirst().length];
        for (int i = 0; i < supplyArr.length; i++) {
            supplyArr[i] = Double.parseDouble(inputList.getFirst()[i]);
        }
        Vector Supply = new Vector(supplyArr);

        // Demand Vector
        double[] destArr = new double[inputList.getLast().length];
        for (int i = 0; i < destArr.length; i++) {
            destArr[i] = Double.parseDouble(inputList.getLast()[i]);
        }
        Vector Demand = new Vector(destArr);

        // Cost Matrix C
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


        // Output of the program
        System.out.println("Initial table:");
        printTable(inputList.getFirst(), costs, inputList.getLast());

        // Check applicability
        if (!isApplicable(Supply, Demand, Coefficients_Of_Costs)) {
            System.out.println("The method is not applicable!");
            return;
        }

        // Check balance
        if (!isBalanced(Supply, Demand)) {
            System.out.println("The problem is not balanced!");
            return;
        }

        // Result of North-West Corner Method
        System.out.println();
        System.out.println("Result of North-West Corner Method:");
        North_West_Corner_Method(Supply, Demand, Coefficients_Of_Costs);

        // Result of Vogel's Approximation Method
        System.out.println();
        System.out.println("Result of Vogel's Approximation Method:");
        Vogel_s_Approximation_Method(Supply, Demand, Coefficients_Of_Costs);

        // Result of Russell's Approximation Method
        System.out.println();
        System.out.println("Result of Russell's Approximation Method:");
        Russell_s_Approximation_Method(Supply, Demand, Coefficients_Of_Costs);
    }

    /**
     * Method that print table with supply, demand and cost values
     *
     * @param supply array with supply values
     * @param costs 2D array with cost values
     * @param demand array with demand values
     */
    private static void printTable(String[] supply, String[][] costs, String[] demand) {
        int S = supply.length; // supply size
        int D = costs[0].length; // demand size

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

    /**
     * Method that check that all 3 algorithms are applicable for this problem
     *
     * @param Supply vector with supply values
     * @param Demand vector with demand values
     * @param Coefficients_Of_Costs matrix with cost values
     */
    public static boolean isApplicable(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs) {
        // if num of supply values != num of rows in the table and num of demand values != num of columns in the table
        if (Supply.size() != Coefficients_Of_Costs.getNumRows() || Demand.size() != Coefficients_Of_Costs.getNumCols()) {
            return false;
        }
        return true;
    }

    /**
     * Method that check that all 3 algorithms are applicable for this problem
     *
     * @param Supply vector with supply values
     * @param Demand vector with demand values
     */
    public static boolean isBalanced(Vector Supply, Vector Demand) {
        // Sum of supply must be equal to sum of demand
        double sumSupply = Supply.getSumValues();
        double sumDemand = Demand.getSumValues();
        return sumSupply == sumDemand;
    }

    /**
     * Method that apply North-West Corner Method to current transportation problem
     *
     * @param Supply vector with supply values
     * @param Demand vector with demand values
     * @param Coefficients_Of_Costs matrix with cost values
     */
    public static void North_West_Corner_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs) {
        int n = Coefficients_Of_Costs.getNumRows(); // num of rows
        int m = Coefficients_Of_Costs.getNumCols(); // num of columns

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

        // Variables for saving result
        Matrix resultPath = new Matrix(n, m); // matrix of chosen costs
        double result = 0; // optimal path
        int i = 0; // current row
        int j = 0; // current column
        // Duplicates for output
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];

        // Haven't reached the end point yet
        while (i < n && j < m) {
            // Maximum available cost from supply and demand
            double currentValue = Math.min(supplyCopy.get(i), demandCopy.get(j));
            resultPath.setValue(i, j, currentValue);
            result += currentValue * costsCopy.getValue(i, j);
            // Decrease supply and demand by used cost
            supplyCopy.set(i, supplyCopy.get(i) - currentValue);
            demandCopy.set(j, demandCopy.get(j) - currentValue);
            // Move on the next cell
            if (supplyCopy.get(i) == 0) {
                i++;
            }
            if (demandCopy.get(j) == 0) {
                j++;
            }
        }
        // Preparation for output
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < m; l++) {
                // if we do not use this cell
                if (resultPath.getValue(k, l) == -1) {
                    resultPath.setValue(k, l, 0);
                }
                // Convert to string
                resultsCost[k][l] = DecimalFormat.getInstance().format(resultPath.getValue(k, l));
                resultsDemand[l] = DecimalFormat.getInstance().format(demandCopy.get(l));
            }
            resultsSupply[k] = DecimalFormat.getInstance().format(supplyCopy.get(k));
        }
        // Print it
        printTable(resultsSupply, resultsCost, resultsDemand);
        // Output optimal path
        System.out.println("Optimal path is " + DecimalFormat.getInstance().format(result));

    }

    /**
     * Method that apply Vogel's approximation Method to current transportation problem
     *
     * @param Supply vector with supply values
     * @param Demand vector with demand values
     * @param Coefficients_Of_Costs matrix with cost values
     */
    public static void Vogel_s_Approximation_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs){
        int n = Coefficients_Of_Costs.getNumRows(); // num of rows
        int m = Coefficients_Of_Costs.getNumCols(); // num of columns

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

        // Variables for saving result
        Matrix resultPath = new Matrix(n, m); // matrix of chosen costs
        double result = 0; // optimal path
        // Duplicates for output
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];

        // Haven't reached the end point yet
        while (true) {
            boolean notEnd = false; // flag to finish algorithm execution
            // Searching the difference between two minimal elements in rows and columns
            double[] differenceRows = new double[supplyCopy.size()];
            double[] differenceCols = new double[demandCopy.size()];
            // Differences in rows
            for (int i = 0; i < n; i++) {
                differenceRows[i] = costsCopy.sub_two_min_Elements_From_Row(i);
                // if there are values for iterations
                if (differenceRows[i] != -1) {
                    notEnd = true;
                }
            }
            // Differences in columns
            for (int i = 0; i < m; i++) {
                differenceCols[i] = costsCopy.sub_two_min_Elements_From_Column(i);
                // if there are values for iteration
                if (differenceCols[i] != -1) {
                    notEnd = true;
                }
            }
            // There is no any differences for iteration
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
            // Define current suitable line
            int curRow, curColumn;
            if (isCol) {
                curRow = miniIndexLine;
                curColumn = maxIndexInCol;
            } else {
                curRow = maxIndexInRow;
                curColumn = miniIndexLine;
            }
            // Maximum available cost from supply and demand
            double currentValue = Math.min(supplyCopy.get(curRow), demandCopy.get(curColumn));
            resultPath.setValue(curRow, curColumn, currentValue);
            result += currentValue * Coefficients_Of_Costs.getValue(curRow, curColumn);
            // Decrease supply and demand by used cost
            supplyCopy.set(curRow, supplyCopy.get(curRow) - currentValue);
            demandCopy.set(curColumn, demandCopy.get(curColumn) - currentValue);
            // Nullify cells that are not available because supply value = 0
            if (supplyCopy.get(curRow) == 0) {
                for (int i = 0; i < m; i++) {
                    costsCopy.setValue(curRow, i, Integer.MAX_VALUE);
                }
            }
            // Nullify cells that are not available because demand value = 0
            if (demandCopy.get(curColumn) == 0) {
                for (int i = 0; i < n; i++) {
                    costsCopy.setValue(i, curColumn, Integer.MAX_VALUE);
                }
            }
        }
        // Preparation for output
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < m; l++) {
                // if we do not use this cell
                if (resultPath.getValue(k, l) == -1) {
                    resultPath.setValue(k, l, 0);
                }
                // Convert to string
                resultsCost[k][l] = DecimalFormat.getInstance().format(resultPath.getValue(k, l));
                resultsDemand[l] = DecimalFormat.getInstance().format(demandCopy.get(l));
            }
            resultsSupply[k] = DecimalFormat.getInstance().format(supplyCopy.get(k));
        }
        // Print it
        printTable(resultsSupply, resultsCost, resultsDemand);
        // Output optimal path
        System.out.println("Optimal path is " + DecimalFormat.getInstance().format(result));
    }

    /**
     * Method that apply Russell's approximation Method to current transportation problem
     *
     * @param Supply vector with supply values
     * @param Demand vector with demand values
     * @param Coefficients_Of_Costs matrix with cost values
     */
    public static void Russell_s_Approximation_Method(Vector Supply, Vector Demand, Matrix Coefficients_Of_Costs){
        // Copies the supply and demand vectors to avoid modifying the original values during calculation
        Vector supplyCopy = new Vector(Supply.size());
        supplyCopy.copy(Supply);
        Vector demandCopy = new Vector(Demand.size());
        demandCopy.copy(Demand);

        // Get dimensions of the cost matrix
        int n = Coefficients_Of_Costs.getNumRows();
        int m = Coefficients_Of_Costs.getNumCols();

        // Initialize matrices to store the result path and an intermediate result path with accumulated costs
        Matrix resultPath = new Matrix(n, m);
        Matrix resultPathCopy = new Matrix(n, m);

        // Initialize arrays to hold the formatted results for output
        String[][] resultsCost = new String[n][m];
        String[] resultsSupply = new String[n];
        String[] resultsDemand = new String[m];
        double result = 0;// Variable to store the total transportation cost

        // Create a copy of the cost matrix for manipulation during calculations
        Matrix costsCopy = new Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                costsCopy.setValue(i, j, Coefficients_Of_Costs.getValue(i, j));
            }
        }

        // Main loop to determine optimal paths until demand is fulfilled
        while (true){
            int x = 0, y = 0;
            double max = Double.MIN_VALUE;
            for (int i = 0; i < n; i++) {
                double row_max = costsCopy.max_Elements_From_Row(i, costsCopy);
                for (int j = 0; j < m; j++) {
                    double col_max = costsCopy.max_Elements_From_Column(j, costsCopy);

                    // Calculate the opportunity cost difference and update x, y if max difference found
                    if (costsCopy.getValue(i,j) != -1 && max < Math.abs(costsCopy.getValue(i, j) - col_max - row_max)) {
                        x = i;
                        y = j;
                        max = Math.abs(costsCopy.getValue(i, j) - col_max - row_max);

                    }
                }
            }

            // Break loop if the demand has been fully satisfied
            if (demandCopy.getSumValues() == 0){
                break;
            }

            // Initialize resultPath if not set
            if (resultPath.getValue(x, y) == -1) {
                resultPath.setValue(x, y, 0);
            }

            // Allocate supply and demand values based on the selected cell (x, y)
            if ((demandCopy.get(y) - supplyCopy.get(x)) == 0) {
                // Case where supply exactly meets demand
                resultPathCopy.setValue(x, y, resultPath.getValue(x, y) + costsCopy.getValue(x, y) * demandCopy.get(y));
                resultPath.setValue(x, y, Math.min(demandCopy.get(y),supplyCopy.get(x)));
                supplyCopy.set(x, 0); // Set supply to zero as it's used up
                demandCopy.set(y, 0); // Set demand to zero as it's fulfilled
                for (int i = 0; i < n; i++) {
                    costsCopy.setValue(i, y, -1); // Set column values to -1 to prevent further allocation to this demand
                }
            } else if ((demandCopy.get(y) - supplyCopy.get(x)) > 0) {
                // Case where demand is greater than supply
                resultPathCopy.setValue(x, y, resultPath.getValue(x, y) + costsCopy.getValue(x, y) * supplyCopy.get(x));
                resultPath.setValue(x, y, Math.min(demandCopy.get(y),supplyCopy.get(x)));
                demandCopy.set(y, demandCopy.get(y) - supplyCopy.get(x)); // Reduce demand by the supply used
                supplyCopy.set(x, 0); // Supply is depleted for this source
                costsCopy.setValue(x, y, -1); // Set the cost to -1 to prevent re-allocation to this cell
            } else if ((demandCopy.get(y) - supplyCopy.get(x)) < 0) {
                // Case where supply is greater than demand
                resultPathCopy.setValue(x, y, resultPath.getValue(x, y) + costsCopy.getValue(x, y) * demandCopy.get(y));
                resultPath.setValue(x, y, Math.min(demandCopy.get(y),supplyCopy.get(x)));
                supplyCopy.set(x, supplyCopy.get(x) - demandCopy.get(y)); // Reduce supply by the demand used
                demandCopy.set(y, 0); // Demand is fulfilled
                for (int i = 0; i < n; i++) {
                    costsCopy.setValue(i, y, -1); // Set column to -1 to prevent further allocation
                }
            }
            costsCopy.setValue(x, y, -1); // Reset cost for the processed cell
        }

        // Final calculations and formatting for output
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < m; l++) {
                if (resultPathCopy.getValue(k, l) == -1) {
                    resultPathCopy.setValue(k, l, 0); // Set unassigned cells to zero in result path copy
                } else {
                    result += resultPathCopy.getValue(k, l); // Sum up the total transportation cost
                }
                if (resultPath.getValue(k, l) == -1) {
                    resultPath.setValue(k, l, 0); // Set unassigned cells to zero in result path
                }
                resultsCost[k][l] = DecimalFormat.getInstance().format(resultPath.getValue(k, l)); // Format result for each cell
                resultsDemand[l] = DecimalFormat.getInstance().format(demandCopy.get(l)); // Format remaining demand for each destination
            }
            resultsSupply[k] = DecimalFormat.getInstance().format(supplyCopy.get(k)); // Format remaining supply for each source
        }

        // Print the final results in a table format
        printTable(resultsSupply, resultsCost, resultsDemand);
        System.out.println("Optimal path is " + DecimalFormat.getInstance().format(result)); // Print total cost
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
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.values[i][j] = -1;
            }
        }
    }

    /**
     * Constructor for creating a matrix with specified dimensions and values.
     *
     * @param values 2D array containing matrix elements
     */
    public Matrix(double[][] values)  {
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
     * Returns the row at this position
     *
     * @param i row index
     * @return the row from the table at the specified position
     */
    public double[] getRow(int i) {
        return values[i];
    }

    /**
     * Returns the column at this position
     *
     * @param i column index
     * @return the column from the table at the specified position
     */
    public double[] getColumn(int i) {
        double[] column = new double[rows];
        for (int row = 0; row < rows; row++) {
            column[row] = values[row][i];
        }
        return column;
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
     * Returns the difference between two minimal values from current row.
     *
     * @param row index of the row that are considered
     * @return the difference between two minimal values
     */
    public double sub_two_min_Elements_From_Row(int row) {
        double[] Row = getRow(row);
        // Copy values not to change them in the table
        double[] rowLine = new double[Row.length];
        for (int i = 0; i < Row.length; i++) {
            rowLine[i] = Row[i];
        }
        // Sort values in increasing order
        Arrays.sort(rowLine);
        // if there is no any available values on the table
        if (rowLine[0] == Integer.MAX_VALUE && rowLine[1] == Integer.MAX_VALUE) {
            return -1;
        }
        // if there is only one available value
        if (rowLine[1] == Integer.MAX_VALUE) {
            return rowLine[0];
        }
        if (rowLine[0] == Integer.MAX_VALUE) {
            return rowLine[1];
        }
        return Math.abs(rowLine[1]-rowLine[0]);
    }

    /**
     * Returns the maximum element from current row.
     *
     * @param row index of the row that are considered
     * @param matrix matrix that are considered
     * @return the maximum element from row
     */
    public double max_Elements_From_Row(int row, Matrix matrix) {
        int n = getNumCols();
        double[] Row = getRow(row);
        // Copy values not to change them in the table
        double[] newRow = new double[Row.length];
        for (int i = 0; i < Row.length; i++) {
            newRow[i] = Row[i];
        }
        // Sort values in increasing order
        Arrays.sort(newRow);
        return newRow[n-1];
    }

    /**
     * Returns the difference between two minimal values from current column.
     *
     * @param column index of the column that are considered
     * @return the difference between two minimal values
     */
    public double sub_two_min_Elements_From_Column(int column) {
        double[] Column = getColumn(column);
        // Copy values not to change them in the table
        double[] columnLine = new double[Column.length];
        for (int i = 0; i < Column.length; i++) {
            columnLine[i] = Column[i];
        }
        // Sort values in increasing order
        Arrays.sort(columnLine);
        // if there is no any available values on the table
        if (columnLine[0] == Integer.MAX_VALUE && columnLine[1] == Integer.MAX_VALUE) {
            return -1;
        }
        // if there is only one available value
        if (columnLine[1] == Integer.MAX_VALUE) {
            return columnLine[0];
        }
        if (columnLine[0] == Integer.MAX_VALUE) {
            return columnLine[1];
        }
        return Math.abs(columnLine[1]-columnLine[0]);
    }

    /**
     * Returns the maximum element from current column.
     *
     * @param column index of the column that are considered
     * @param matrix matrix that are considered
     * @return the maximum element from column
     */
    public double max_Elements_From_Column(int column, Matrix matrix) {
        int n = getNumRows();
        double[] Column = getColumn(column);
        // Copy values not to change them in the table
        double[] newColumn = new double[Column.length];
        for (int i = 0; i < Column.length; i++) {
            newColumn[i] = Column[i];
        }
        // Sort values in increasing order
        Arrays.sort(newColumn);
        return newColumn[n-1];
    }
}

/**
 * Vector class representing a mathematical vector.
 */
class Vector {
    private double[] values;

    /**
     * Constructor for creating a vector without values.
     *
     * @param n size of vector
     */
    public Vector(int n) { this.values = new double[n]; }

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
     * Set a new value at a specific index.
     *
     * @param index index of the element to be changed
     * @param value new value at this position
     */
    public void set(int index, double value) { values[index] = value; }

    /**
     * Returns the size of the vector.
     *
     * @return the size of the vector
     */
    public int size() {
        return values.length;
    }

    /**
     * Make a copy of vector.
     *
     * @param vector the vector that will be copied
     */
    public void copy(Vector vector) {
        for (int i = 0; i < vector.size(); i++) {
            values[i] = vector.get(i);
        }
    }

    /**
     * Return sum of all values of current vector
     *
     * @return sum of all values
     */
    public double getSumValues() {
        double sum = 0.0;
        for (double v : values) {
            sum += v;
        }
        return sum;
    }
}
