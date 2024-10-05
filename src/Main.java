import java.sql.SQLOutput;
import java.util.Scanner;
import java.math.BigDecimal;
import java.math.RoundingMode;

/*
 * Формат ввода:
 * 0. min или max: что мы хотим сделать
 * 1. Первая строка: количество переменных целевой функции (functionSize).
 * 2. Вторая строка: коэффициенты целевой функции (coeffFunction), перечисленные через пробел.
 * 3. Третья строка: количество ограничений (constraintSize).
 * 4. Следующие строки (ConstraintSize штук): коэффициенты каждого ограничения
 *    (coeffConstraints), по одному ограничению в строке.
 * 5. Далее строка с правыми частями ограничений (rightHandValues), перечисленными через пробел.
 * 6. Последняя строка: точность (accuracy)
 *
 * Пример ввода:
 * max
 * 3                      // количество переменных
 * 9 10 16                // коэффициенты целевой функции F(x1, x2, x3) = 9x1 + 10x2 + 16x3
 * 3                      // количество ограничений
 * 18 15 12               // первое ограничение: 18x1 + 15x2 + 12x3 <= 360
 * 6 4 8                  // второе ограничение: 6x1 + 4x2 + 8x3 <= 192
 * 5 3 3                  // третье ограничение: 5x1 + 3x2 + 3x3 <= 180
 * 360 192 180            // правые части ограничений
 * 6                      // точность
 *
 * Формат вывода:
 * 1. Выводится оптимальное решение в виде значений переменных x1, x2, ..., xn.
 * 2. Каждая переменная выводится в формате "xN = значение", где N — номер переменной.
 *
 * Пример вывода:
 * Optimal solution:
 * x1 = 0.00
 * x2 = 8.00
 * x3 = 20.00
 */


/**
 * Main class that handles user input and applies the Simplex method for solving
 * linear programming problems.
 */
public class Main {

    /**
     * The main method to read user input and start the Simplex method.
     *
     * @param args command-line arguments (not used)
     */
    public static void main(String[] args) {

        // Input
        Scanner scanner = new Scanner(System.in);
        String action = scanner.nextLine();
        int functionSize = scanner.nextInt();
        int[] coeffFunction = new int[functionSize];
        for (int i = 0; i < functionSize; i++) {
            coeffFunction[i] = scanner.nextInt();
        }
        int constraintSize = scanner.nextInt();
        int[][] coeffConstraints = new int[constraintSize][functionSize];
        for (int i = 0; i < coeffConstraints.length; i++) {
            for (int j = 0; j < coeffConstraints[0].length; j++) {
                coeffConstraints[i][j] = scanner.nextInt();
            }
        }
        int[] rightHandValues = new int[constraintSize];
        for (int i = 0; i < rightHandValues.length; i++) {
            rightHandValues[i] = scanner.nextInt();
        }
        int accuracy = scanner.nextInt();

        // Creating elements
        Vector function = new Vector(functionSize, coeffFunction);
        Matrix constraints = new Matrix(constraintSize, functionSize, coeffConstraints);
        Vector rightConstraints = new Vector(functionSize, rightHandValues);

        // Print the initial problem
        System.out.print(action + " z = ");
        for (int i = 0; i < function.getSize(); i++) {
            System.out.printf("%s%d•x%d",
                    i == 0 ? (function.getValue(i) < 0 ? "-" : "") : (function.getValue(i) < 0 ? " - " : " + "),
                    Math.abs(function.getValue(i)),
                    i + 1);
        }
        System.out.println();//make new line
        System.out.println("Subject to the constrains:");
        for (int i = 0; i < constraints.getNumRows(); i++) {
            for (int j = 0; j < constraints.getNumCols(); j++) {
                System.out.printf("%s%d•x%d",
                        j == 0 ? (constraints.getValue(i, j) < 0 ? "-" : "") :
                                (constraints.getValue(i, j) < 0 ? " - " : " + "),
                        Math.abs(constraints.getValue(i, j)), j + 1);
            }
            System.out.printf("≤ %d\n", rightConstraints.getValue(i));
        }


        // Applying simplex method
        simplexMethod(function, constraints, rightConstraints, accuracy, action);
    }

    /**
     * Executes the Simplex method for solving linear programming problems.
     *
     * @param C        the vector representing the objective function's coefficients
     * @param A        the matrix representing the constraints' coefficients
     * @param b        the vector representing the right-hand side of the constraints
     * @param accuracy the accuracy to be applied during calculations
     * @param action   specifies whether the problem is a maximization ("max") or minimization ("min")
     */
    public static void simplexMethod(Vector C, Matrix A, Vector b, int accuracy, String action) {
        int numRows = A.getNumRows();
        int numCols = A.getNumCols();

        // Tableau creation
        double[][] tableau = new double[numRows + 1][numCols + numRows + 1];

        // Fill tableau with constraints
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                tableau[i][j] = A.getValue(i, j);
            }
            tableau[i][numCols + i] = 1; // Slack variable
            tableau[i][numCols + numRows] = b.getValue(i);
        }

        // Fill last row with objective function
        for (int j = 0; j < numCols; j++) {
            tableau[numRows][j] = -C.getValue(j);
        }

        // Perform the simplex method iteration
        while (true) {
            //System.out.println("Current tableau:");
            //printTableau(tableau, numRows, numCols + numRows, accuracy);

            int pivotColumn = findPivotColumn(tableau, numRows, numCols + numRows, action);
            if (pivotColumn == -1) {
                System.out.println("Optimal solution found.");
                break; // Solution is optimal
            }

            int pivotRow = findPivotRow(tableau, pivotColumn, numRows);
            if (pivotRow == -1) {
                System.out.println("The solution is unbounded.");
                return;
            }

            // Perform pivot
            //System.out.println("Pivoting on row " + pivotRow + " and column " + pivotColumn);
            pivot(tableau, pivotRow, pivotColumn, numRows, numCols + numRows, accuracy);

            // Print tableau after each iteration
            //System.out.println("Tableau after pivoting:");
            //printTableau(tableau, numRows, numCols + numRows, accuracy);
        }

        // Output the solution
        printSolution(tableau, numRows, numCols, C, accuracy);
    }

    /**
     * Finds the pivot column based on the action ("max" or "min").
     *
     * @param tableau the tableau to be evaluated
     * @param numRows the number of rows in the tableau
     * @param numCols the number of columns in the tableau
     * @param action  the action specifying if the problem is a maximization ("max") or minimization ("min")
     * @return the index of the pivot column, or -1 if no valid pivot column is found
     */
    public static int findPivotColumn(double[][] tableau, int numRows, int numCols, String action) {
        int pivotColumn = -1;
        if (action.equals("max")) {
            double minValue = 0;
            for (int j = 0; j < numCols; j++) {
                if (tableau[numRows][j] < minValue) {
                    minValue = tableau[numRows][j];
                    pivotColumn = j;
                }
            }
        } else if (action.equals("min")) {
            double maxValue = 0;
            for (int j = 0; j < numCols; j++) {
                if (tableau[numRows][j] > maxValue) {
                    maxValue = tableau[numRows][j];
                    pivotColumn = j;
                }
            }
        }
        return pivotColumn;
    }

    /**
     * Finds the pivot row for the given pivot column.
     *
     * @param tableau     the tableau to be evaluated
     * @param pivotColumn the index of the pivot column
     * @param numRows     the number of rows in the tableau
     * @return the index of the pivot row, or -1 if the solution is unbounded
     */
    public static int findPivotRow(double[][] tableau, int pivotColumn, int numRows) {
        int pivotRow = -1;
        double minRatio = Double.MAX_VALUE;
        for (int i = 0; i < numRows; i++) {
            if (tableau[i][pivotColumn] > 0) {
                double ratio = tableau[i][tableau[0].length - 1] / tableau[i][pivotColumn];
                if (ratio < minRatio) {
                    minRatio = ratio;
                    pivotRow = i;
                }
            }
        }
        return pivotRow;
    }

    /**
     * Performs the pivot operation on the tableau.
     *
     * @param tableau     the tableau to be updated
     * @param pivotRow    the index of the pivot row
     * @param pivotColumn the index of the pivot column
     * @param numRows     the number of rows in the tableau
     * @param numCols     the number of columns in the tableau
     * @param accuracy    the decimal precision to apply to calculations
     */
    public static void pivot(double[][] tableau, int pivotRow, int pivotColumn, int numRows, int numCols,
                             int accuracy) {
        double pivotValue = tableau[pivotRow][pivotColumn];

        // Normalize the pivot row
        for (int j = 0; j < numCols + 1; j++) {
            tableau[pivotRow][j] = round(tableau[pivotRow][j] / pivotValue, accuracy);
        }

        // Make other rows zero in pivot column
        for (int i = 0; i <= numRows; i++) {
            if (i != pivotRow) {
                double factor = tableau[i][pivotColumn];
                for (int j = 0; j < numCols + 1; j++) {
                    tableau[i][j] = round(tableau[i][j] - factor * tableau[pivotRow][j], accuracy);
                }
            }
        }
    }

    /**
     * Prints the current state of the tableau.
     *
     * @param tableau  the tableau to be printed
     * @param numRows  the number of rows in the tableau
     * @param numCols  the number of columns in the tableau
     * @param accuracy the decimal precision to apply when printing
     */
    public static void printTableau(double[][] tableau, int numRows, int numCols, int accuracy) {
        for (int i = 0; i <= numRows; i++) {
            for (int j = 0; j <= numCols; j++) {
                System.out.printf("%." + accuracy + "f ", tableau[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    /**
     * Prints the optimal solution after the Simplex method is completed.
     *
     * @param tableau  the final tableau containing the solution
     * @param numRows  the number of rows in the tableau
     * @param numCols  the number of columns representing variables
     * @param C        the vector representing the objective function's coefficients
     * @param accuracy the decimal precision to apply when printing
     */
    public static void printSolution(double[][] tableau, int numRows, int numCols, Vector C, int accuracy) {
        double[] solution = new double[numCols];
        for (int i = 0; i < numCols; i++) {
            boolean isBasic = true;
            int rowIndex = -1;
            for (int j = 0; j < numRows; j++) {
                if (tableau[j][i] == 1) {
                    if (rowIndex == -1) {
                        rowIndex = j;
                    } else {
                        isBasic = false;
                        break;
                    }
                } else if (tableau[j][i] != 0) {
                    isBasic = false;
                    break;
                }
            }
            if (isBasic && rowIndex != -1) {
                solution[i] = tableau[rowIndex][tableau[0].length - 1];
            } else {
                solution[i] = 0;
            }
        }

        // Output solution
        System.out.println("Optimal solution:");
        for (int i = 0; i < solution.length; i++) {
            System.out.printf("x%d = %." + accuracy + "f\n", i + 1, solution[i]);
        }

        // Calculate optimal value
        double optimalValue = 0;
        for (int i = 0; i < numCols; i++) {
            optimalValue += C.getValue(i) * solution[i];
        }
        System.out.printf("Optimal value: %." + accuracy + "f\n", optimalValue);
    }

    /**
     * Rounds a given value to the specified number of decimal places.
     *
     * @param value    the value to be rounded
     * @param accuracy the number of decimal places
     * @return the rounded value
     */
    public static double round(double value, int accuracy) {
        BigDecimal bd = new BigDecimal(Double.toString(value));
        bd = bd.setScale(accuracy, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }
}

/**
 * Vector class representing a mathematical vector.
 */
class Vector {
    private int n; // Size of the vector
    private int[] values; // Array storing vector elements

    /**
     * Constructor for creating a vector with specified size and values.
     *
     * @param n      size of the vector
     * @param values array containing vector elements
     */
    Vector(int n, int[] values) {
        this.n = n;
        this.values = values;
    }

    /**
     * Returns the value of the vector at a specific index.
     *
     * @param i index of the element to be returned
     * @return the value at the specified index
     */
    public int getValue(int i) {
        return values[i];
    }

    /**
     * Returns the size of the vector.
     *
     * @return the size of the vector
     */
    public int getSize() {
        return n;
    }
}

/**
 * Matrix class representing a mathematical matrix.
 */
class Matrix {
    private int n; // Number of rows
    private int m; // Number of columns
    private int[][] values; // 2D array storing matrix elements

    /**
     * Constructor for creating a matrix with specified dimensions and values.
     *
     * @param n      number of rows
     * @param m      number of columns
     * @param values 2D array containing matrix elements
     */
    Matrix(int n, int m, int[][] values) {
        this.n = n;
        this.m = m;
        this.values = values;
    }

    /**
     * Returns the value of the matrix at a specific position.
     *
     * @param i row index
     * @param j column index
     * @return the value at the specified position
     */
    public int getValue(int i, int j) {
        return values[i][j];
    }

    /**
     * Returns the number of rows in the matrix.
     *
     * @return the number of rows
     */
    public int getNumRows() {
        return n;
    }

    /**
     * Returns the number of columns in the matrix.
     *
     * @return the number of columns
     */
    public int getNumCols() {
        return m;
    }
}


