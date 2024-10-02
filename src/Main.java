import java.util.Scanner;

/*
 * Формат ввода:
 * 0. min или max: что мы хотим сделать
 * 1. Первая строка: количество переменных целевой функции (functionSize).
 * 2. Вторая строка: коэффициенты целевой функции (coeffFunction), перечисленные через пробел.
 * 3. Третья строка: количество ограничений (rowConstraintSize).
 * 4. Четвертая строка: количество переменных в каждом ограничении (columnConstraintSize).
 * 5. Следующие строки (rowConstraintSize штук): коэффициенты каждого ограничения
 *    (coeffConstraints), по одному ограничению в строке.
 * 6. Строка после ограничений: количество правых частей ограничений (rightHandSize).
 * 7. Далее строка с правыми частями ограничений (rightHandValues), перечисленными через пробел.
 * 8. Последняя строка: точность (accuracy) — на данном этапе точность может быть произвольной и игнорироваться.
 *
 * Пример ввода:
 * 3                      // количество переменных
 * 9 10 16                // коэффициенты целевой функции F(x1, x2, x3) = 9x1 + 10x2 + 16x3
 * 3                      // количество ограничений
 * 3                      // количество переменных в каждом ограничении
 * 18 15 12               // первое ограничение: 18x1 + 15x2 + 12x3 <= 360
 * 6 4 8                  // второе ограничение: 6x1 + 4x2 + 8x3 <= 192
 * 5 3 3                  // третье ограничение: 5x1 + 3x2 + 3x3 <= 180
 * 3                      // количество правых частей ограничений
 * 360 192 180            // правые части ограничений
 * 6                      // точность (можно игнорировать)
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

        // input
        Scanner scanner = new Scanner(System.in);
        String action = scanner.nextLine();
        int functionSize = scanner.nextInt();
        int[] coeffFunction = new int[functionSize];
        for (int i = 0; i < functionSize; i++) {
            coeffFunction[i] = scanner.nextInt();
        }
        int rowConstraintSize = scanner.nextInt();
        int columnConstraintSize = scanner.nextInt();
        int[][] coeffConstraints = new int[rowConstraintSize][columnConstraintSize];
        for (int i = 0; i < rowConstraintSize; i++) {
            for (int j = 0; j < columnConstraintSize; j++) {
                coeffConstraints[i][j] = scanner.nextInt();
            }
        }
        int rightHandSize = scanner.nextInt();
        int[] rightHandValues = new int[rightHandSize];
        for (int i = 0; i < rightHandSize; i++) {
            rightHandValues[i] = scanner.nextInt();
        }
        int accuracy = scanner.nextInt();

        // creating elements
        Vector function = new Vector(functionSize, coeffFunction);
        Matrix constraints = new Matrix(rowConstraintSize, columnConstraintSize, coeffConstraints);
        Vector rightConstraints = new Vector(rightHandSize, rightHandValues);

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
            System.out.println("Current tableau:");
            printTableau(tableau, numRows, numCols + numRows);

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
            System.out.println("Pivoting on row " + pivotRow + " and column " + pivotColumn);
            pivot(tableau, pivotRow, pivotColumn, numRows, numCols + numRows);

            // Print tableau after each iteration
            System.out.println("Tableau after pivoting:");
            printTableau(tableau, numRows, numCols + numRows);
        }

        // Output the solution
        printSolution(tableau, numRows, numCols);
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
     */
    public static void pivot(double[][] tableau, int pivotRow, int pivotColumn, int numRows, int numCols) {
        double pivotValue = tableau[pivotRow][pivotColumn];

        // Normalize the pivot row
        for (int j = 0; j < numCols + 1; j++) {
            tableau[pivotRow][j] /= pivotValue;
        }

        // Make other rows zero in pivot column
        for (int i = 0; i <= numRows; i++) {
            if (i != pivotRow) {
                double factor = tableau[i][pivotColumn];
                for (int j = 0; j < numCols + 1; j++) {
                    tableau[i][j] -= factor * tableau[pivotRow][j];
                }
            }
        }
    }

    /**
     * Prints the current state of the tableau.
     *
     * @param tableau the tableau to be printed
     * @param numRows the number of rows in the tableau
     * @param numCols the number of columns in the tableau
     */
    public static void printTableau(double[][] tableau, int numRows, int numCols) {
        for (int i = 0; i <= numRows; i++) {
            for (int j = 0; j <= numCols; j++) {
                System.out.printf("%.2f ", tableau[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    /**
     * Prints the optimal solution after the Simplex method is completed.
     *
     * @param tableau the final tableau containing the solution
     * @param numRows the number of rows in the tableau
     * @param numCols the number of columns representing variables
     */
    public static void printSolution(double[][] tableau, int numRows, int numCols) {
        double[] solution = new double[numCols]; // numCols refers to the number of variables (without slack)

        // Find basic variables, among which only x1, x2, x3 are of interest
        for (int j = 0; j < numCols; j++) {
            boolean isBasic = false;
            double value = 0.0;
            for (int i = 0; i < numRows; i++) {
                if (tableau[i][j] == 1) {
                    isBasic = true;
                    value = tableau[i][tableau[0].length - 1];
                } else if (tableau[i][j] != 0) {
                    isBasic = false;
                    break;
                }
            }
            if (isBasic) {
                solution[j] = value;
            } else {
                solution[j] = 0.0; // If the variable is non-basic, it is set to 0
            }
        }

        // Output the optimal solution
        System.out.println("Optimal solution:");
        for (int j = 0; j < numCols; j++) {
            System.out.printf("x%d = %.2f\n", j + 1, solution[j]);
        }
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


