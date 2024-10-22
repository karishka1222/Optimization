import java.util.Arrays;

// Data for test from test 4 for Simplex Method
// нам надо унифицировать ввод для 2-х методов и поправить структуру кода

public class InteriorPointMethod {

    public static void main(String[] args) {
        // Coefficients for the objective function: Maximize 9x1 + 10x2 + 16x3
        double[] c = {9, 10, 16, 0, 0, 0}; // Slack variables have zero coefficients Test 4

//        double[] c = {3, 1, 3, 0, 0, 0}; // Slack variables have zero coefficients Test 4

        // Coefficients for the constraints (including slack variables) Test 4
        double[][] A = {
                {18, 15, 12, 1, 0, 0}, // First constraint
                {6, 4, 8, 0, 1, 0},    // Second constraint
                {5, 3, 3, 0, 0, 1}     // Third constraint
        };

//        double[][] A = {
//                {2, 1, 1, 1, 0, 0}, // First constraint
//                {1, 2, 3, 0, 1, 0},    // Second constraint
//                {2, 2, 1, 0, 0, 1}     // Third constraint
//        };

        // Right-hand side of the constraints Test 4
        double[] b = {360, 192, 180};

//        double[] b = {2, 5, 6};

        // Initial solution: (x1, x2, x3, x4, x5, x6) = (1, 1, 1, 315, 174, 169)
        double[] x = {1, 1, 1, 315, 174, 169};

        double alpha = 0.5;
        int iteration = 1;

        System.out.println("For alpha = " + alpha);

        while (true) {
            double[] v = Arrays.copyOf(x, x.length); // Save a copy of x
            double[][] D = diag(x);                  // Create diagonal matrix from x
            double[][] AA = dot(A, D);               // A * D
            double[] cc = dot(D, c);                 // D * c
            double[][] I = eye(c.length);            // Identity matrix
            double[][] F = dot(AA, transpose(AA));   // AA * AA^T
            double[][] FI = inverse(F);              // Inverse of F
            double[][] H = dot(transpose(AA), FI);   // AA^T * FI
            double[][] P = subtract(I, dot(H, AA));  // P = I - H * AA
            double[] cp = dot(P, cc);                // cp = P * cc

            double nu = Math.abs(min(cp));           // Calculate nu
            if (nu == 0) {
                System.out.println("Warning: nu is zero, cp: " + Arrays.toString(cp));
                break; // Handle zero nu case
            }

//             y = 1 + (alpha / nu) * cp
            double[] y = add(ones(cp.length), multiply(cp, alpha / nu));
            double[] yy = dot(D, y);                 // yy = D * y
            x = yy;                                  // Update x

            // Display iteration steps
            if (iteration <= 4) {
                System.out.println("In iteration " + iteration + " we have x = " + Arrays.toString(x) + "\n");
            }
            iteration++;

//             Termination condition
            if (norm(subtract(yy, v), 2) < 0.0001) {
                break;
            }
        }

        //results
        System.out.println("In the last iteration we have x = " + Arrays.toString(x) + "\n");
        System.out.println(9 * x[0] + 10 * x[1] + 16 * x[2]);
    }

    // желательно оформить в методы класса Matrix и/или Vector(но это если по-хорошему сделать)

    // D matrix
    public static double[][] diag(double[] vector) {
        double[][] result = new double[vector.length][vector.length];
        for (int i = 0; i < vector.length; i++) {
            result[i][i] = vector[i];
        }
        return result;
    }

    // MUX for matrix
    public static double[][] dot(double[][] A, double[][] B) {
        int rowsA = A.length;           // number of row in A
        int colsA = A[0].length;       // number of columns in A
        int rowsB = B.length;           // number of row in B
        int colsB = B[0].length;       // number of columns in B

        // check of validity of matrix size
        if (colsA != rowsB) {
            System.out.println("Invalid multiplication: colsA (" + colsA + ") != rowsB (" + rowsB + ")");
            throw new IllegalArgumentException("Number of columns in A must be equal to the number of rows in B.");
        }

        double[][] result = new double[rowsA][colsB];

        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsB; j++) {
                for (int k = 0; k < colsA; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }


    // MUX for vector
    public static double[] dot(double[][] A, double[] B) {
        int rowsA = A.length;
        int colsA = A[0].length;
        double[] result = new double[rowsA];
        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsA; j++) {
                result[i] += A[i][j] * B[j];
            }
        }
        return result;
    }

    // T for matrix
    public static double[][] transpose(double[][] A) {
        int rowsA = A.length;
        int colsA = A[0].length;
        double[][] result = new double[colsA][rowsA];
        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsA; j++) {
                result[j][i] = A[i][j];
            }
        }
        return result;
    }

    // Inverse
    public static double[][] inverse(double[][] A) {
        int n = A.length;
        // Проверка, что матрица является квадратной
        if (n == 0 || A[0].length != n) {
            throw new IllegalArgumentException("Матрица должна быть квадратной.");
        }

        // Matrix [A | I]
        double[][] augmentedMatrix = new double[n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = A[i][j];
            }
            augmentedMatrix[i][i + n] = 1; // Заполнение единичной матрицы
        }

        for (int i = 0; i < n; i++) {
            // Max item in column
            double maxEl = Math.abs(augmentedMatrix[i][i]);
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(augmentedMatrix[k][i]) > maxEl) {
                    maxEl = Math.abs(augmentedMatrix[k][i]);
                    maxRow = k;
                }
            }

            double[] temp = augmentedMatrix[maxRow];
            augmentedMatrix[maxRow] = augmentedMatrix[i];
            augmentedMatrix[i] = temp;

            // Upper triangle form
            for (int k = i + 1; k < n; k++) {
                double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
                for (int j = 0; j < 2 * n; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }

        for (int i = n - 1; i >= 0; i--) {

            double leadingElement = augmentedMatrix[i][i];
            for (int j = 0; j < 2 * n; j++) {
                augmentedMatrix[i][j] /= leadingElement;
            }
            // pivoting
            for (int k = i - 1; k >= 0; k--) {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }

        double[][] inverseMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inverseMatrix[i][j] = augmentedMatrix[i][j + n];
            }
        }

        return inverseMatrix;
    }

    // Identity matrix
    public static double[][] eye(int size) {
        double[][] result = new double[size][size];
        for (int i = 0; i < size; i++) {
            result[i][i] = 1;
        }
        return result;
    }

    // subtraction for matrix
    public static double[][] subtract(double[][] A, double[][] B) {
        int rows = A.length;
        int cols = A[0].length;
        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = A[i][j] - B[i][j];
            }
        }
        return result;
    }

    // Identity vector
    public static double[] ones(int size) {
        double[] result = new double[size];
        Arrays.fill(result, 1);
        return result;
    }

    // MUX for vector
    public static double[] multiply(double[] vector, double scalar) {
        double[] result = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            result[i] = vector[i] * scalar;
        }
        return result;
    }

    // instruction for vector
    public static double[] subtract(double[] A, double[] B) {
        double[] result = new double[A.length];
        for (int i = 0; i < A.length; i++) {
            result[i] = A[i] - B[i];
        }
        return result;
    }

    // norm of vector
    public static double norm(double[] vector, int ord) {
        double sum = 0;
        for (double v : vector) {
            sum += Math.pow(v, ord);
        }
        return Math.pow(sum, 1.0 / ord);
    }

    // min item of vector
    public static double min(double[] vector) {
        double minValue = vector[0];
        for (double v : vector) {
            if (v < minValue) {
                minValue = v;
            }
        }
        return minValue;
    }

    // adding for vector
    public static double[] add(double[] A, double[] B) {
        double[] result = new double[A.length];
        for (int i = 0; i < A.length; i++) {
            result[i] = A[i] + B[i];
        }
        return result;
    }
}
















// по-идее на всё это ООП можно забить конечно

//class Vector {
//    private int n;          // Размер вектора
//    private double[] values; // Массив для хранения элементов вектора
//
//    /**
//     * Конструктор для создания вектора заданного размера.
//     * Все элементы вектора инициализируются нулем.
//     *
//     * @param n размер вектора
//     */
//    public Vector(int n) {
//        this.n = n;
//        this.values = new double[n];
//    }
//
//    /**
//     * Конструктор для создания вектора с заданными значениями.
//     *
//     * @param values массив значений вектора
//     */
//    public Vector(double[] values) {
//        this.n = values.length;
//        this.values = Arrays.copyOf(values, n);
//    }
//
//    /**
//     * Получает значение вектора по заданному индексу.
//     *
//     * @param i индекс элемента
//     * @return значение элемента по указанному индексу
//     */
//    public double getValue(int i) {
//        if (i < 0 || i >= n) {
//            throw new IndexOutOfBoundsException("Индекс вне диапазона.");
//        }
//        return values[i];
//    }
//
//    /**
//     * Устанавливает значение элемента по указанному индексу.
//     *
//     * @param i индекс элемента
//     * @param value новое значение элемента
//     */
//    public void setValue(int i, double value) {
//        if (i < 0 || i >= n) {
//            throw new IndexOutOfBoundsException("Индекс вне диапазона.");
//        }
//        values[i] = value;
//    }
//
//    /**
//     * Возвращает размер вектора.
//     *
//     * @return размер вектора
//     */
//    public int getSize() {
//        return n;
//    }
//
//    /**
//     * Складывает два вектора.
//     *
//     * @param other другой вектор
//     * @return результат сложения в виде нового вектора
//     */
//    public Vector add(Vector other) {
//        if (n != other.getSize()) {
//            throw new IllegalArgumentException("Размеры векторов должны совпадать.");
//        }
//        double[] result = new double[n];
//        for (int i = 0; i < n; i++) {
//            result[i] = this.values[i] + other.getValue(i);
//        }
//        return new Vector(result);
//    }
//
//    /**
//     * Вычитает другой вектор из текущего.
//     *
//     * @param other другой вектор
//     * @return результат вычитания в виде нового вектора
//     */
//    public Vector subtract(Vector other) {
//        if (n != other.getSize()) {
//            throw new IllegalArgumentException("Размеры векторов должны совпадать.");
//        }
//        double[] result = new double[n];
//        for (int i = 0; i < n; i++) {
//            result[i] = this.values[i] - other.getValue(i);
//        }
//        return new Vector(result);
//    }
//
//    /**
//     * Умножает вектор на скаляр.
//     *
//     * @param scalar значение скаляра
//     * @return новый вектор, умноженный на скаляр
//     */
//    public Vector multiply(double scalar) {
//        double[] result = new double[n];
//        for (int i = 0; i < n; i++) {
//            result[i] = this.values[i] * scalar;
//        }
//        return new Vector(result);
//    }
//
//    /**
//     * Возвращает евклидову норму (длину) вектора.
//     *
//     * @return евклидова норма вектора
//     */
//    public double norm() {
//        double sum = 0;
//        for (double v : values) {
//            sum += v * v;
//        }
//        return Math.sqrt(sum);
//    }
//
//    /**
//     * Возвращает минимальное значение вектора.
//     *
//     * @return минимальное значение вектора
//     */
//    public double min() {
//        double minValue = values[0];
//        for (double v : values) {
//            if (v < minValue) {
//                minValue = v;
//            }
//        }
//        return minValue;
//    }
//
//    /**
//     * Возвращает строковое представление вектора.
//     *
//     * @return строковое представление вектора
//     */
//    @Override
//    public String toString() {
//        return Arrays.toString(values);
//    }
//}
//
//class Matrix {
//    private int n;           // Количество строк
//    private int m;           // Количество столбцов
//    private double[][] values; // 2D массив для хранения элементов матрицы
//
//    /**
//     * Конструктор для создания матрицы заданного размера.
//     * Все элементы матрицы инициализируются нулем.
//     *
//     * @param n количество строк
//     * @param m количество столбцов
//     */
//    public Matrix(int n, int m) {
//        this.n = n;
//        this.m = m;
//        this.values = new double[n][m];
//    }
//
//    /**
//     * Конструктор для создания матрицы с заданными значениями.
//     *
//     * @param values 2D массив значений матрицы
//     */
//    public Matrix(double[][] values) {
//        this.n = values.length;
//        this.m = values[0].length;
//        this.values = new double[n][m];
//        for (int i = 0; i < n; i++) {
//            this.values[i] = Arrays.copyOf(values[i], m);
//        }
//    }
//
//    public Matrix(int size) {
//        this.values = new double[size][size];
//        for (int i = 0; i < size; i++) {
//            this.values[i][i] = 1;
//        }
//    }
//
//    /**
//     * Получает значение элемента матрицы по указанным индексам.
//     *
//     * @param i индекс строки
//     * @param j индекс столбца
//     * @return значение элемента матрицы
//     */
//    public double getValue(int i, int j) {
//        if (i < 0 || i >= n || j < 0 || j >= m) {
//            throw new IndexOutOfBoundsException("Индексы вне диапазона.");
//        }
//        return values[i][j];
//    }
//
//    /**
//     * Устанавливает значение элемента матрицы по указанным индексам.
//     *
//     * @param i     индекс строки
//     * @param j     индекс столбца
//     * @param value новое значение элемента
//     */
//    public void setValue(int i, int j, double value) {
//        if (i < 0 || i >= n || j < 0 || j >= m) {
//            throw new IndexOutOfBoundsException("Индексы вне диапазона.");
//        }
//        values[i][j] = value;
//    }
//
//    /**
//     * Возвращает количество строк в матрице.
//     *
//     * @return количество строк
//     */
//    public int getNumRows() {
//        return n;
//    }
//
//    /**
//     * Возвращает количество столбцов в матрице.
//     *
//     * @return количество столбцов
//     */
//    public int getNumCols() {
//        return m;
//    }
//
//    /**
//     * Транспонирует текущую матрицу.
//     *
//     * @return транспонированная матрица
//     */
//    public Matrix transpose() {
//        double[][] transposed = new double[m][n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < m; j++) {
//                transposed[j][i] = values[i][j];
//            }
//        }
//        return new Matrix(transposed);
//    }
//
//    public Matrix inverse() {
//        // Check that the matrix is square
//        if (n != m) {
//            throw new IllegalArgumentException("Matrix must be square.");
//        }
//
//
//        double[][] augmentedMatrix = new double[n][2 * n];
//
//        // Create augmented matrix [A | I]
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                augmentedMatrix[i][j] = values[i][j];
//            }
//            augmentedMatrix[i][i + n] = 1; // Fill identity matrix
//        }
//
//        // Forward elimination
//        for (int i = 0; i < n; i++) {
//            // Find the maximum element in the current column
//            double maxEl = Math.abs(augmentedMatrix[i][i]);
//            int maxRow = i;
//            for (int k = i + 1; k < n; k++) {
//                if (Math.abs(augmentedMatrix[k][i]) > maxEl) {
//                    maxEl = Math.abs(augmentedMatrix[k][i]);
//                    maxRow = k;
//                }
//            }
//
//            // Swap maximum row with current row
//            double[] temp = augmentedMatrix[maxRow];
//            augmentedMatrix[maxRow] = augmentedMatrix[i];
//            augmentedMatrix[i] = temp;
//
//            // Convert to upper triangular form
//            for (int k = i + 1; k < n; k++) {
//                double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
//                for (int j = 0; j < 2 * n; j++) {
//                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
//                }
//            }
//        }
//
//        // Backward elimination
//        for (int i = n - 1; i >= 0; i--) {
//            // Normalize the current row
//            double leadingElement = augmentedMatrix[i][i];
//            for (int j = 0; j < 2 * n; j++) {
//                augmentedMatrix[i][j] /= leadingElement;
//            }
//
//            // Eliminate above rows
//            for (int k = i - 1; k >= 0; k--) {
//                double factor = augmentedMatrix[k][i];
//                for (int j = 0; j < 2 * n; j++) {
//                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
//                }
//            }
//        }
//
//        // Extract the inverse matrix from the augmented matrix
//        double[][] inverseMatrix = new double[n][n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                inverseMatrix[i][j] = augmentedMatrix[i][j + n];
//            }
//        }
//
//        // Return the inverse as a new Matrix object
//        return new Matrix(inverseMatrix);
//    }
//
//    /**
//     * Складывает текущую матрицу с другой матрицей.
//     *
//     * @param other другая матрица
//     * @return результат сложения в виде новой матрицы
//     */
//    public Matrix add(Matrix other) {
//        if (n != other.getNumRows() || m != other.getNumCols()) {
//            throw new IllegalArgumentException("Размеры матриц должны совпадать.");
//        }
//        double[][] result = new double[n][m];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < m; j++) {
//                result[i][j] = this.values[i][j] + other.getValue(i, j);
//            }
//        }
//        return new Matrix(result);
//    }
//
//    /**
//     * Вычитает другую матрицу из текущей.
//     *
//     * @param other другая матрица
//     * @return результат вычитания в виде новой матрицы
//     */
//    public Matrix subtract(Matrix other) {
//        if (n != other.getNumRows() || m != other.getNumCols()) {
//            throw new IllegalArgumentException("Размеры матриц должны совпадать.");
//        }
//        double[][] result = new double[n][m];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < m; j++) {
//                result[i][j] = this.values[i][j] - other.getValue(i, j);
//            }
//        }
//        return new Matrix(result);
//    }
//
//    /**
//     * Умножает текущую матрицу на скаляр.
//     *
//     * @param scalar значение скаляра
//     * @return новая матрица, умноженная на скаляр
//     */
//    public Matrix multiply(double scalar) {
//        double[][] result = new double[n][m];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < m; j++) {
//                result[i][j] = this.values[i][j] * scalar;
//            }
//        }
//        return new Matrix(result);
//    }
//
//    /**
//     * Умножает текущую матрицу на другой вектор.
//     *
//     * @param vector другой вектор
//     * @return новый вектор, являющийся результатом умножения
//     */
//    public Vector multiply(Vector vector) {
//        if (m != vector.getSize()) {
//            throw new IllegalArgumentException("Число столбцов матрицы должно совпадать с размером вектора.");
//        }
//        double[] result = new double[n];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < m; j++) {
//                result[i] += this.values[i][j] * vector.getValue(j);
//            }
//        }
//        return new Vector(result);
//    }
//
//    /**
//     * Умножает текущую матрицу на другую матрицу.
//     *
//     * @param other другая матрица
//     * @return новая матрица, являющаяся результатом умножения
//     */
//    public Matrix multiply(Matrix other) {
//        if (m != other.getNumRows()) {
//            throw new IllegalArgumentException("Число столбцов первой матрицы должно совпадать с числом строк второй матрицы.");
//        }
//        double[][] result = new double[n][other.getNumCols()];
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < other.getNumCols(); j++) {
//                for (int k = 0; k < m; k++) {
//                    result[i][j] += this.values[i][k] * other.getValue(k, j);
//                }
//            }
//        }
//        return new Matrix(result);
//    }
//
//    /**
//     * Возвращает строковое представление матрицы.
//     *
//     * @return строковое представление матрицы
//     */
//    @Override
//    public String toString() {
//        StringBuilder sb = new StringBuilder();
//        for (double[] row : values) {
//            sb.append(Arrays.toString(row)).append("\n");
//        }
//        return sb.toString();
//    }
//}