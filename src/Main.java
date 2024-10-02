import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;

public class Main {

    public static void main(String[] args) throws IOException {
        FileReader reader = new FileReader("./tests/test1.txt");
        // input
        Scanner scanner = new Scanner(reader);
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

        // output results
        if (true) {
            System.out.println("The method is not applicable!");
        } else {

        }
    }

    public static void simplexMethod(Vector C, Matrix A, Vector b, int accuracy) {

    }
}


class Vector {
    private int n;
    private int[] values;
    Vector(int n, int[] values) {
        this.n = n;
        this.values = values;
    }
}

class Matrix {
    private int n;
    private int m;
    private int[][] values;
    Matrix(int n, int m, int[][] values) {
        this.n = n;
        this.m = m;
        this.values = values;
    }
}