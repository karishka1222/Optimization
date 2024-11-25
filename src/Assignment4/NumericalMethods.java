package Assignment4;

public class NumericalMethods {

    /** Demo examples
     * Bisection Method for Root-Finding - Function 1: f(x) = x^3 - 6x^2 + 11x - 6
     *  Golden Section Method for Unimodal Function Optimization - Function 2: f(x) = (x - 2)^2 + 3
     *  Gradient Ascent Method for Maximizing a Function - Function 3: f(x) = -x^2 + 4x + 1**/

    public static void main(String[] args) {
        //  Bisection Method for Root-Finding
        System.out.println(" Bisection Method for Root-Finding");
        double a1 = 0, b1 = 5, epsilon1 = 1e-6;
        double root = bisection(a1, b1, epsilon1);
        System.out.println("Approximate root: " + root);

        // Golden Section Method for Unimodal Function Optimization
        System.out.println("\nGolden Section Method for Unimodal Function Optimization");
        double a2 = 0, b2 = 5, epsilon2 = 1e-4;
        double xmin = goldenSection(a2, b2, epsilon2);
        System.out.println("Approximate xmin: " + xmin);
        System.out.println("f(xmin): " + f2(xmin));

        // Gradient Ascent Method for Maximizing a Function
        System.out.println("\nGradient Ascent Method for Maximizing a Function");
        double x0 = 0, alpha = 0.1;
        int iterations = 100;
        double xmax = gradientAscent(x0, alpha, iterations);
        System.out.println("Approximate xmax: " + xmax);
        System.out.println("f(xmax): " + f3(xmax));
    }



    //  Bisection Method for Root-Finding
    public static double f1(double x) {
        // Function f(x) = x^3 - 6x^2 + 11x - 6
        return Math.pow(x, 3) - 6 * Math.pow(x, 2) + 11 * x - 6;
    }

    public static double bisection(double a, double b, double epsilon) {
        if (f1(a) == 0) {
            return a; // a is root yet
        }
        if (f1(b) == 0) {
            return b; // b is root yet
        }
        if (f1(a) * f1(b) >= 0) {
            throw new IllegalArgumentException("f(a) and f(b) must not have the same sign.");
        }
        double c = a; // Средняя точка
        while ((b - a) >= epsilon) {
            c = (a + b) / 2;
            if (Math.abs(f1(c)) < epsilon) { // check stop condition
                break;
            }
            if (f1(c) * f1(a) < 0) {
                b = c;
            } else {
                a = c;
            }
        }
        return c;
    }



    // Golden Section Method for Unimodal Function Optimization
    public static double f2(double x) {
        // Function f(x) = (x - 2)^2 + 3
        return Math.pow(x - 2, 2) + 3;
    }

    public static double goldenSection(double a, double b, double epsilon) {
        double gr = (1 + Math.sqrt(5)) / 2;
        double c = b - (b - a) / gr;
        double d = a + (b - a) / gr;
        while (Math.abs(b - a) > epsilon) {
            if (f2(c) < f2(d)) {
                b = d;
            } else {
                a = c;
            }
            c = b - (b - a) / gr;
            d = a + (b - a) / gr;
        }
        return (a + b) / 2;
    }



    // Gradient Ascent Method for Maximizing a Function
    public static double f3(double x) {
        // Function f(x) = -x^2 + 4x + 1
        return -Math.pow(x, 2) + 4 * x + 1;
    }

    // нужно сделать динамический подсчёт производных
    public static double f3Prime(double x) {
        // Derivative of function - f'(x) = -2x + 4
        return -2 * x + 4;
    }

    public static double gradientAscent(double x0, double alpha, int iterations) {
        double x = x0;
        for (int i = 0; i < iterations; i++) {
            x += alpha * f3Prime(x); // update step
        }
        return x;
    }
}