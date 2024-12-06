package Assignment4;

import java.util.function.Function;

public class NumericalMethods {

    public static void main(String[] args) {
        // Bisection Method for Root-Finding
        System.out.println("Bisection Method for Root-Finding");

        // Define the interval [a, b], tolerance (epsilon), and the function f(x)
        double a1 = 1, b1 = 2, epsilon1 = 1e-6;
        Function<Double, Double> f1 = x -> Math.pow(x, 5) - 6 * Math.pow(x, 2) + 11 * x - 6;

        // Call the bisection method to find the root
        double root = bisection(f1, a1, b1, epsilon1);
        System.out.println("Approximate root: " + root);
        // Golden Section Method for Unimodal Function Optimization
        System.out.println("\nGolden Section Method for Unimodal Function Optimization");

        // Define the interval [a, b], tolerance (epsilon), and the unimodal function f(x)
        double a2 = 0, b2 = 5, epsilon2 = 1e-4;
        Function<Double, Double> f2 = x -> Math.pow(x - 2, 2) + 3;

        // Call the golden section method to find the minimum
        double xmin = goldenSection(f2, a2, b2, epsilon2);
        System.out.println("Approximate xmin: " + xmin);
        System.out.println("f(xmin): " + f2.apply(xmin));
        // Gradient Ascent Method for Maximizing a Function
        System.out.println("\nGradient Ascent Method for Maximizing a Function");
        // Define the initial guess (x0), learning rate (alpha), number of iterations,
        // and the function's derivative f'(x)
        double x0 = 0, alpha = 0.1;
        int iterations = 100;

        Function<Double, Double> f3 = x -> -Math.pow(x, 2) + 4 * x + 1; // Original function
        Function<Double, Double> f3Prime = x -> -2 * x + 4; // Derivative of the function

        // Call the gradient ascent method to find the maximum
        double xmax = gradientAscent(f3Prime, x0, alpha, iterations);
        System.out.println("Approximate xmax: " + xmax);
        System.out.println("f(xmax): " + f3.apply(xmax));
    }

    // Bisection Method for Root-Finding
    public static double bisection(Function<Double, Double> f, double a, double b, double epsilon) {
        // Check if the endpoints are roots
        if (f.apply(a) == 0) {
            System.out.println("I = [" + a + ", " + b + "]");
            return a; // a is a root
        }
        if (f.apply(b) == 0) {
            System.out.println("I = [" + a + ", " + b + "]");
            return b; // b is a root
        }
        // Ensure the function has opposite signs at a and b
        if (f.apply(a) * f.apply(b) >= 0) {
            throw new IllegalArgumentException("f(a) and f(b) must not have the same sign.");
        }
        double c = a; // Middle point of the interval
        while ((b - a) >= epsilon) { // Iterate until the interval is smaller than epsilon
            c = (a + b) / 2; // Calculate the midpoint
            if (Math.abs(f.apply(c)) < epsilon) { // Stop if the root is found
                break;
            }
            if (f.apply(c) * f.apply(a) < 0) { // Narrow down the interval
                b = c;
            } else {
                a = c;
            }
        }
        System.out.println("I = [" + a + ", " + b + "]");
        return c; // Return the approximate root
    }

    // Golden Section Method for Unimodal Function Optimization
    public static double goldenSection(Function<Double, Double> f, double a, double b, double epsilon) {
        // Define the golden ratio
        double gr = (1 + Math.sqrt(5)) / 2;
        // Initialize points c and d inside the interval
        double c = b - (b - a) / gr;
        double d = a + (b - a) / gr;
        while (Math.abs(b - a) > epsilon) { // Iterate until the interval is smaller than epsilon
            if (f.apply(c) < f.apply(d)) { // Update the interval based on function comparisons
                b = d; // Narrow down to the left side
            } else {
                a = c; // Narrow down to the right side
            }
            // Recalculate points c and d
            c = b - (b - a) / gr;
            d = a + (b - a) / gr;
        }
        System.out.println("I = [" + a + ", " + b + "]");
        return (a + b) / 2; // Return the approximate minimum
    }

    // Gradient Ascent Method for Maximizing a Function
    public static double gradientAscent(Function<Double, Double> fPrime, double x0, double alpha, int iterations) {
        double x = x0; // Start at the initial guess
        for (int i = 0; i < iterations; i++) { // Iterate for a fixed number of steps
            x += alpha * fPrime.apply(x); // Update the value of x using the gradient
        }
        return x; // Return the approximate maximum
    }
}
