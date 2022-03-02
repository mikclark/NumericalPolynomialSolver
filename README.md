# NumericalPolynomialSolver
For a set of polynomial coefficients, `polynomial_solver.solve_poly` finds all real roots of the polynomial.

## HOW TO

Any *N*-degree polynomial *P(x)* can be expressed as follows:

*P(x) = a₀ + a₁x + a₂x² + a₃x³ + ... aₙxⁿ*

Therefore, a polynomial can be completely described by the array [*a₀, a₁, a₂, a₃, ..., aₙ*]. This is the input to the function `solve_poly` in this file. `solve_poly` will numerically search for all real roots of the polynomial (all *x* where *P(x)* = 0).

## Code Example

```
from polynomial_solver import solve_poly

linear_solutions = solve_poly([1,1])
# linear_solutions = [-1]

quadratic_solutions = solve_poly([1,-3,2])
# quadratic_solutions = [-1, -2]
```

## Methodology

`solve_poly` is a recursive function. For linear or quadratic inputs, the solutions are found explicitly (using either *x = -b/m* or the quadratic formula). For polynomials of degree 3 or larger, the first step is to determine all extreme values (i.e., [maxima and minima](https://en.wikipedia.org/wiki/Maxima_and_minima)) of the polynomial, which is found by calculating the array representing the [derivative](https://en.wikipedia.org/wiki/Derivative) of *P(x)* and calling `solve_poly` on the derivative. Then, solutions are sought around the extreme values: at most one solution to *P(x)* might lie to the left of the smallest-*x* extreme value, at most one value might lie to the right of the largest-*x* extreme value, and there is at most one solution between each pair of extreme values. The [bisection method](https://en.wikipedia.org/wiki/Bisection_method) is used to find solutions between extreme values, down to [machine epsilon](https://en.wikipedia.org/wiki/Machine_epsilon).

## Other functions that come for free

- bisection_method : accepts a function and two x-values
- secant_method : accepts a function and two x-values
- solve_quadratic : accepts 3 floats, representing the coefficients of a quadratic polynomial, and uses the quadratic formula to find the real solutions
- eval_poly : accepts a list of polynomial coefficients and a x-value. Evaluates the polynomial at x.
- add_polys : accepts two lists of polynomial coefficients and returns the coefficients of the polynomial that represents their sum
- multiply_polys : accepts two lists of polynomial coefficients and returns the coefficients of the polynomial that represents their product
- differentiate_poly : accepts a list of polynomial coefficients and returns the coefficients of the polynomial that represents its derivative
