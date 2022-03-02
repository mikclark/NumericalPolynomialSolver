from math import sqrt, prod, exp, log
from random import random, gauss
import time

EPS = 1e-13

def is_near(a:float, b:float, eps:float=EPS) -> bool:
    if b == 0.0:
        return abs(a) < eps
    elif a == 0.0:
        return abs(b) < eps
    return abs(float(a)/b - 1.0) < eps

def bisection_method(f, x0:float, x1:float) -> float:
    f0 = f(x0)
    f1 = f(x1)
    assert f0*f1 < 0, "Input x values must straddle the solution"
    x2 = x0
    dx = x1 - x0
    for i in range(52):
        dx *= 0.5
        x3 = x2 + dx
        if x3 == x2:
            break
        f3 = f(x3)
        if f3 == 0:
            return x3
        if f0*f3 > 0:
            x2 = x3
    return x2

def secant_method(f, x0:float, x1:float, eps:float=EPS) -> float:
    f0 = f(x0)
    f1 = f(x1)
    assert f0*f1 < 0, "Input x values must straddle the solution"
    if f0 - f1 == 0.0:
        x1 = 0.5*(x0+x1)
        f1 = f(x1)
        assert f0 - f1 != 0.0, "No solution for this function"
    iterations = []    
    for i in range(1000):
        x2 = x0 - f0 * (x0 - x1)/(f0 - f1)
        f2 = f(x2)
        iterations.append([(x0,f0),(x1,f1),':',(x2,f2)])
        if f2 == 0.0:
            return x2
        elif x2 == x1:
            x0 = 0.5*(x0 + x1)
            f0 = f(x0)
        elif x2 == x0:
            x1 = 0.5*(x0 + x1)
            f1 = f(x1)
        elif f2 * f0 < 0:
            x1 = x2
            f1 = f2
        else:
            x0 = x2
            f0 = f2
        if is_near(x1,x0,eps):
            return x2
    iter_str = '\n'.join(str(x) for x in iterations)
    raise Exception(f"Could not find solution. Iterations = \n{iter_str}")

def solve_quadratic(c : float, b : float, a : float) -> list[float]:
    if a == 0.0:
        if b == 0.0:
            return []
        return [-float(c)/b]
    determ = float(b*b) - 4.0*a*c
    const = -float(b)/2.0/a
    if determ < 0.0:
        return []
    elif determ == 0.0:
        return [const,const]
    sq_determ = sqrt(determ)/2.0/a
    return [
        const + sq_determ,
        const - sq_determ
        ]

def eval_poly(coefficients, x):
    assert isinstance(coefficients,list) and isinstance(x, (int,float)), \
           "Inputs must be a list and a value"
    xpow = 1.0
    total = 0.0
    for c in coefficients:
        total += c*xpow
        xpow *= x
    return total

def add_polys(p1 : list[float], p2 : list[float]) -> list[float]:
    p3 = [a1+a2 for a1,a2 in zip(p1,p2)]
    if len(p1) < len(p2):
        p3.extend( p2[len(p1):] )
    elif len(p1) > len(p2):
        p3.extend( p1[len(p2):] )
    return p3
        

def multiply_polys(p1 : list[float], p2 : list[float], *args : list[list[float]]) -> list[float]:
    result_poly = []
    if not p1 or not p2:
        return result_poly
    for i, coeff1 in enumerate(p1):
        temp_poly = [0 for ii in range(i)]
        temp_poly.extend(coeff1*coeff2 for coeff2 in p2)
        result_poly = add_polys(result_poly, temp_poly)

    for i in args:
        result_poly = multiply_polys(result_poly, i)
    return result_poly


def differentiate_poly(coefficients):
    return [coefficients[i]*i for i in range(1,len(coefficients))]

def solve_poly(coefficients : list[float], eps:float=EPS) -> list[float]:
    assert isinstance(coefficients, list), "Input must be a list of numbers"

    # End points of the recursion
    if len(coefficients) < 2:
        return []
    elif len(coefficients) == 2:
        b, m = coefficients
        if m == 0.0:
            return []
        return [-float(b)/m]
    elif len(coefficients) == 3:
        return solve_quadratic(*coefficients)

    # Recursive algorithm: find extrema first, then solve numerically
    f = lambda xx : eval_poly(coefficients, xx)
    limit_pos_inf = 1.0 if coefficients[-1] > 0.0 else -1.0
    limit_neg_inf = 1.0 if len(coefficients) % 2 == 1 or coefficients[-1] < 0.0 else -1.0
    extrema = sorted(solve_poly(differentiate_poly(coefficients)))
    extreme_values = [f(x) for x in extrema]
    #print("Extrema:\n" + "\n".join(
    #    f"   f({x}) = {y}" for x,y in zip(extrema,extreme_values)
    #    ))
    
    solutions = []
    solutions.extend(x for x,y in zip(extrema, extreme_values)
                     if y == 0.0)

    # Solve (or abandon) left-most and right-most values
    x_extr_left = extrema[0]
    y_extr_left = extreme_values[0]
    if y_extr_left != 0.0 and limit_neg_inf * y_extr_left < 0:
        dx = 2.0
        for i in range(1000):
            if f(x_extr_left - dx) * y_extr_left < 0:
                break
            dx *= 2.0
        x0 = x_extr_left - dx
        #print(f"Left ({x0},{x_extr_left})")
        x_solution_left = bisection_method(f, x0, x_extr_left)
        solutions.append(x_solution_left)
    
    x_extr_right = extrema[-1]
    y_extr_right = extreme_values[-1]
    if y_extr_right != 0.0 and limit_pos_inf * y_extr_right < 0:
        dx = 2.0
        for i in range(1000):
            if f(x_extr_right + dx) * y_extr_right < 0:
                break
            dx *= 2.0
        x0 = x_extr_right + dx
        #print(f"Right ({x0},{x_extr_right})")
        x_solution_right = bisection_method(f, x0, x_extr_right)
        solutions.append(x_solution_right)
    
    # Solve (or abandon) values between extrema
    for i in range(0, len(extrema)-1):
        x0 = extrema[i]
        x1 = extrema[i+1]
        if x0 == 0.0 or x1 == 0.0:
            continue
        #print(f"Mid {i+1}: ({x0},{x1})")
        x_soln = bisection_method(f, x0, x1)
        solutions.append(x_soln)

    return solutions


################################################################
# TESTING FUNCTIONS ONLY
################################################################

def test_is_near():
    print("...test_is_near")
    assert is_near(1.0, 1.0)
    assert is_near(1.0 + 1.0e-14, 1.0 + 3.0e-14)
    assert is_near(1.0 - 1.0e-14, 1.0 - 3.0e-14)
    assert not is_near(1.0 + 1.0e-13, 1.0 + 3.0e-13)
    assert not is_near(1.0 - 1.0e-13, 1.0 - 3.0e-13)
    assert is_near(3.0 + 1.0e-13, 3.0 + 3.0e-13)
    assert is_near(3.0 - 1.0e-13, 3.0 - 3.0e-13)
    assert not is_near(3.0 + 1.0e-13, 3.0 + 5.0e-13)
    assert not is_near(3.0 - 1.0e-13, 3.0 - 5.0e-13)

    assert is_near(0.0, 9e-14)
    assert is_near(0.0, -9e-14)
    assert is_near(9e-14, 0.0)
    assert is_near(-9e-14, 0.0)
    assert not is_near(0.0, 1.1e-13)
    assert not is_near(0.0, -1.1e-13)
    assert not is_near(1.1e-13, 0.0)
    assert not is_near(-1.1e-13, 0.0)

    from math import sin, cos, pi
    for i in range(10):
        assert is_near(sin(float(i)*pi), 0.0)
        assert is_near(cos(float(i)*pi), 1.0 if i%2==0 else -1.0)

def test_eval_poly():
    print("...test_eval_poly")
    for i in range(20):
        a = gauss(0.0, 100.0)
        x = gauss(0.0, 100.0)
        assert is_near(eval_poly([a], x), a)
        assert is_near(eval_poly([a, 0], x), a)
        assert is_near(eval_poly([a, 0,0,0], x), a)
        assert is_near(eval_poly([a, 0,0,0,0,0,0,0], x), a)
    for i in range(20):
        m = gauss(1.0, 100.0)
        b = gauss(-1.0, 100.0)
        assert is_near(eval_poly([b,m], 0.0), b)
        assert is_near(eval_poly([b,m,0,0], 0.0), b)
        assert is_near(eval_poly([b,m,0,0,0,0,0,0], 0.0), b)
        x = gauss(0.0, 100.0)
        assert is_near(eval_poly([b,m], x), b+m*x)
        assert is_near(eval_poly([b,m,0,0], x), b+m*x)
        assert is_near(eval_poly([b,m,0,0,0,0,0,0], x), b+m*x)
    for i in range(20):
        a = gauss(0.0, 100.0)
        b = gauss(0.0, 100.0)
        c = gauss(0.0, 100.0)
        assert is_near(eval_poly([c,b,a], 0.0), c)
        assert is_near(eval_poly([c,b,a,0,0,0], 0.0), c)
        for j in range(10):
            x = gauss(0.0, 10.0)
            y = c + b*x + a*x*x
            assert is_near(eval_poly([c,b,a], x), y)
            assert is_near(eval_poly([c,b,a,0,0,0], x), y)
    for i in range(20):
        a = gauss(0.0, 100.0)
        b = gauss(0.0, 100.0)
        c = gauss(0.0, 100.0)
        d = gauss(0.0, 100.0)
        e = gauss(0.0, 100.0)
        f = gauss(0.0, 100.0)
        assert is_near(eval_poly([a,b,c,d,e,f], 0.0), a)
        for j in range(10):
            x = gauss(0.0, 1.0)
            y = a + b*x + c*x*x + d*x*x*x + e*x*x*x*x + f*x*x*x*x*x
            assert is_near(eval_poly([a,b,c,d,e,f], x), y)
            assert is_near(eval_poly([a,b,c,d,e,f,0,0,0,0,0,0], x), y)

    x = gauss(3.0,1.0)
    for i in range(20):
        coeffs = [0]*i + [1.0]
        assert is_near(eval_poly(coeffs, x), pow(x,i))


def check_answers(list1, list2, eps=EPS):
        value = (
            len(list1) == len(list2) and
            all(
                sum(1 for j in list2 if is_near(i,j,eps)) ==
                sum(1 for k in list1 if is_near(i,k,eps))
                for i in list1
                )
            )
        return value

def test_differentiate_poly():
    print("...test_differentiate_poly")
    assert check_answers(differentiate_poly([1]), [])
    assert check_answers(differentiate_poly([0,1]), [1])
    assert check_answers(differentiate_poly([0,0,1]), [0,2])
    assert check_answers(differentiate_poly([3,4,2]), [4,4])
    assert check_answers(differentiate_poly([1,1,1,1,1,1]), [1,2,3,4,5])


def test_add_polys():
    def verify(a, b, c):
        c0 = add_polys(a,b)
        assert check_answers(c, c0), f"Mismatch: {c} != {c0}"
    print("...test_add_polys")
    verify([1],[2], [3])
    verify([0],[0], [0])
    verify([1,1],[2], [3,1])
    verify([1],[2,2], [3,2])
    verify([0,1,2,3,4,5,6],[6,5,4,3,2,1,0], [6]*7)


def test_multiply_polys_pairs():
    def verify(a, b, c):
        c0 = multiply_polys(a,b)
        assert check_answers(c, c0), f"Mismatch: {c} != {c0}"
    print("...test_multiply_polys_pairs")
    verify([1],[1], [1])
    verify([2],[2], [4])
    verify([0],[1], [0])
    verify([1],[0], [0])
    
    verify([1,1],[2], [2,2])
    verify([2],[1,1], [2,2])
    verify([1,1,1],[3], [3,3,3])
    verify([3],[1,1,1], [3,3,3])
    verify([1,2],[0], [0,0])
    verify([0],[1,2], [0,0])
    verify([1,2,3],[0], [0,0,0])
    verify([0],[1,2,3], [0,0,0])

    verify([],[1,2,3], [])
    verify([1,2,3],[], [])

    verify([1,1],[1,1], [1,2,1])
    verify([1,-1],[1,-1], [1,-2,1])
    verify([1,1],[1,-1], [1,0,-1])
    verify([1,-1],[1,1], [1,0,-1])
    verify([1,2],[1,2], [1,4,4])
    verify([1,-2],[1,-2], [1,-4,4])
    verify([1,2],[1,-2], [1,0,-4])
    verify([1,-2],[1,2], [1,0,-4])
    verify([1,3],[1,3], [1,6,9])
    verify([2,2],[2,2], [4,8,4])
    verify([2,3],[5,1], [10,17,3])

    verify([1,2,1],[1,1], [1,3,3,1])
    verify([1,2,1],[1,2,1], [1,4,6,4,1])
    verify([1,-2,1],[1,-1], [1,-3,3,-1])
    verify([1,-2,1],[1,-2,1], [1,-4,6,-4,1])


def test_multiply_polys_args():
    def verify(s1, s2):
        assert check_answers(s1,s2), f"Mismatch: {s1} != {s2}"
    verify(multiply_polys([1,1], [1,1], [1,1]),
           [1,3,3,1])
    verify(multiply_polys([1,-1], [1,-1], [1,-1]),
           [1,-3,3,-1])
    verify(multiply_polys([1,1], [1,1], [1,1], [1,1]),
           [1,4,6,4,1])
    verify(multiply_polys([1,-1], [1,-1], [1,-1], [1,-1]),
           [1,-4,6,-4,1])
    verify(multiply_polys([1,1], [1,1], [1,1], [1,1], [1,1]),
           [1,5,10,10,5,1])
    verify(multiply_polys([1,-1], [1,-1], [1,-1], [1,-1], [1,-1]),
           [1,-5,10,-10,5,-1])

    verify(multiply_polys([1,2,1], [1,1], [1,1]),
           [1,4,6,4,1])
    verify(multiply_polys([1,1], [1,2,1], [1,1]),
           [1,4,6,4,1])
    verify(multiply_polys([1,1], [1,1], [1,2,1]),
           [1,4,6,4,1])
    
    
def test_check_answers():
    print("...test_check_answers")
    x = 1e-14
    y = 1.1e-13
    assert check_answers([1],[1])
    assert check_answers([1,1],[1,1])
    assert check_answers([1,1,1,1],[1,1,1,1])
    assert check_answers([1,3],[1,3])
    assert check_answers([1,3],[3,1])
    assert check_answers([1,3,3],[3,1,3])
    assert check_answers([1,3,3],[3,3,1])
    assert check_answers([1,2,3,4,5,6],[2,4,6,1,3,5])
    assert check_answers([1,1,2,2,2],[2,1,2,1,2])

    assert not check_answers([1],[])
    assert not check_answers([], [1])
    assert not check_answers([1,2],[1,3])
    assert not check_answers([2,4],[3,5])
    assert not check_answers([1,1,1,1],[1,1,1])
    assert not check_answers([1,3,3],[1,1,3])
    assert not check_answers([1,1,3],[1,3,3])
    assert not check_answers([1,2,3,4,5,6],[1,2,3,4,5,7])
    
    assert check_answers([1.0, 2.0+x, 3.0-x],[1.0,2.0-x,3.0-x])
    assert not check_answers([1.0, 2.0+x, 3.0-x],[1.0+y,2.0-x,3.0-x])
    assert check_answers([1.0, 1.0+x, 1.0+x+x], [1.0, 1.0+x+x, 1.0+x+x+x+x])
    assert not check_answers([1.0, 1.0+x, 1.0+x+x], [1.0, 1.0+x+x, 1.0+y])


def test_bisection_method():
    print("...test_bisection_method")
    for i in range(1,100,2):
        f = lambda x: exp(x) - float(i)
        calc = bisection_method(f, -1.0, 10.0)
        real = log(float(i))
        assert is_near(calc, real), f"Mismatch: {real} != {calc}"
        

def test_secant_method():
    print("...test_secant_method")
    for i in range(20):
        b = 10.0*(random() - 0.5)
        m = 10.0*(random() - 0.5)
        x_soln = secant_method(lambda x: m*x + b, -50000.0, 50000.0)
        assert is_near(x_soln, -b/m)
    f = lambda x: -2.0 - x + x*x
    x_soln_1 = secant_method(f, 0.0, 50.0)
    assert is_near(x_soln_1, 2.0)
    x_soln_2 = secant_method(f, 0.0, -50.0)
    assert is_near(x_soln_2, -1.0)
        

def test_quadratic():
    print("...test_quadratic")
    #null
    assert check_answers(solve_quadratic(gauss(1.0,100.0),0.0,0.0), [])

    #linear
    assert check_answers(solve_quadratic(1.0,-1.0,0.0), [1.0])
    for i in range(20):
        c = gauss(1.0, 10.0)
        b = gauss(-1.0, 10.0)
        assert check_answers(solve_quadratic(c, b, 0.0), [-c/b])

    #quadratic
    assert check_answers(solve_quadratic(-9.0, 0.0, 1.0), [3.0, -3.0])
    assert check_answers(solve_quadratic(9.0, 0.0, 1.0), [])
    
    assert check_answers(solve_quadratic(-2.0, -1.0, 1.0), [2.0, -1.0])
    assert check_answers(solve_quadratic(-2.0, 1.0, 1.0), [-2.0, 1.0])
    assert check_answers(solve_quadratic(2.0, 1.0, -1.0), [2.0, -1.0])
    assert check_answers(solve_quadratic(2.0, -1.0, -1.0), [-2.0, 1.0])
    assert check_answers(solve_quadratic(4.0, -4.0, 1.0), [2.0,2.0])
    assert check_answers(solve_quadratic(4.0, 4.0, 1.0), [-2.0,-2.0])
    assert check_answers(solve_quadratic(1.0,1.0,1.0), [])
    assert check_answers(solve_quadratic(1.0,-1.0,1.0), [])
    phi = 1.618033988749895
    assert check_answers(solve_quadratic(-1.0,1.0,1.0), [-phi, 1.0/phi])
    assert check_answers(solve_quadratic(-1.0,-1.0,1.0), [phi, -1.0/phi])

def test_solve_poly_null():
    print("...test_solve_poly_null")
    assert len(solve_poly([])) == 0
    assert len(solve_poly([gauss(1.0,100.0)])) == 0

def test_solve_poly_linear():
    print("...test_solve_poly_linear")
    assert check_answers(solve_poly([1.0,-1.0]), [1.0])
    for i in range(20):
        b = gauss(1.0, 100.0)
        m = gauss(-1.0, 100.0)
        check_answers(solve_poly([b,m]), [-b/m])

def test_solve_poly_quad():
    print("...test_solve_poly_quad")
    assert check_answers(solve_poly([-9.0, 0.0, 1.0]), [3.0, -3.0])
    assert check_answers(solve_poly([9.0, 0.0, 1.0]), [])
    
    assert check_answers(solve_poly([-2.0, -1.0, 1.0]), [2.0, -1.0])
    assert check_answers(solve_poly([-2.0, 1.0, 1.0]), [-2.0, 1.0])
    assert check_answers(solve_poly([2.0, 1.0, -1.0]), [2.0, -1.0])
    assert check_answers(solve_poly([2.0, -1.0, -1.0]), [-2.0, 1.0])
    assert check_answers(solve_poly([4.0, -4.0, 1.0]), [2.0,2.0])
    assert check_answers(solve_poly([4.0, 4.0, 1.0]), [-2.0,-2.0])
    assert check_answers(solve_poly([1.0,1.0,1.0]), [])
    assert check_answers(solve_poly([1.0,-1.0,1.0]), [])
    phi = 1.618033988749895
    assert check_answers(solve_poly([-1.0,1.0,1.0]), [-phi, 1.0/phi])
    assert check_answers(solve_poly([-1.0,-1.0,1.0]), [phi, -1.0/phi])

def test_solve_poly_cubic():
    print("...test_solve_poly_cubic")
    # (x-1)*(x-2)*(x-3)
    calc = solve_poly([-6.0, 11.0, -6.0, 1.0])
    real = [1.0, 2.0, 3.0]
    assert check_answers(calc, real), f"Mismatch:\n{calc}\n{real}"

    # 2*(x-4.5)*(x+1.2)*(x+10.9)
    calc = solve_poly([-117.72, -82.74, 15.2, 2.0])
    real = [4.5, -1.2, -10.9]
    assert check_answers(calc, real), f"Mismatch:\n{calc}\n{real}"

def test_solve_poly_general():
    print("...test_solve_poly_general")
    n_polys = 0
    for n in range(5,16):
        print(f"    ... order {n}")
        t0 = time.time()
        for i in range(10):
            n_polys += 1
            solutions = [10.0*(random()-0.5) for i in range(n)]
            monomials = [[-s,1] for s in solutions]
            polynomial = multiply_polys(*monomials)
            calc = solve_poly(polynomial)
            assert check_answers(solutions, calc, eps=1e-6), \
                   f"""Mismatch:
real={sorted(solutions)},
calc={sorted(calc)},
poly={polynomial}"""
        t1 = time.time()
        print(f"    Time to solve 10 order-{n} polynomials: {t1-t0}")
        
    
if __name__ == "__main__":
    test_is_near()
    test_eval_poly()
    test_check_answers()
    test_differentiate_poly()
    test_add_polys()
    test_multiply_polys_pairs()
    test_multiply_polys_args()
    test_bisection_method()
    #test_secant_method()
    test_quadratic()
    test_solve_poly_null()
    test_solve_poly_linear()
    test_solve_poly_quad()
    test_solve_poly_cubic()
    test_solve_poly_general()
    print("All tests passed")
