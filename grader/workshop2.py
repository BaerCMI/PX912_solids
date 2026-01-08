import sympy as sym
import numpy as np
from sympy.matrices import Matrix, eye, zeros, ones, diag, MatMul

class Workshop2:
    print("Library Loaded!")
    def __init__(self, ):
        self.q1a = False
        self.q1b = False

        self.q2a = False
        self.q2b = False
        self.q2c = False
        
        self.q3a = False
        self.q3b = False
        
    # QUESTION 1 ------------------------------------------------------
    def hint1a(self,):
        print("Consider where the shape function has to take the values 1 and 0.")
        print("Inbetween the nodel the shape function is linear.")
    
    def check1a(self, N1, N2):
        X = sym.symbols('X')
        h = sym.Symbol('h', positive=True)
        if N1==sym.Piecewise((X/h, (0 <= X) & (X < h)), (2-X/h, (h <= X) & (X < 2*h)), (0, True)):
            if N2==sym.Piecewise((X/h-1, (h <= X) & (X < 2*h)), (3-X/h, (2*h <= X) & (X < 3*h)), (0, True)):
                print("\033[1;32m Correct!")
                self.q1a = True
            else:
                print("\033[0;31m Incorrect. Check N2.")
                self.q1a = False
        else:
            print("\033[0;31m Incorrect. Check N1.")
            self.q1a = False
    
    def hint1b(self,):
        print("Use function.diff(X)to dufferentiate with respect to X.")

    def check1b(self, dN1, dN2):
        X = sym.symbols('X')
        h = sym.Symbol('h', positive=True)
        if dN1==sym.Piecewise((1/h, (0 <= X) & (X < h)), (-1/h, (h <= X) & (X < 2*h)), (0, True)):
            if dN2==sym.Piecewise((1/h, (h <= X) & (X < 2*h)), (-1/h, (2*h <= X) & (X < 3*h)), (0, True)):
                print("\033[1;32m Correct!")
                self.q1b = True
            else:
                print("\033[0;31m Incorrect. Check dN2.")
                self.q1b = False
        else:
            print("\033[0;31m Incorrect. Check dN1.")
            self.q1b = False       
    
    # QUESTION 2 ------------------------------------------------------
    def hint2a(self,):
        print("You need to define a shape function for each of the 6 nodes (5 elements +1).")
        print("Differentiate the shape functions as before, and use sym.integrate() for the integral.")
        print("Fill the 6x6 stiffness matrix K according to the equation provided. ")
    
    def check2a(self, K):
        h = sym.Symbol('h', positive=True)
        E = sym.Symbol('E', positive=True)

        if K==sym.Matrix([[E/h, -E/h, 0, 0, 0, 0],
                          [-E/h, 2*E/h, -E/h, 0, 0, 0],
                          [0, -E/h, 2*E/h, -E/h, 0, 0],
                          [0, 0, -E/h, 2*E/h, -E/h, 0],
                          [0, 0, 0, -E/h, 2*E/h, -E/h],
                          [0, 0, 0, 0, -E/h, E/h]]):
            print("\033[1;32m Correct!")
            self.q2a = True
        else:
            print("\033[0;31m Incorrect.")
            self.q2a = False
    
    def hint2b(self,):
        print("Again use sym.integrate() for the integral and the shape functions defined in a).")
        print("f is a vector of length 6.")

    def check2b(self, f):
        h = sym.Symbol('h', positive=True)     # Element length
        rho = sym.Symbol('rho', positive=True) # Density
        g = sym.Symbol('g', positive=True)     # Gravitational acceleration

        if f==sym.Matrix([g*h*rho/2,
                          g*h*rho,
                          g*h*rho,
                          g*h*rho,
                          g*h*rho,
                          g*h*rho/2]):
            print("\033[1;32m Correct!")
            self.q2b = True
        else:
            print("\033[0;31m Incorrect.")
            self.q2b = False
    
    def hint2c(self,):
        print("Pass the reduced stiffness matrix and force vector to sym.linsolve()")
        print("together with a list of 5 symbols, one for each unknown displacement.")

    def check2c(self, d):
        h = sym.Symbol('h', positive=True)     # Element length
        E = sym.Symbol('E', positive=True)     # Young's modulus
        rho = sym.Symbol('rho', positive=True)
        g = sym.Symbol('g', positive=True)     # Gravitational acceleration

        if d==sym.FiniteSet(sym.Tuple(9*g*h**2*rho/(2*E), 8*g*h**2*rho/E, 21*g*h**2*rho/(2*E), 12*g*h**2*rho/E, 25*g*h**2*rho/(2*E))):
            print("\033[1;32m Correct!")
            self.q2c = True
        elif d==sym.Matrix([9*g*h**2*rho/(2*E), 8*g*h**2*rho/E, 21*g*h**2*rho/(2*E), 12*g*h**2*rho/E, 25*g*h**2*rho/(2*E)]):
            print("\033[1;32m Correct!")
            self.q2c = True
        else:
            print("\033[0;31m Incorrect.")
            self.q2c = False
        

    # QUESTION 3 ------------------------------------------------------
    def hint3a(self,):
        print("Integrate the polynomials from hand or using sym.integrate().")

    def check3a(self, p_exact_integral, q_exact_integral, r_exact_integral):
        if p_exact_integral==sym.Integer(26) or p_exact_integral==26:
            if q_exact_integral==sym.Integer(22) or q_exact_integral==22:
                if r_exact_integral==sym.Integer(393982) or r_exact_integral==393982:
                    print("\033[1;32m Correct!")
                    self.q3a = True
                else:
                    print("\033[0;31m Incorrect. Check r_exact_integral.")
                    self.q3a = False
            else:
                print("\033[0;31m Incorrect. Check q_exact_integral.")
                self.q3a = False
        else:
            print("\033[0;31m Incorrect. Check p_exact_integral.")
            self.q3a = False
    
    def hint3b(self,):
        print("Use scipy.integrate.quad() to perform the numerical integration.")
        print("Then compare the numerical and exact integrals to compute the error.")

    def check3b(self, p_absolute_errors, q_absolute_errors, r_absolute_errors):
        p_errors = np.array([6.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
        q_errors = np.array([16.0, 0.888888888888889, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
        r_errors = np.array([393834.0, 318511.466, 119088.932, 21985.7130, 2025.20597, 83.3415732, 1.17036129, 0.00237396493, 0.00])

        if np.allclose(p_absolute_errors, p_errors, rtol=1e-5, atol=1e-8):
            if np.allclose(q_absolute_errors, q_errors, rtol=1e-5, atol=1e-8):
                if np.allclose(r_absolute_errors, r_errors, rtol=1e-5, atol=1e-8):
                    print("\033[1;32m Correct!")
                    self.q3b = True
                else:
                    print("\033[0;31m Incorrect. Check r_absolute_error.")
                    self.q3b = False
            else:
                print("\033[0;31m Incorrect. Check q_absolute_error.")
                self.q3b = False
        else:
            print("\033[0;31m Incorrect. Check p_absolute_error.")
            self.q3b = False
    

    # RESULTS    ------------------------------------------------------
    
    def results(self):
        performance = [self.q1a,
                       self.q1b,  
                       self.q2a,
                       self.q2b,
                       self.q2c, 
                       self.q3a, 
                       self.q3b]
        
        names = ["Question 1a",
                 "Question 1b",
                 "Question 2a",
                 "Question 2b",
                 "Question 2c",
                 "Question 3a",
                 "Question 3b"]
        
        coupled = [(performance[i], names[i]) for i in range(len(performance))]
        
        total_score = sum(performance)
        
        print(f"Your score for this assignment: {total_score}/7")
        if total_score == 7:
            print("Excellent work!")
        else:
            print("The following questions have yet to be answered correctly:")
            for i in coupled:
                if i[0] == False:
                    print(i[1])
            print("Please feel free to reach out if you need any help.")
