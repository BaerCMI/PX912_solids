import sympy as sym
import numpy as np
import math
from sympy.matrices import Matrix, eye, zeros, ones, diag, MatMul
import copy

class Workshop3:
    def __init__(self, ):
        print("Grader Library Loaded!")
        self.q1a = False
        self.q1c = False

        self.q2a = False
        self.q2b = False
        self.q2c = False
        self.q2d = False
        
    def hint1a(self,):
        print("Take y = A*x*x + B*x + C and apply these constraints to get the values of A, B, and C.")
        print("It helps to work this out by hand first!")
        
    def check1a(self, y1, y2,y3):
        x = sym.symbols('x')        
        name = ["y1", "y2", "y3"]
        answer = [y1,y2,y3]
        solution = [0.5*x**2 - 0.5*x, 
                   1.0 - x**2, 
                   0.5*x**2 + 0.5*x]
        altsolution = [sym.Rational(1,2)*x**2 - sym.Rational(1,2)*x, 
                   1 - x**2, 
                   sym.Rational(1,2)*x**2 + sym.Rational(1,2)*x]
        
        incorrect = []
        for i in range(3):
            if answer[i] != solution[i] and answer[i] != altsolution[i]:
                incorrect.append(name[i])
                
        if len(incorrect) == 0:
            print("\033[1;32m Correct!")
            self.q1a = True
        else:
            incorrect_str = ", ".join(incorrect)
            print(f"\033[0;31m Incorrect. Check your answers for: {incorrect_str}.")
            self.q1a = False

            
    def hint1b(self,):
        print("The correct format for the plotting is: p = sym.plotting.plot(y_1, y_2, y_3, (x,min,max))") 
    
    def hint1c(self,):
        print("The correct format for the integration is: integrand.integrate((x,0,L)) or sym.integrate(integrand, (x, 0, L)).")
    def check1c(self, k):
        A, E, L = sym.symbols('A E L')
        solution = Matrix([[0.333333333333333*A*E*L**3 - 0.5*A*E*L**2 + 0.25*A*E*L,
                            -0.666666666666667*A*E*L**3 + 0.5*A*E*L**2,
                            0.333333333333333*A*E*L**3 - 0.25*A*E*L],
                           [-0.666666666666667*A*E*L**3 + 0.5*A*E*L**2,
                            4*A*E*L**3/3,
                            -0.666666666666667*A*E*L**3 - 0.5*A*E*L**2],
                            [0.333333333333333*A*E*L**3 - 0.25*A*E*L,
                             -0.666666666666667*A*E*L**3 - 0.5*A*E*L**2,
                             0.333333333333333*A*E*L**3 + 0.5*A*E*L**2 + 0.25*A*E*L]])
        korig = copy.deepcopy(k)
        for a in sym.preorder_traversal(k):
            if isinstance(a, sym.Float):
                k = k.subs(a, round(a, 1))
        for a in sym.preorder_traversal(solution):
            if isinstance(a, sym.Float):
                solution = solution.subs(a, round(a, 1))
                
        if k == solution:
            print("\033[1;32m Correct!")
            self.q1c = True
        else:
          solution = Matrix([[sym.Rational(1,3)*A*E*L**3 - sym.Rational(1,2)*A*E*L**2 + sym.Rational(1,4)*A*E*L,
                              -sym.Rational(2,3)*A*E*L**3 + sym.Rational(1,2)*A*E*L**2,
                              sym.Rational(1,3)*A*E*L**3 - sym.Rational(1,4)*A*E*L],
                             [-sym.Rational(2,3)*A*E*L**3 + sym.Rational(1,2)*A*E*L**2,
                              sym.Rational(4,3)*A*E*L**3,
                              -sym.Rational(2,3)*A*E*L**3 - sym.Rational(1,2)*A*E*L**2],
                              [sym.Rational(1,3)*A*E*L**3 - sym.Rational(1,4)*A*E*L,
                               -sym.Rational(2,3)*A*E*L**3 - sym.Rational(1,2)*A*E*L**2,
                               sym.Rational(1,3)*A*E*L**3 + sym.Rational(1,2)*A*E*L**2 + sym.Rational(1,4)*A*E*L]])
          if korig == solution:
              print("\033[1;32m Correct!")
              self.q1c = True
          else: 
            print(f"\033[0;31m Incorrect.")
            self.q1c = False


    
    # -----------------------------------------------------------------------
    # Problem 2
    def hint2a(self,):
        print("Use np.ix_(indices, indices) to find the correct indices of the")
        print("global K matrix that need to be filled with the sum of entries of the smaller K matrices.")
        
    def check2a(self, K_global):
        answer = [[ 200000.,       0., -200000.],
                  [      0.,   50000.,  -50000.], 
                  [-200000.,  -50000.,  250000.]]
        if np.array_equal(K_global, answer):
            print("\033[1;32m Correct!")
            self.q2a = True
        else:
            print(f"\033[0;31m Incorrect.")
            self.q2a = False


    def hint2b(self,):
        print("What portion of the matrix do you need? That is, what index of K_global corresponds to the node where a force is acting? What do you know about walls?")
    def check2b(self, d_3):
        if d_3 == 4e-05:
            print("\033[1;32m Correct!")
            self.q2b = True
        else:
            print(f"\033[0;31m Incorrect.")
            self.q2b = False           
 

    def hint2c(self,):
        print("For this we have to multiply the portion of the global stiffness matrix that we need and the displacement vector that we have already found.") 
    def check2c(self, R_F):
        solution = [-8., -2.]
        if np.array_equal(solution, R_F):
            print("\033[1;32m Correct!")
            self.q2c = True
        else:
            print(f"\033[0;31m Incorrect.")
            self.q2c = False           

    def hint2d(self,):
        print("For this, you need to find the derivatives of the 2-nodal shape functions for each individual element.") 
        print("You can use those to find the elemental strains. These, in turn, can be used to find the elemental stresses.") 
    def check2d(self, eps_1, eps_2, eps_3, sigma_1, sigma_2, sigma_3):
        name = ["eps_1",
                "eps_2",
                "eps_3",
                "sigma_1",
                "sigma_2",
                "sigma_3"]
        
        answer = [eps_1, eps_2, eps_3, sigma_1, sigma_2, sigma_3]
        solution = [0.0004, 0.0004, -0.0008, 4000000.0, 4000000.0, -4000000.0]
        incorrect = []
        
        for i in range(6):
            if answer[i] != solution[i]:
                incorrect.append(name[i])
                
        if len(incorrect) == 0:
            print("\033[1;32m Correct!")
            self.q2d = True
        else:
            incorrect_str = ", ".join(incorrect)
            print(f"\033[0;31m Incorrect. Check your answers for: {incorrect_str}.")
            self.q2d = False
            
            
    def results(self):
        performance = [self.q1a,
                       self.q1c,
                       self.q2a, 
                       self.q2b,
                       self.q2c,
                       self.q2d]
        
        names = ["Question 1a",
                 "Question 1c",
                 "Question 2a",
                 "Question 2b",
                 "Question 2c",
                 "Question 2d"]

        coupled = [(performance[i], names[i]) for i in range(6)]
        
        total_score = sum(performance)
        
        print(f"Your score for this assignment: {total_score}/6")
        if total_score == 6:
            print("Excellent work!")
        else:
            print("The following questions have yet to be answered correctly:")
            for i in coupled:
                if i[0] == False:
                    print(i[1])
            print("Please feel free to reach out if you need any help.")
