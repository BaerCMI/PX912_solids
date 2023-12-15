import sympy as sym
import numpy as np
import math
from sympy.matrices import Matrix, eye, zeros, ones, diag, MatMul

class Workshop2:
    def __init__(self, ):
        print("Grader Library Loaded!")
        self.q1a = False
        self.q1b = False

        self.q2 = False
        
        self.q3a = False
        self.q3b = False
        
        self.q4a = False
        self.q4b = False
        self.q4c = False
        self.q4d = False
        
    # Question 1
    def hint1a(self,):
        print("Write each traction as a vector, then use them to build the stress tensor. What should zero traction along a direction look like?")  
        
    def check1a(self, stress):
        T_1, T_2 = sym.symbols('T_1 T_2')
        if stress == Matrix([[T_1, 0, 0], [0, T_2, 0], [0, 0, 0]]):
            print("\033[1;32m Correct!")
            self.q1a = True
        else:
            print("\033[0;31m Incorrect.")
            self.q1a = False    

    def hint1b(self,):
        print("This is a simpler case of the previous question. Use the same approach.")
        
    def check1b(self,stress):
        T = sym.symbols('T')
        if stress == Matrix([[T, 0, 0], [0, T, 0], [0, 0, 0]]):
            print("\033[1;32m Correct!")
            self.q1b = True
        else:
            print("\033[0;31m Incorrect.")
            self.q1a = False   
            
    # Question 2
    def hint2(self,):
        print("Tensor divergence is calculated by calculating the divergence of each row. How can this be implemented? You'll need a for loop. Your answer should return a vector.")
        
    def check2(self, div_sigma, b):
        rho, g, x_1, x_2, x_3 = sym.symbols('rho g x_1 x_2 x_3')
        if (div_sigma + b) == Matrix([[0], [0], [0]]):
            print("\033[1;32m Correct!")
            self.q2 = True
        else:
            print("\033[0;31m Incorrect.")
            self.q2 = False
        
    # Question 3
    def hint3a(self,):
        print("If you haven't written a routine for the infinitesimal strain yet, feel free to use mine: tools.infinitesimal_strain(phi,[X_1, X_2].")
    
    def check3a(self, voigt_eps):
        cond1 = math.isclose(voigt_eps[0], 0.000499999999999945)
        cond2 = math.isclose(voigt_eps[1], 0.000199999999999978)
        cond3 = math.isclose(voigt_eps[2], 0.00100000000000000)
                
        if cond1 and cond2 and cond3:
            print("\033[1;32m Correct!")
            self.q3a = True
        else:
            print("\033[0;31m Incorrect.")
            self.q3a = False

    def hint3b(self,):
        print("Check the last page of the lecture notes to remember plane")
        print("strain and plane stress.) You can represent the multiplication")
        print("of two matrices using the sym.MatMul(A, B) function. To actually")
        print("carry out the calculation? Write sym.MatMul(A, B).doit()")
        
    def check3b(self, cpe, cps):
        cond1 = math.isclose(cpe[0], 3901234.56790080)
        cond2 = math.isclose(cpe[1], 3012345.67901201)
        cond3 = math.isclose(cpe[2], 1481481.48148148)
        cond4 = math.isclose(cps[0], 2598290.59829031)
        cond5 = math.isclose(cps[1], 1709401.70940152)
        cond6 = math.isclose(cps[2], 1481481.48148148)        
                
        if cond1 and cond2 and cond3:
            if cond4 and cond5 and cond6:
                print("\033[1;32m Correct!")
                self.q3b = True
            else:
                print("\033[0;31m Incorrect. Check your answer for plane stress.")
                self.q3b = False

        else:
            print("\033[0;31m Incorrect. Check your answer for plane strain.")
            self.q3b = False

    
    # Question 4
    def hint4a(self,):
        print("Use xi.diff() like usual. What variable are you differentiating by this time?")
    def check4a(self,v1, v2, v3):
        t = sym.symbols('t')
        X = sym.Function('X')(t)
        Y = sym.Function('Y')(t)
        Z = sym.Function('Z')(t)
        cond1 = (v1 == 4*sym.Derivative(X, t))
        cond2 = (v2 == sym.Derivative(Y, t))
        cond3 = (v3 == sym.Derivative(Z, t))
        
        if cond1:
            if cond2:
                if cond3:
                    print("\033[1;32m Correct!")
                    self.q4a = True
                else:
                    print("\033[0;31m Incorrect. Check your answer for v3.")
                    self.q4a = False
            else:
                print("\033[0;31m Incorrect. Check your answer for v2.")
                self.q4a = False
        else:
            print("\033[0;31m Incorrect. Check your answer for v1.")
            self.q4a = False


    
    def hint4b(self,):
        print("Do you really need to calculate the second component of this expression? Why or why not?")
    def check4b(self, a1, a2, a3):
        t = sym.symbols('t')
        X = sym.Function('X')(t)
        Y = sym.Function('Y')(t)
        Z = sym.Function('Z')(t)
        cond1 = (a1 == 4*sym.Derivative(X, (t, 2)))
        cond2 = (a2 == sym.Derivative(Y, (t, 2)))
        cond3 = (a3 == sym.Derivative(Z, (t, 2)))
        
        if cond1:
            if cond2:
                if cond3:
                    print("\033[1;32m Correct!")
                    self.q4b = True
                else:
                    print("\033[0;31m Incorrect. Check your answer for a3.")
                    self.q4b = False
            else:
                print("\033[0;31m Incorrect. Check your answer for a2.")
                self.q4b = False
        else:
            print("\033[0;31m Incorrect. Check your answer for a1.")
            self.q4b = False
    
    def hint4c(self,):
        print("This should be similar to Question 2. You may need to make a change of symbol.")
    def check4c(self, divt):
        t = sym.symbols('t')
        X = sym.Function('X')(t)
        Y = sym.Function('Y')(t)
        Z = sym.Function('Z')(t)
        
        if divt == Matrix([[0], [10*Y], [-18*Z]]):
            print("\033[1;32m Correct!")
            self.q4c = True
        else:
            print("\033[0;31m Incorrect.")
            self.q4c = False


    
    def hint4d(self,):
        print("Check the form of Newton's second law and arrange all parts of the expression to equal zero.")
    def check4d(self,funx, funy, funz):
        t = sym.symbols('t')
        X = sym.Function('X')(t)
        Y = sym.Function('Y')(t)
        Z = sym.Function('Z')(t)
        cond1 = (funx == sym.Derivative(X, t) + 4*sym.Derivative(X, (t, 2)))
        cond2 = (funy == -10*Y + 2*sym.Derivative(Y, t) + sym.Derivative(Y, (t, 2)))
        cond3 = (funz == 18*Z + 5*sym.Derivative(Z, t) + sym.Derivative(Z, (t, 2)))
        
        if cond1:
            if cond2:
                if cond3:
                    print("\033[1;32m Correct!")
                    print("Well done for completing the question. Excellent work!")
                    self.q4d = True
                else:
                    print("\033[0;31m Incorrect. Check your answer for funx.")
                    self.q4d = False
            else:
                print("\033[0;31m Incorrect. Check your answer for funy.")
                self.q4d = False
        else:
            print("\033[0;31m Incorrect. Check your answer for funy.")
            self.q4d = False        
    
    
    def results(self):
        performance = [self.q1a,
                       self.q1b,  
                       self.q2,
                       self.q3a, 
                       self.q3b,
                       self.q4a,
                       self.q4b,
                       self.q4c,
                       self.q4d]
        
        names = ["Question 1a",
                 "Question 1b",
                 "Question 2",
                 "Question 3a",
                 "Question 3b",
                 "Question 4a",
                 "Question 4b",
                 "Question 4c",
                 "Question 4d"]
        
        coupled = [(performance[i], names[i]) for i in range(9)]
        
        total_score = sum(performance)
        
        print(f"Your score for this assignment: {total_score}/9")
        if total_score == 9:
            print("Excellent work!")
        else:
            print("The following questions have yet to be answered correctly:")
            for i in coupled:
                if i[0] == False:
                    print(i[1])
            print("Please feel free to reach out if you need any help.")