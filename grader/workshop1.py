import sympy as sym
from sympy.matrices import Matrix, eye, zeros, ones, diag, MatMul

class Workshop1:
    print("Library Loaded!")
    def __init__(self, ):
        self.q1a = False
        self.q1b = False

        self.q2a = False
        self.q2b = False
        self.q2c = False
        
        self.q3a = False
        self.q3b = False

        self.q4b = False
        self.q4c = False
        self.q4d = False
        
    # QUESTION 1 ------------------------------------------------------
    def hint1a(self,):
        print("To fill the deformation gradient matrix, you can use two")
        print("loops, one over phi(x) and one over the list of symbols (X)")
        print("Phi contains the definition of x1, x2 and x3.")
        
    def check1a(self, F_a, homogeneous):
        X_1, X_2, X_3 = sym.symbols('X_1 X_2 X_3')
        if F_a == Matrix([[1, 0, 0],
                          [0, X_3, X_2],
                          [0, 0, 1]]):
            if homogeneous==False:
                print("\033[1;32m Correct!")
                self.q1a = True
            else:
                print("\033[0;31m Incorrect. Recall the definition of homogeneity.")
                self.q1a = False    
        else:
            print("\033[0;31m Incorrect. Check your answer for the deformation gradient.")
            self.q1a = False
        
    def hint1b(self,):
        print("This is very similar to the first question. Chance phi according to")
        print("the defintion provided. The symbold X_1, X_2 and X_3 stay the same.")        
    def check1b(self, F_b, homogeneous):
        if F_b == Matrix([[0, 2, 0], [0, 0, 1], [5, 0, 0]]):
            if homogeneous==True:
                print("\033[1;32m Correct!")
                self.q1b = True
            else:
                print("\033[0;31m Incorrect. Recall the definition of homogeneity.")
                self.q1b = False  
        else:
            print("\033[0;31m Incorrect. Check your answer for the deformation gradient.")
            self.q1b = False           


    # QUESTION 2 ------------------------------------------------------
    def hint2a(self, ):
        print("The main difference here from above is the number of entries of phi.")
        print("You will have the make an additional change to the input of the")
        print("deformation gradient function as only two material coordinates are needed.")

    def check2a(self, F_2):
        if F_2 == Matrix([[1.00000000000000, -1.50000000000000], 
                          [0, -1.50000000000000]]):
            print("\033[1;32m Correct!")
            self.q2a = True
        else:
            print("\033[0;31m Incorrect.")
            self.q2a = False             
    
    def hint2b(self):
        print("The deformation expression can be thought of as a translation")
        print("rule. The input in undeformed coordinates will be the output")
        print("in deformed coordinates.")
        print("Applying the deformation gradient will account for shape changes,")
        print("but the translation of the origin during the deformation also needs")
        print("to be taken into account.")

    def check2b(self, deformed_corners):
        if deformed_corners == sym.Matrix([[4.0, 2.0],
                                             [7.0, 5.0],
                                             [5.0, 5.0],
                                             [2.0, 2.0]]).T:
            print("\033[1;32m Correct!")
            self.q2b = True
        else:
            print("\033[0;31m Incorrect.")
            print("Note: if you've modified the order of the square_corners list, you may be getting an error due to the order of your points.")
            self.q2b = False                 

    
    def hint2c(self,):
        print("This question requires application of the deformation gradient to E_1 and E_2,")
        print("but e_1 and e_2 require the inverse treatment.")
        print("Sympy provides the .inv() method.")

        
    def check2c(self, e1, e2, E1, E2):
        ce1, ce2, cE1, cE2 = False, False, False, False
        if e1 == sym.Matrix([1.0,0.0]):
            ce1=True
        else:
            ce1=False
            print("\033[0;31m def_e1 is incorrect.")
            self.q2b = False 
           
        if e2 == sym.Matrix([-1.0,-2/3]):
            ce2=True
        else:
            ce2=False
            print("\033[0;31m def_e2 is incorrect.")
            self.q2b = False 
                    
        if E1 == sym.Matrix([1.0,0.0]):
            cE1=True
        else:
            cE1=False   
            print("\033[0;31m def_E1 is incorrect.")
            self.q2b = False 
            
        if E2 == sym.Matrix([-1.5,-1.5]):
            cE2=True
        else:
            cE2=False
            print("\033[0;31m def_E2 is incorrect.")
            self.q2b = False 
                        
        if ce1 and ce2 and cE1 and cE2:
            print("\033[1;32m Correct!")
            self.q2c = True
            
            
    # QUESTION 3 ------------------------------------------------------
    def hint3a(self):
        print("Recall the equation for the infinitesimal strain. What are the ")
        print("separate operations it demands?")
        print("One way of building an identity matrix in SymPy is using the eye")
        print("method, Matrix.eye(N). N is the size of the identity matrix.")
        print("You can calculate the transpose of the matrix A by simply")
        print("writing A.T")

    
    def check3a(self, ifts):
        alpha = sym.symbols('alpha')
        if ifts == Matrix([[0, 0.5*alpha, 0], 
                           [0.5*alpha, 0, 0],
                           [0, 0, 0]]):
            print("\033[1;32m Correct!")
            self.q3a = True           
        elif ifts == Matrix([[0, sym.Rational(1,2)*alpha, 0], 
                           [sym.Rational(1,2)*alpha, 0, 0],
                           [0, 0, 0]]):
            print("\033[1;32m Correct!")
            self.q3a = True           
        else:
            print("\033[0;31m Incorrect.")
            self.q3a = False   
     
    def hint3b(self):
        print("Isochoric = isos + khora = equal + space.")
        print("What does it mean for a a transformation to have equal space?")
        print("Note that in this situation, space is volume.")

    def check3b(self, isochoric):
        if isochoric == True:
            print("\033[1;32m Correct!")
            self.q3b = True           
        else:
            print("\033[0;31m Incorrect.")
            self.q3b = False   
            

    # def hint3c(self):
    #     print("This should be very similar to a previous problem.")
    #     print("Simply translate the points that you are given.")
      
    # def check3c(self, deformed_corners):
    #     if deformed_corners == sym.Matrix([[2.5, 1.0],
    #                                         [-0.5, -1.0],
    #                                         [-2.5, -1.0],
    #                                         [0.5, 1.0]]).T:
    #         print("\033[1;32m Correct!")
    #         self.q3c = True
    #     else:
    #         print("\033[0;31m Incorrect.")
    #         print("Note: if you've modified the order of the square_corners list, you may be getting an error due to the order of your points.")
    #         self.q3c = False    


    # QUESTION 4 ------------------------------------------------------ 
    def hint4a(self):
        print("Use Matrix.trace() and Matrix.eye(N), where N is ")
        print("the size of the identity matrix.")
    
    def hint4b(self):
        print("Apply the infinitesimal strain tensor from Q3 to the ")
        print("function created in 4a.")

    def check4b(self, S):
        E, alpha, nu = sym.symbols('E alpha nu')
        if S == sym.Matrix([[0, (0.5*E*alpha)/(nu+1), 0],
                            [0.5*E*alpha/(nu+1), 0, 0],
                            [0, 0, 0]]):
            print("\033[1;32m Correct!")
            self.q4b = True
        elif S == sym.Matrix([[0, sym.Rational(1,2)*E*alpha/(nu+1), 0],
                            [sym.Rational(1,2)*E*alpha/(nu+1), 0, 0],
                            [0, 0, 0]]):
            print("\033[1;32m Correct!")
            self.q4b = True
        else:
            print("\033[0;31m Incorrect.")
            self.q4b = False
    
    def hint4c(self):
        print("Create a sym.Matrix of the rotation matrix R.")
        print("Use sym.cos() and sym.sin() for cosine and sine functions.")

    def check4c(self, S_rotated):
        E, nu, theta = sym.symbols('E nu theta')
        if S_rotated == sym.Matrix([[E*nu*(2*sym.cos(theta)+1)/((1-2*nu)*(nu+1)) + E*sym.cos(theta)/(nu+1), - (E*sym.sin(theta))/(nu+1), 0],
                                    [E*sym.sin(theta)/(nu+1), E*nu*(2*sym.cos(theta)+1)/((1-2*nu)*(nu+1)) + E*sym.cos(theta)/(nu+1), 0],
                                    [0, 0, E*nu*(2*sym.cos(theta)+1)/((1-2*nu)*(nu+1)) + E/(nu+1)]]):
            print("\033[1;32m Correct!")
            self.q4c = True
        else:
            print("\033[0;31m Incorrect.")
            self.q4c = False

    def hint4d(self):
        print("Apply the expr.Expr.series(x, xo, n) method to S11.")
        print("Here, x is the variable you want to expand about, ")
        print("xo is the value around which x is calculated, ")
        print("and n is the order of the expansion.")
    
    def check4d(self, S11_exp):
        E, nu, theta = sym.symbols('E nu theta')
        if S11_exp == E/(nu+1)+(3*E*nu)/((1-2*nu)*(nu+1)) + sym.O(theta**2):
            print("\033[1;32m Correct!")
            self.q4d = True
        else:
            print("\033[0;31m Incorrect.")
            self.q4d = False


    # RESULTS    ------------------------------------------------------
    
    def results(self):
        performance = [self.q1a,
                       self.q1b,  
                       self.q2a,
                       self.q2b,
                       self.q2c, 
                       self.q3a, 
                       self.q3b, 
                       #self.q3c,
                       self.q4b,
                       self.q4c,
                       self.q4d]
        
        names = ["Question 1a",
                 "Question 1b",
                 "Question 2a",
                 "Question 2b",
                 "Question 2c",
                 "Question 3a",
                 "Question 3b",
                # "Question 3c",
                 "Question 4b",
                 "Question 4c",
                 "Question 4d"]
        
        coupled = [(performance[i], names[i]) for i in range(len(performance))]
        
        total_score = sum(performance)
        
        print(f"Your score for this assignment: {total_score}/10")
        if total_score == 10:
            print("Excellent work!")
        else:
            print("The following questions have yet to be answered correctly:")
            for i in coupled:
                if i[0] == False:
                    print(i[1])
            print("Please feel free to reach out if you need any help.")
