import sympy as sym
import numpy as np
from sympy.matrices import Matrix, eye, zeros, ones, diag, MatMul

class Workshop1:
    print("Library Loaded!")
    def __init__(self, ):
        self.q1a = False
        self.q1b = False
        self.q1c = False

        self.q2a = False
        self.q2b = False
        self.q2c = False
        
        self.q3a = False
        self.q3b = False
        
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

    def hint1c(self,):
        print("Remember to represent the exponential using SymPy. Otherwise it")
        print("should be the same as the two previous questions.")        
    def check1c(self,F_c, homogeneous):
        X_1, X_2, X_3 = sym.symbols('X_1 X_2 X_3')
        if F_c == Matrix([[sym.exp(X_1), 0, 0], [0, 0, -1], [0, 1, 0]]):
            if homogeneous==False:
                print("\033[1;32m Correct!")
                self.q1c = True
            else:
                print("\033[0;31m Incorrect. Recall the definition of homogeneity.")
                self.q1c = False  
        else:
            print("\033[0;31m Incorrect. Check your answer for the deformation gradient.")
            self.q1c = False             

    
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
        print("Another way to go about it is to apply the deformation gradient,")
        print("but the change of the origin during the deformation also needs")
        print("to be taken into account.")

    def check2b(self, deformed_corners):
        if deformed_corners == [[4.0, 2.0],
                                [7.0, 5.0],
                                [5.0, 5.0],
                                [2.0, 2.0]]:
            print("\033[1;32m Correct!")
            self.q2b = True
        elif deformed_corners == sym.Matrix([[4.0, 2.0],
                                             [7.0, 5.0],
                                             [5.0, 5.0],
                                             [2.0, 2.0]]):
            print("\033[1;32m Correct!")
            self.q2b = True
        else:
            print("\033[0;31m Incorrect.")
            print("Note: if you've modified the order of the square_corners list, you may be getting an error due to the order of your points.")
            self.q2b = False                 

    
    def hint2c(self,):
        print("This question requires linear algebra. NumPy is your friend here.")
        print("Array multiplication can be performed with np.dot(A, x), where A ")
        print("is a matrix and x is a vector.")

        
    def check2c(self, e1, e2, E1, E2):
        ce1, ce2, cE1, cE2 = False, False, False, False
        if e1[0]==1 and e1[1]==0:
            ce1=True
        else:
            ce1=False
            print("\033[0;31m def_e1 is incorrect.")
            self.q2b = False 
           
        if e2[0]==-1 and e2[1]==-2/3:
            ce2=True
        else:
            ce2=False
            print("\033[0;31m def_e2 is incorrect.")
            self.q2b = False 
                    
        if E1[0]==1 and E1[1]==0:
            cE1=True
        else:
            cE1=False   
            print("\033[0;31m def_E1 is incorrect.")
            self.q2b = False 
            
        if E2[0]==-1.5 and E2[1]==-1.5:     
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
        alpha, L = sym.symbols('alpha, L')
        if ifts == Matrix([[0, 0.5*alpha/L, 0], 
                           [0.5*alpha/L, 0, 0],
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
            

    def hint3c(self):
        print("This should be very similar to a previous problem.")
        print("Simply translate the points that you are given.")
    # RESULTS    ------------------------------------------------------
    
    def results(self):
        performance = [self.q1a,
                       self.q1b, 
                       self.q1c,  
                       self.q2a,
                       self.q2b,
                       self.q2c, 
                       self.q3a, 
                       self.q3b]
        
        names = ["Question 1a",
                 "Question 1b",
                 "Question 1c",
                 "Question 2a",
                 "Question 2b",
                 "Question 2c",
                 "Question 3a",
                 "Question 3b"]
        
        coupled = [(performance[i], names[i]) for i in range(8)]
        
        total_score = sum(performance)
        
        print(f"Your score for this assignment: {total_score}/8")
        if total_score == 8:
            print("Excellent work!")
        else:
            print("The following questions have yet to be answered correctly:")
            for i in coupled:
                if i[0] == False:
                    print(i[1])
            print("Please feel free to reach out if you need any help.")
