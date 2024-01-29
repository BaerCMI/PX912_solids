import sympy as sym
import numpy as np
import math
from sympy.matrices import Matrix, eye, zeros, ones, diag, MatMul

class Workshop4:
    def __init__(self, ):
        print("Grader Library Loaded!")

        self.q1a = False
        self.q1b = False

        self.q2b  = False
        self.q2c = False
        self.q3a = False
        self.q4 = False
        self.q5 = False
        self.q6 = False
        self.q7 = False

    # Problem 1 ----------------------------------------------------------------------------------
    def hint1a(self,):
        print("How many Gauss points do you need in this case?")
        print("Normally, weights are coded within FEA software. For now, use the weight values")
        print("that you can find in the lecture notes.")
        
    def check1a(self, numerical, exact):
        if np.isclose(numerical,26/15):
            if np.isclose(numerical,26/15):
                print("\033[1;32m Correct!")
                self.q1a = True
            else:
                print(f"\033[0;31m Incorrect. Check your exact calculation.")
                self.q1a = False
        else:
            print(f"\033[0;31m Incorrect. Check your numerical calculation.")
            self.q1a = False

    def hint1b(self,):
        print("Here you also need to consider the conversion from physical to natural")
        print("coordinates and multiply the integral by the Jacobian. The equations ")
        print("for this can be found in the lecture notes.")           
    def check1b(self, numerical, exact):
        if np.isclose(numerical,76/3):
            if np.isclose(exact,76/3):
                print("\033[1;32m Correct!")
                self.q1b = True
            else:
                print(f"\033[0;31m Incorrect. Check your exact calculation.")
                self.q1b = False
        else:
            print(f"\033[0;31m Incorrect. Check your numerical calculation.")
            self.q1b = False
    
    # Problem 2 ----------------------------------------------------------------------------------
    def hint2ai(self,):
        print("The 4 shape functions have the same form, but differ by the signs of xi and eta")
        print("depending on the position of the nodes.")        

    def hint2aii(self,):
        print("Fill a 4x2 matrix with the physical coordinates of the nodes, which")
        print("can be found in the picture above.")

    def hint2aiii(self,):
        print("Write the derivatives of the shape functions in terms of eta and xi.")
        print("Use the M.diff(x) function, with M being the matrix and x the variable that is being")
        print("differentiated with resepect to.")

    def hint2aiv(self,):
        print("Determine the inverse and determinant using M.inv() and M.det()")
    
    def hint2av(self,):
        print("A general form of the strain-displacement matric for a four-node 2D element")
        print("is in the lecture notes. The final matrix should have dimensions (3x8).")
            
    def hint2b(self,):
        print("You can use the plane stress function you defined in Workshop 2")  
    
    def check2b(self, C):
        if isinstance(C, Matrix):
            C = np.array(C).astype(np.float64)
        if C.shape == (3,3):
            if np.array_equal(np.round(C, 6), np.round(np.array([[32967032.96703297, 
            9890109.89010989, 0.],
            [ 9890109.89010989, 32967032.96703297, 0.],
            [0., 0., 11538461.53846154]]), 6)):
                print("\033[1;32m Correct!")
                self.q2b = True
            else:
                print(f"\033[0,31m Incorrect. Check the values of your plane stress matrix")
                self.q2b = False
        else:
            print(f"\033[0;31m Incorrect. Check the shape of your plane stress matrix")
            self.q2b = False


    def hint2c(self,):
        print("The format when evaluating a function M at the integration points is:")
        print("M.evalf(subs = {variable1:value1, variable2:value2})")
    
    def check2c(self, K_e):
        if isinstance(K_e, Matrix):
            K_e = np.array(K_e).astype(np.float64)
        
        if K_e.shape == (8,8):
            if np.array_equal(np.round(K_e,6), np.round(np.array([[ 14898562.97548605,  
            -7417582.41758242,  -6656804.73372781,
            1648351.64835165,  -9763313.60946745,   6593406.59340659,
            1521555.36770921,   -824175.82417582],
            [ -7417582.41758242,  27467244.29416737,   2472527.47252747,
            -24582628.90955198,   6593406.59340659, -16768808.11496196,
            -1648351.64835165,  13884192.73034657],
            [ -6656804.73372781,   2472527.47252747,  10777683.85460693,
            3296703.2967033 ,   1521555.36770921,  -1648351.64835165,
            -5642434.48858833,  -4120879.12087912],
            [  1648351.64835165, -24582628.90955198,   3296703.2967033 ,
            26024936.60185968,   -824175.82417582,  13884192.73034657,
            -4120879.12087912, -15326500.42265427],
            [ -9763313.60946745,   6593406.59340659,   1521555.36770921,
            -824175.82417582,  20033812.34150464,  -8241758.24175824,
            -11792054.09974641,   2472527.47252747],
            [  6593406.59340659, -16768808.11496196,  -1648351.64835165,
            13884192.73034657,  -8241758.24175824,  38165680.47337277,
            3296703.2967033 , -35281065.08875739],
            [  1521555.36770921,  -1648351.64835165,  -5642434.48858833,
            -4120879.12087912, -11792054.09974641,   3296703.2967033 ,
            15912933.22062553,   2472527.47252747],
            [  -824175.82417582,  13884192.73034657,  -4120879.12087912,
            -15326500.42265427,   2472527.47252747, -35281065.08875739,
            2472527.47252747,  36723372.78106508]]),6)):
                print("\033[1;32m Correct!")
                self.q2c = True
            else:
                print(f"\033[0;31m Incorrect. Check the values of your elemental stiffness matrix")
                self.q2c = False
        else:
            print(f"\033[0;31m Incorrect. Check the shape of your elemental stiffness matrix")
            self.q2c = False
    
    def hint3a(self,):
        print("All the degrees of freedom that don't experience any external force should be set")
        print("to zero. The non-zero degrees of freedoms contain the respective shape function.")
        print("The traction vector only needs to to have length 2. Consider in which direction the")
        print("force is acting, the other direction will have 0 traction.")
    
    def check3a(self, f_e):
        if isinstance(f_e, Matrix):
            f_e = np.array(f_e).astype(np.float64)
        if f_e.shape == (8,1):
            if np.array_equal(f_e, np.array([[0], [-20.],[0], [0], 
            [0], [0], [0], [-20.0]])):
                print("\033[1;32m Correct!")
                self.q3a = True
            else:
                print(f"\033[0;31m Incorrect. Check the values of your nodal force vector")
                self.q3a = False
        else:
            print(f"\033[0;31m Incorrect. Check the shape of your nodal force vector")
            self.q3a = False
        
    def hint3b(self,):
        print("Consider the individual degrees of freedom and be sure to add zero")
        print("for fixed degrees of freedom.")

    def hint4(self,):
        print("Only take the entries in the matrix/vector that respond to free nodes.")
        print("The reduced elemental stiffness matrix should be a [4x4] matrix and the")
        print("reduced nodal force vector should only have length 4.")
    
    def check4(self, K_eR, f_eR):
        if isinstance(K_eR, Matrix):
            K_eR = np.array(K_eR).astype(np.float64)
        f_eR = np.array(f_eR).astype(np.float64)
        
        if K_eR.shape == (4,4):
            if f_eR.shape == (4,):
                if np.array_equal(np.round(K_eR,6), np.round(np.array([[ 20033812.34150464,  
                -8241758.24175824, -11792054.09974641, 2472527.47252747],
                [ -8241758.24175824,  38165680.47337277,   3296703.2967033 ,
                -35281065.08875739],
                [-11792054.09974641,   3296703.2967033 ,  15912933.22062553,
                2472527.47252747],
                [  2472527.47252747, -35281065.08875739,   2472527.47252747,
                36723372.78106508]]),6)):
                    if np.array_equal(f_eR, np.array([  0.,   0.,   0., -20.])):
                        print("\033[1;32m Correct!")
                        self.q4 = True
                    else: 
                        print(f"\033[0;31m Incorrect. Check the values of f_eR")
                        self.q4 = False
                else:
                    print(f"\033[0;31m Incorrect. Check the values of K_eR")
                    self.q4 = False
            else: 
                print(f"\033[0;31m Incorrect. Check the shape of f_eR")
                self.q4 = False
        else:
            print(f"\033[0;31m Incorrect. Check the shape of K_eR")
            self.q4 = False
    
    def hint5(self,):
        print("Use np.dot() or @ for matrix multiplication and np.linalg.inv() for")
        print("matrix inversion in numpy.")
        print("Alternatively, you can use np.linalg.solve to solve a linear system.")
        print("Your solution should be a vector of 4 values")
    
    def check5(self, d_eR):
        if isinstance(d_eR, Matrix):
            d_eR = np.array(d_eR).astype(np.float64)
        if d_eR.shape == (4,1):
            if np.array_equal(np.round(d_eR, 14), np.round(np.array([[-1.17777210e-06], 
            [-9.66972494e-06], [2.67425251e-06], [-9.93531521e-06]]), 14)):
                print("\033[1;32m Correct!")
                self.q5 = True
            else:
                print(f"\033[0;31m Incorrect. Check the values of d_eR")
                self.q5 = False
        else:
            print(f"\033[0;31m Incorrect. Check the shape of d_eR")
            self.q5 = False

    def hint6(self,):
        print("Evaluate the elemental stiffness matrix at the integration points you defined above")
        print("Use: M.evalf(subs = {variable1:value1, variable2:value2})")
    
    def check6(self, eps1, eps2, eps3, eps4):
        if isinstance(eps1, Matrix):
            eps1 = np.array(eps1).astype(np.float64)
        if isinstance(eps2, Matrix):
            eps2 = np.array(eps2).astype(np.float64)
        if isinstance(eps3, Matrix):
            eps3 = np.array(eps3).astype(np.float64)
        if isinstance(eps4, Matrix):
            eps4 = np.array(eps4).astype(np.float64)
        
        if np.array_equal(np.round(eps1,14), np.round(np.array([[8.82024842e-07], 
        [-6.27568703e-08], [-4.02607634e-06]]),14)):
            if np.array_equal(np.round(eps2,14), np.round(np.array([[-3.61335342e-07], 
            [-6.27568703e-08], [-3.94034886e-06]]),14)):
                if np.array_equal(np.round(eps3,14), np.round(np.array([[-1.17086818e-06], 
                [-3.45843536e-07], [2.21253041e-07]]),14)):
                    if np.array_equal(np.round(eps4,14), np.round(np.array([[6.65111171e-07], 
                    [-3.45843536e-07], [9.46655214e-08]]),14)):
                        print("\033[1;32m Correct!")
                        self.q6 = True 
                    else:
                        print(f"\033[0;31m Incorrect. Check the values for eps4")
                        self.q6 = False
                else:
                    print(f"\033[0;31m Incorrect. Check the values for eps3")
                    self.q6 = False
            else:
                print(f"\033[0;31m Incorrect. Check the values for eps2")
                self.q6 = False
        else:
            print(f"\033[0;31m Incorrect. Check the values for eps1")
            self.q6 = False

    def hint7(self,):
        print("Use the strains computed in the previous question")
        print("Evaluate the elemental stiffness matrix at the integration points you defined above")

    
    def check7(self, sig1, sig2, sig3, sig4):
        if isinstance(sig1, Matrix):
            sig1 = np.array(sig1).astype(np.float64)
        if isinstance(sig2, Matrix):
            sig2 = np.array(sig2).astype(np.float64)
        if isinstance(sig3, Matrix):
            sig3 = np.array(sig3).astype(np.float64)
        if isinstance(sig4, Matrix):
            sig4 = np.array(sig4).astype(np.float64)
        
        if np.array_equal(np.round(sig1,6), np.round(np.array([[28.4570697], [6.6544148], 
        [-46.45472703]]),6)):
            if np.array_equal(np.round(sig2,6), np.round(np.array([[-12.53282647],  
            [-5.64255405], [-45.46556381]]),6)):
                if np.array_equal(np.round(sig3,6), np.round(np.array([[-42.02048058], 
                [-22.98145024], [2.55291971]]),6)):
                    if np.array_equal(np.round(sig4,6), np.round(np.array([[18.50631133],  
                    [-4.82341267], [1.09229448]]),6)):
                        print("\033[1;32m Correct!")
                        self.q7 = True 
                    else:
                        print(f"\033[0;31m Incorrect. Check the values for sig4")
                        self.q7 = False
                else:
                    print(f"\033[0;31m Incorrect. Check the values for sig3")
                    self.q7 = False
            else:
                print(f"\033[0;31m Incorrect. Check the values for sig2")
                self.q7 = False
        else:
            print(f"\033[0;31m Incorrect. Check the values for sig1")
            self.q7 = False

    def hint2bi(self,):
        print("1. Specifiy the quadrilateral in natural coordinates.")
        print("2. Write the derivatives of the shape functions in terms of eta and xi.")
        print("   This is expression 2.50 in the notes.")
        print("3. Calculate the element Jacobian matrix. You'll need np.dot() and np.linalg.det.")
        print("4. The final matrix should have dimensions (3x8). Check against the notes!")
        
        
    def hint2biii(self,):
        print("Steps:")
        print("1. Build a boundary conditions vector, BC. You only need this to modify the")
        print("   elemental stiffness matrix.")
        print("2. Be sure to set the fixed nodes to zero.")
        print("3. Obtain the indices of the non-zero nodes.")
        print("4. Build the elemental stiffness matrix: easy because there's only one.")
        print("   Remember, you only need the component of ke that isn't fixed.")
        print("   Keep in mind that this is an integration. You will have to implement")
        print("   Gauss quadrature in order to carry this out.")        

    def hint2biv(self, ):
        print("Which nodes actually have a force acting on them? Your answer should be multiplicable with the expression for the elemental stiffness matrix.")        
    
    def hint2ci(self,):
        print("Use np.linalg.solve to solve a linear system.")
        print("Your solution should be a vector of 4 values. Remember to add the displacements")
        print("for nodes that were fixed! (basically, add a bunch of zeros).")            

    def hint2cii(self,):
        print("All required components are present. Use a for-loop to calculate each strain.")
    
    def hint2ciii(self,):
        print("Remember to use the Jacobian when summing the stresses!")

            
    def results(self):
        performance = [self.q1a,
                       self.q1b,
                       self.q2b,
                       self.q2c,
                       self.q3a,
                       self.q4,
                       self.q5,
                       self.q6,
                       self.q7]
        
        names = ["Question 1a",
                 "Question 1b",
                 "Qestion 2b",
                 "Question 2c",
                 "Question 3a",
                 "Question 4",
                 "Question 5",
                 "Question 6",
                 "Question 7"]

        coupled = [(performance[i], names[i]) for i in range(len(performance))]
        
        total_score = sum(performance)
        
        print(f"Your score for this assignment: {total_score}/{len(performance)}")
        if total_score == len(performance):
            print("Excellent work!")
        else:
            print("The following questions have yet to be answered correctly:")
            for i in coupled:
                if i[0] == False:
                    print(i[1])
            print("Please feel free to reach out if you need any help.")
