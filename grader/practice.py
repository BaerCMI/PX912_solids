# Program name: practice.py
# Author: Aravinthen Rajkumar
# Description: proposes practice problems and solutions based on a tag.
#              all problems and answers are provided in a text file.
import os
import random 

class Problems:
    def __init__(self, ):
        print("Practice Questions Loaded!")
        # these should be loaded with 
        self.bank = "questions.txt"
        self.selected_question = None
        self.selected_solution = None
        
        self.variables = {}
        
        self.tries = 0
        
    def question(self, qtype=None):
        # qtype: the tutorial question that the practice problem is based around.
        # if qtype=None, all questions in the question bank are 
        with open("grader/questions.txt", 'r') as f:
            all_questions = f.readlines()
            
            valid_questions = []
            for i in all_questions:
                q = i.split(",")
                if qtype != None:
                    if q[0] == qtype:
                        valid_questions.append(q)
                else:
                    valid_questions.append(q)
                    
        pick_q = random.choice(valid_questions)
        
        # collect variables
        for word in pick_q[1].split(" "):
            if len(word)!=0:
                if word[0] == '$':
                    varname = word[1]+word[2]
                    varrange = word.split('[', 1)[1].split(']')[0].split('-')
                    
                    varval = 0
                    if word[1] == 'I':
                        varval = random.randint(int(varrange[0]), int(varrange[1]))
                    elif word[1] == 'F':
                        varval = round(random.uniform(int(varrange[0]), int(varrange[1])),3)
                    else:
                        raise TypeError
                
                    self.variables[varname] = varval
        
        # build question
        question_string = ""
        for word in pick_q[1].split(" "):
            if len(word)!=0:
                if word[0] == '$':
                    question_string+= str(self.variables[ word[1]+word[2] ]) + " "
                else:
                    question_string+= word + " "
        
        self.selected_question = question_string
        
        # build and evaluate solution
        evaluate = False
        for i in pick_q[2]:
            if i == '$':
                evaluate = True
        
        if evaluate == True:
            solution_string = ""
            for word in pick_q[2].split(" "):
                if len(word)!=0:
                    if word[0] == '$':
                        solution_string += str(self.variables[word[1]+word[2]])
                    else:
                        solution_string += word
            self.selected_solution = str(eval(solution_string))
        else:
            self.selected_solution = pick_q[2]
        
        
        
        print(self.selected_question)
        
    def answer(self, answ):
        if str(answ) == self.selected_solution:
            print("Correct!")
        else:
            print("Incorrect.")
            self.tries += 1
            
            if self.tries > 3:
                print("Solution: " + self.selected_solution)
                self.tries = 0