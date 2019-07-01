# -*- coding: utf-8 -*-
"""
Created on Fri May 24 12:46:50 2019

@author: grans
"""
#Constants in the program
import subprocess
import random

a = "0"

b = "1"

alpha = "1/2"

beta = "1/3"

n = 2

#function in format bertini will understand
p = "2 * y^3"

def writeD1():#v stands for variables, funcs stands for the amount of functions

    inputFile = open("d1Python.txt", "w")
    
    inputFile.write("INPUT\n\n")
    
    inputFile.write("variable_group y;\n\n")
    
    inputFile.write("function f;\n\n")
    
    inputFile.write("constant alpha, beta, a, b, n;\n\n")
    
    inputFile.write("alpha = " + alpha + ";\n\n")
    
    inputFile.write("beta = " + beta + ";\n\n")
    
    inputFile.write("a = " + a + ";\n\n")
    
    inputFile.write("b = " + b + ";\n\n")
    
    inputFile.write("n = 1;\n\n")
    
    inputFile.write("h = ((b-a)/(n+1));\n\n")
    
    inputFile.write("f = alpha - 2 * y + beta - h^2 * (" + p + ");\n\n")
    
    inputFile.write("END;")

    inputFile.close()

def writeSS(yN, nValue): #yN is the answer that we know in the equation

    startSystemFile = open("startSystem.txt", "w")

    startSystemFile.write("INPUT\n\n")

    startSystemFile.write("variable_group y;\n\n")

    startSystemFile.write("function f;\n\n")

    startSystemFile.write("constant alpha, beta, a, b, yN, n;\n\n")

    startSystemFile.write("alpha = " + alpha + ";\n\n")

    startSystemFile.write("beta = " + beta + ";\n\n")
   
    startSystemFile.write("a = " + a + ";\n\n")
   
    startSystemFile.write("b = " + b + ";\n\n")
   
    startSystemFile.write("n = " + str(nValue) + ";\n\n")
    
    startSystemFile.write("yN = " + yN + ";\n\n")
   
    startSystemFile.write("h = ((b-a)/(n+1));\n\n")
   
    startSystemFile.write("f = " + yN + " - 2 * y + beta - h^2 * (" + p + ");\n\n")

    startSystemFile.close()

def makeHomotopyFile(nValue):

    hFile = open("homotopy.txt", "w")
    
    hFile.write("CONFIG\n\n")

    hFile.write("USERHOMOTOPY: 1;\n\n")

    hFile.write("END;\n\n")

    hFile.write("INPUT\n\n")

    hFile.write("variable " + makeVariableString(nValue) + "\n")

    hFile.write("function " + makeFunctionString(nValue) + "\n")

    hFile.write("pathvariable t;\n\n")

    hFile.write("parameter s;\n\n")

    hFile.write("constant a, b, alpha, beta, n;\n\n")

    hFile.write("random gamma;\n\n")

    hFile.write("alpha = " + alpha + ";\n\n")
                                                     
    hFile.write("beta = " + beta + ";\n\n")
   
    hFile.write("a = " + a + ";\n\n")
   
    hFile.write("b = " + b + ";\n\n")
   
    hFile.write("n = " + str(nValue) + ";\n\n")

    hFile.write("s = t;\n\n")

    hFile.write("gamma_t = (gamma ^ 2) * s + (1 -s);\n\n")

    hFile.write("h_t =  gamma * s * ( ( b - a ) / (n + 1) ) + ( 1 - s ) * ( ( b - a ) / (n + 2) );\n\n")
    
    hFile.write("y_1 = ( 1 - s ) * y" + str(nValue + 1) + " + (gamma ^ 2) * beta * s;\n\n")

    hFile.write("y_2 = beta * ( 1 - s );\n\n")

    hFile.write("y0 = alpha;\n\n")
    
    writeEquations(nValue, hFile)

    hFile.write("END;")

    hFile.close()


def writeStart(finalSol):
    
    startFile = open("start.txt", "w")

    startFile.write(str(len(finalSol) / 2)+ "\n\n")

    for i in finalSol:
        startFile.write(i);

    startFile.close()


def runBertini(txtFile, startFile = "none"): #will run Bertini with the file as input
    if (startFile == "none"):
        bertiniPath = r'./bertini'
        subprocess.call([bertiniPath, txtFile], stdout=subprocess.PIPE)
    else: 
        bertiniPath = r'./bertini'
        subprocess.call([bertiniPath, txtFile, startFile])

def readInSolutionsFinite(listOfSolutions, nValue):
    inputFile = open("finite_solutions", "r")
    if inputFile.mode == 'r':
        lines = inputFile.readlines()
    for j in lines:
        if j == "\n":
            lines.remove(j)
    for k in range(1, len(lines)):
        listOfSolutions.append(lines[k])

def readInSolutionsNon(listOfSolutions, nValue):
    inputFile = open("nonsingular_solutions", "r")
    if inputFile.mode == 'r':
        lines = inputFile.readlines()
    for j in lines:
        if j == "\n":
            lines.remove(j)
    for k in range(1, len(lines)):
        listOfSolutions.append(lines[k])
        
def combineSolutions(yN, ssSol, finalSol, nValue):
    if(nValue == 1):
        for j in ssSol:
            finalSol.append(yN + j)
            finalSol.append("\n")
    else:
        for k in ssSol:
            solutionString = ""
            beginningIndex = dNSol.index(yN) - (nValue - 1)
            while(beginningIndex <= dNSol.index(yN)):
                solutionString += dNSol[beginningIndex]
                beginningIndex += 1
            solutionString += k
            finalSol.append(solutionString)
            finalSol.append("\n")
    return finalSol
  

def fixSolutions(listOfSolutions):
    fixedSolutions = []
    for i in range(0, len(listOfSolutions)):
        listOfSolutions[i].replace("\n", '')
        splitSolutions = listOfSolutions[i].split()
        newSolution = splitSolutions[0] + " + (" + splitSolutions[1] + ") * I \n"
        fixedSolutions.append(newSolution)
    return fixedSolutions

def makeVariableString(nValue):
    variableString = ""
    for i in range(1, nValue + 1):
        variableString += "y" + str(i) + ", "
    variableString += "y" + str(nValue + 1) + ";\n"
    return variableString

def makeFunctionString(nValue):
    functionString = ""
    for i in range(0, nValue):
        functionString += "f" + str(i) + ", "
    functionString += "f" + str(nValue) + ";\n"
    return functionString

def writeEquations(nValue, fileObject):
    tempValue = nValue
    while(tempValue > 1):
        newEquation = p.replace("y", "y" + str(tempValue - 1))
        fileObject.write("f" + str(tempValue - 2) + " = (gamma_ t) * (y" + str(tempValue -2) + " - 2 * y" + str(tempValue - 1) + " + y" + str(tempValue) + ") - (h_t) ^ 2 *(" + newEquation + ");\n")
        tempValue = tempValue - 1
    Y_1Equation = p.replace("y", "y" + str(nValue))
    fileObject.write("f" + str(nValue -1) + " = (gamma_t) * (y" + str(nValue - 1) + " - 2 * y" + str(nValue) + ") + y_1 - (h_t) ^ 2 *(" + Y_1Equation + ");\n")
    Y_2Equation = p.replace("y", "y" + str(nValue + 1))
    fileObject.write("f" + str(nValue) + " = (gamma_t) * (y" + str(nValue) + " - 2 * y" + str(nValue + 1) + " + beta) - (h_t) ^ 2 * (" + Y_2Equation + ");\n\n") 

"""
main

"""

n = 3

dNSol = []
ssSol=[]
finalSol = []
writeD1()
runBertini("d1Python.txt")
readInSolutionsFinite(dNSol, 1)
fixedSolutionsDN = fixSolutions(dNSol)
for j in range(0, len(fixedSolutionsDN)):
    writeSS(fixedSolutionsDN[j], 1)
    runBertini("startSystem.txt")
    readInSolutionsFinite(ssSol, 1)
    finalSol = combineSolutions(dNSol[j], ssSol, finalSol, 1)
    ssSol = []
writeStart(finalSol)
makeHomotopyFile(1)
runBertini("homotopy.txt", "start.txt")
print("This is iteration 1")

#When n is more than 1"

for nValue in range(2, n+1):

    dNSol = []
    ssSol = []
    finalSol = []
    readInSolutionsNon(dNSol, nValue)
    fixedSolutionsDN = fixSolutions(dNSol)
    for k in range(nValue -1, len(fixedSolutionsDN), nValue):
        writeSS(fixedSolutionsDN[k], nValue)
        runBertini("startSystem.txt")
        readInSolutionsFinite(ssSol, nValue)
        finalSol = combineSolutions(dNSol[k], ssSol, finalSol, nValue)
        ssSol = []
    writeStart(finalSol)
    makeHomotopyFile(nValue)
    runBertini("homotopy.txt", "start.txt")
    print("Iteration: ", nValue)
    
    #subprocess.call(['rm', 'start.txt', 'startSystem.txt', 'homotopy.txt'])
