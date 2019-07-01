# -*- coding: utf-8 -*-
"""
Created on Fri May 24 12:46:50 2019

@author: grans
"""
#Constants in the program
import subprocess
import random
import cmath
import os.path

a = "0"

b = "1"

alpha = "0"

beta = "0"

#function in format bertini will understand
p = "-.5 * (y - ((y^3)/6) + ((y^5)/120))"

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

    hFile.write("TRACKTOLBEFOREEG: 1e-12;\n\n")

    hFile.write("DEGREEBOUND: 2;\n\n")

    hFile.write("COEFFBOUND: 2.5;\n\n")

    hFile.write("END;\n\n")

    hFile.write("INPUT\n\n")

    hFile.write("variable " + makeVariableString(nValue) + "\n")

    hFile.write("function " + makeFunctionString(nValue) + "\n")

    hFile.write("pathvariable t;\n\n")

    hFile.write("parameter s;\n\n")

    hFile.write("constant gamma, a, b, alpha, beta, n;\n\n")
    
    hFile.write("gamma = .5 + .5 * I;\n\n") 

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

def readInSolutions(listOfSolutions, numSol):
    if(os.path.exists("main_data")):
        inputFile = open("main_data", "r")
        if inputFile.mode == 'r':
            temp = "null"
            while(temp != ""):
                temp = inputFile.readline()
                if(temp.find("Cycle number") != -1):
                    for i in range(0, numSol):
                        temp = inputFile.readline()
                        listOfSolutions.append(temp)

def readInSolutionsFinite(listOfSolutions, nValue):
    inputFile = open("finite_solutions", "r")
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

def filterSolutions(listOfSolutions, nValue):
    newList = []
    for i in range(0, len(listOfSolutions), nValue):
        lastIndex = i + nValue - 1
        if (abs(calculatingMag(listOfSolutions[i]) - calculatingMag(listOfSolutions[lastIndex])) < 10 ** -16):
            for j in range(i, lastIndex + 1):
                newList.append(listOfSolutions[j])
    return newList

def calculatingMag(string):
    string.replace("\n", "")
    splitSolutions = string.split()
    real = float(splitSolutions[0])
    imaginary = float(splitSolutions[1])
    return cmath.sqrt((real ** 2) + (imaginary ** 2))

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
derivative filter methods
"""
def derivativeFiltering(listOfStrings, listOfComplex,nValue):
    sumOfComparison = 0
    newList = []
    pDerivative = lambda x: -.5 * (1 - (.5 * x ** 2) + (.05 * x ** 5))
    for j in range(0, len(listOfComplex), nValue):
        temp = list(listOfComplex[j:j+nValue])
        for i in range(2, nValue - 2):
            complexNumber = thirdDerivative(temp[i-2], temp[i-1], temp[i+1], temp[i+2]) - (pDerivative(temp[i]) * firstDerivative(temp[i+1], temp[i-1]))
            instance = cmath.sqrt(complexNumber.real ** 2 + complexNumber.imag ** 2)
            sumOfComparison += instance.real
        if(sumOfComparison < 10 ** -16):
            stringsOfSolutions = list(listOfStrings[j:j+nValue])
            newList.extend(stringsOfSolutions)
    return listOfStrings

def convertToComplex(listOfStrings):
    newList = []
    for i in listOfStrings:
        splitString = i.split()
        real = float(splitString[0])
        imaginary = float(splitString[1])
        complexNumber = complex(real, imaginary)
        newList.append(complexNumber)
    return newList 

def thirdDerivative(yplus2, yplus1, yminus1, yminus2):
    numerator = (yplus2 - (2 * yplus1) + (2 * yminus1) - yminus2)
    denominator = (2 * ((float(b) - float(a)) / (nValue + 1)) ** 3)
    return numerator/denominator

def firstDerivative(yplus1, yminus1):
    numerator = (yplus1 - yminus1)
    denominator = (2 * ((float(b) - float(a)) / (nValue + 1)))
    return numerator/denominator

"""
end of set of methods
"""
"""
def saveFiles(nValue):
    subprocess.call(["mkdir", "outputFiles" + str(nValue])
    subprocess.call(["cp", "homotopy.txt", "outputFiles" + str(nValue) + "/homotopy" + str(1) + ".txt"])
    subprocess.call(["cp", "start.txt", "outputFiles" + str(nValue) + "/start" + str(1) + ".txt"])
    subprocess.call(['rm', 'start.txt', 'startSystem.txt', 'homotopy.txt'])
"""
"""
main
"""

n = 2

dNSol = []
ssSol=[]
finalSol = []
writeD1()
runBertini("d1Python.txt")
readInSolutions(dNSol, 1)
fixedSolutionsDN = fixSolutions(dNSol)
for j in range(0, len(fixedSolutionsDN)):
    writeSS(fixedSolutionsDN[j], 1)
    runBertini("startSystem.txt")
    readInSolutions(ssSol, 1)
    finalSol = combineSolutions(dNSol[j], ssSol, finalSol, 1)
    ssSol = []
writeStart(finalSol)
makeHomotopyFile(1)
runBertini("homotopy.txt", "start.txt")
subprocess.call(["cp", "main_data", "main_data_1.txt"])
print("This is iteration 1")
"""
subprocess.call(["mkdir", "outputFiles1"])
subprocess.call(["cp", "homotopy.txt", "outputFiles1/homotopy" + str(1) + ".txt"])
subprocess.call(["cp", "start.txt", "outputFiles1/start" + str(1) + ".txt"])
subprocess.call(['rm', 'start.txt', 'startSystem.txt', 'homotopy.txt'])
"""
#When n is more than 1"

for nValue in range(2, n+1):

    dNSol = []
    ssSol = []
    finalSol = []
    readInSolutions(dNSol, nValue)
    if(nValue > 4):
        dNSol = filterSolutions(dNSol, nValue)
        tempdNSol = convertToComplex(dNSol)
        dNSol = derivativeFiltering(dNSol, tempdNSol, nValue)
    fixedSolutionsDN = fixSolutions(dNSol)
    for k in range(nValue -1, len(fixedSolutionsDN), nValue):
        writeSS(fixedSolutionsDN[k], nValue)
        runBertini("startSystem.txt")
        readInSolutions(ssSol, 1)
        finalSol = combineSolutions(dNSol[k], ssSol, finalSol, nValue)
        ssSol = []
    writeStart(finalSol)
    makeHomotopyFile(nValue)
    runBertini("homotopy.txt", "start.txt")
    print("Iteration: ", nValue)
    subprocess.call(["cp", "homotopy.txt", "homotopy" + str(nValue) + ".txt"])
    subprocess.call(["cp", "start.txt", "start" + str(nValue) + ".txt"])
    subprocess.call(['rm', 'start.txt', 'startSystem.txt', 'homotopy.txt'])
