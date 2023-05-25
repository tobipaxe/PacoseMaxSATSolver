#!/bin/env python3

# Converts MaxSAT Formula from new WDIMACS format to old WDIMACS format
# usage: WDIMACSConversionNew2Old.py OldFormatFile OutputFile

import sys

file = open(sys.argv[1], 'r')

numberOfClauses = 0 
numberOfSoftClauses = 0
sumOfWeights = 0
maxValue = 0

print("Read file: ", sys.argv[1])

for x in file :
#    print("line: ", x)
    if (x == '\n') : 
        continue
    y = x.split()
    if (y[0] == '') or (y[0] == 'c') :
        continue
    elif (y[0].isnumeric()) == True :
        weight = int(y[0])
        sumOfWeights += weight
        numberOfSoftClauses += 1
    elif (y[0] != 'h') :
        print("ERROR - strange first value")
        exit(0)
    numberOfClauses += 1
    for z in y[1:] :
#        print(abs(int(z)), " ", z, "   ")
#        print(maxValue)
        if (abs(int(z)) > maxValue) :
            maxValue = abs(int(z))

print("NoClauses: ", numberOfClauses)
print("HighestVarIndex: ", maxValue)
print("SumOfSoftWeights: ", sumOfWeights)
print("NoSoftClauses: ", numberOfSoftClauses)
top = sumOfWeights + 1

file.close()

readFile = open(sys.argv[1], 'r')
writeFile = open(sys.argv[2], 'w')

writeFile.write("p wcnf " + str(maxValue) + " " + str(numberOfClauses) + " " + str(top) + '\n')

for x in readFile :
#    print("line: ", x)
    if (x == '\n') : 
        continue
    y = x.split()
    if (y[0] == '') or (y[0] == 'c') or ((y[0].isnumeric()) == True):
        writeFile.write(x)
        continue
    elif (y[0] == 'h') :
        writeFile.write(str(top) + " " + str(" ".join(y[1:])) + '\n')
    else :
        print("ERROR - strange first value2")
        exit(0)
readFile.close()
writeFile.close()

print("DONE, file written: ", sys.argv[2])
