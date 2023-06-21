# Import all packages
from sympy import symbols, Piecewise, exp, Matrix, core, diff
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import scipy.integrate as spi
from scipy.optimize import *
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt
import time
import contextlib
import os
import sys
from sympy.utilities.autowrap import autowrap
import csv

import fcntl
import numpy.matlib

# Define functions
def replaceS(stressName,valueName):
    sText = stressName + " = " + valueName
    return sText

def makeReplaceSText(stressName, attach, iRange):
    sTextList = []
    for i in iRange:
        sText = attach + str(i) + " = " + attach + str(i) + ".subs(" + stressName + ",Stress)"
        sTextList = sTextList + [sText]
    return sTextList

def makeExposurePWText(doseNames, doseList, funcType = "stepwise"):
    if funcType == "stepwise":
        # Make the piecewise text
        sText = "exposurePiecewise = Piecewise("
        for i in range(0,len(doseList)):
            sText = sText + "(" + doseNames[i] + ", (dose >= " + str(doseList[i] - 0.00001) + ") & (dose <= " + str(doseList[i] + 0.00001) + ")), "
        sText = sText + "(0, True))"
    elif funcType == "expDecay":
        # Make the piecewise text
        sText = "exposurePiecewise = Piecewise("
        for i in range(0, len(doseList)):
            sText = sText + "(" + doseNames[i] + " * exp(-tsym * tau1)" + ", (dose >= " + str(doseList[i] - 0.00001) + ") & (dose <= " + str(doseList[i] + 0.00001) + ")), "
        sText = sText + "(0, True))"
    return sText

def makeutText(PWfunctionName, tName):
    sText = "Stress = " + PWfunctionName + " * " + "Piecewise((1, (" + tName + " > 0)), (0,True))"
    return sText

def makeutMatrixText(name, attach, iRange, add1 = False):
    # Make a string called 'sText' that can make a matrix when executed with exec(sText).
    # attach is the name of the mathematical function, e.g. f or g
    # iRange is the range passed to the for loop. Usually goes from 0 to NumberOfStateVariables
    # add1 is True or False, depending on whether 1 should be added to i
    # Returns sText
    sText = name + " = Matrix([" # + attach1 + ","
    if add1:
        for i in iRange:
            if i == iRange[-1]:
                sText = sText + attach + str(i + 1) + "])"
            else:
                sText = sText + attach + str(i + 1) + ","
    else:
        for i in iRange:
            if i == iRange[-1]:
                sText = sText + attach + str(i) + "])"
            else:
                sText = sText + attach + str(i) + ","
    return(sText)

def allUnique(x):
    # Check whether all values in list x are unique. Returns False or True
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)

def makeSymbols(symbolList):
    # Make a string called 'sText' that can initialize sympy symbols when executed with exec(sText).
    # symbolList is a list of strings that should be initialized
    # Returns sText
    sText = ""
    sTextspace = ""
    if not allUnique(symbolList):
        symbolList = list(set(symbolList))

    for symbol in symbolList:
        sTextspace = sTextspace + symbol + ' '
        if (symbol == symbolList[-1]):
            sText = sText + symbol
        else:
            sText = sText + symbol + ','
    sText = sText + ' = symbols(\''
    sText = sText + sTextspace + '\')'
    return sText

def makeMatrixText(attach, iRange, add1 = False):
    # Make a string called 'sText' that can make a matrix when executed with exec(sText).
    # attach is the name of the mathematical function, e.g. f or g
    # iRange is the range passed to the for loop. Usually goes from 0 to NumberOfStateVariables
    # add1 is True or False, depending on whether 1 should be added to i
    # Returns sText

    # Test whether there is input
    if len(iRange) == 0:
        print("No parameters given. Building empty matrix.")
        sText = attach + " = Matrix([])"
    else:
        sText = attach +  " = Matrix(["
        if add1:
            for i in iRange:
                if i == iRange[-1]:
                    sText = sText + attach + str(i + 1) + "])"
                else:
                    sText = sText + attach + str(i + 1) + ","
        else:
            for i in iRange:
                if i == iRange[-1]:
                    sText = sText + attach + str(i) + "])"
                else:
                    sText = sText + attach + str(i) + ","
    return(sText)

def makeDerivText(derivName, attach, iRange, List):
    # Make a list of strings called 'textList' that will calculate the derivatives of a function
    # with respect to the parameters defined by List when executed with exec(sText).
    # derivName is the name of the newly made derivative function
    # attach is the name of the mathematical function for which the derivative should be calculated, e.g. f or fSnow
    # iRange is the range passed to the for loop. Usually goes from 0 to NumberOfStateVariables
    # List is the list of parameters with respect to which the derivative should be calculated
    # Returns textList
    textList = []
    for i in iRange:
        sText = attach + derivName + str(i) + " = diff(" + attach + ", " + List[i] + ")"
        textList.append(sText)
    return(textList)

def makeMatrixText2(attach, iRange):
    # Make a string called 'sText' that can join the derivatives per row to form a matrix when executed with exec(sText).
    # attach is the name of the list of derivatives, e.g. fSnow
    # iRange is the range passed to the for loop. Usually goes from 0 to NumberOfStateVariables
    # Returns sText

    # Test whether there is input
    if len(iRange) == 0:
        print("No parameters given. Building empty matrix.")
        sText = attach + "M = Matrix([])"
    else:
        sText = attach + "M = "
        for i in iRange:
            if i == 0:
                sText = sText + attach + str(i)
            else:
                sText = sText + ".row_join(" + attach + str(i) + ")"
    return(sText)

def makeAutowrapText(attach, backend, arguments):
    '''
    Make a string called 'sText' that will 'autowrap' the function when executed with exec(sText).
    :param attach: is the name of the object, e.g. f of fSnowM
    :param backend: specifies the backend programming language
    :param arguments: passes the list of parameters to 'args' of autowrap. The order of the arguments is important. \
    The same order is used for the parameters in map2Jy, to map the parameter names to their calculated values
    :return: a text called sText that will 'autowrap' the function when executed with exec(sText)
    '''
    sText = attach + "autowrap = eval(\'autowrap(" + attach + ", backend = \"" + backend + "\", args = " + str(tuple(arguments)).replace("'","") + ")')"
    return sText

def makeParamDictText(stateList, paraList):
    '''
    to initialize symbolic states and parameters
    :param stateList: list of the states
    :param paraList: list of the parameters
    :return:
    sText: the string to initialize symbolic states and parameters
    '''
    sStart = "dictXe = {"
    sEnd = "}"
    sText = ""
    for i in range(0, len(stateList)):
        if len(paraList)==0:
            sText = sText  + str(stateList[i]) + " : wx[" + str(i) + "]"
        else:
            sText = sText  + str(stateList[i]) + " : wx[" + str(i) + "],"
    for j in range(0, len(paraList)):
        if j == len(paraList)-1:
            sText = sText + str(paraList[j]) + " : p[" + str(j) + "]"
        else:
            sText = sText + str(paraList[j]) + " : p[" + str(j) + "],"
    sText = sStart + sText + sEnd
    return sText

def makeTextRunAutowrap(autowrapName, statelistName, parlistName, statelist, parlist, withTimeDose = True):
    # Make nwvarName
    sText = autowrapName + "("
    for i in range(0,len(statelist)):
        sText = sText + statelistName + "[" + str(i) + "]" + ", "
    for i in range(0,len(parlist)-1):
        sText = sText + parlistName + "[" + str(i) + "]" + ", "
    sText = sText + parlistName + "[" + str(i+1) + "]"

    if withTimeDose:
        sText = sText + ", t, dose)"
    elif not withTimeDose:
        sText = sText + ")"

    return sText

def makeTextRunutAutowrap(autowrapName, doseparlistName, doselist, pkparlistName, pkparmList, withTimeDose = True):
    # Make nwvarName
    sText = autowrapName + "("
    if len(doselist) > 1:
        for i in range(0,len(doselist)-1):
            sText = sText + doseparlistName + "[" + str(i) + "]" + ", "
        sText = sText + doseparlistName + "[" + str(i+1) + "]"
    elif len(doselist) == 1:
        sText = sText + doseparlistName + "[0]"

    if withTimeDose and len(doselist) > 0:
        sText = sText + ", t, dose"
    elif withTimeDose and len(doselist) == 0:
        sText = sText + "t, dose"

    if len(pkparmList) > 1:
        for i in range(0, len(pkparmList) - 1):
            sText = sText + ", " + pkparlistName + "[" + str(i) + "]"
        sText = sText + ", " + pkparlistName + "[" + str(i + 1) + "]" + ")"
    elif len(pkparmList) == 1:
        sText = sText + ", " + pkparlistName + "[0]" + ")"
    else:
        sText = sText + ")"

    return sText

def transformParms(transformParms, useParms, df, matrix):
    '''
    Calculate the star parameters from the (initial or) steady states and other parameters
    :param transformParms: a list of the parameters that should be calculated from others, e.g. paraStarList
    :param useParms: a list of the parameters that should be used to calculate transformParms, e.g. stateOList + paraSnowList
    :param df: a pandas data frame with the names of parameters in useParms as row names and their value, e.g. parmsDf
    :param matrix: a matrix that contains the functions that can be used to calculate the value of transformParms, e.g. fR
    :return: an extended pandas data frame with the names of all parameters in useParms and transformParms as row names and their (calculated) value
    '''
    i = 0
    for tparm in transformParms:
        # Assign the symbolic expressions that can be found in the fR-matrices to the paraStar parameters
        df.loc[tparm] = matrix[i]

        # Replace symbols by values
        for uparm in useParms:
            df.loc[tparm] = df.loc[tparm][0].subs(uparm, df.loc[uparm])

        i += 1

        # Convert the sympy Float to regular Python float
        df.loc[tparm] = float(df.loc[tparm])

    # Convert the df.dtypes to df.float64
    df = df.apply(pd.to_numeric, errors='coerce')

    return df

def assignParametersText(parList, valueListName):
    sText = "["
    for parm in parList:
        if parm == parList[-1]:
            sText = sText + parm + "]"
        else:
            sText = sText + parm + ", "
    sText = sText + " = " + valueListName
    return sText

def xOParam2Rparam(x0, dpksnowParms, Nso, paraStarList, paraSnowList, doseParmsList, PKparms, fR, allIniStatesList):
    '''
    Function to calculate the star parameters from the initial states and snow parameters
    :param xO: initial/steady states in linear scale
    :param thetaSnow: parameter Snow in linear scale
    :return: rParam, the parameters can be determined from the initial states and parameter Snow
    '''

    # Make a data frame
    parmsDf = pd.DataFrame({'value': np.append(x0, dpksnowParms)}, index=allIniStatesList + doseParmsList + PKparms + paraSnowList)

    # Calculate the values of the star parameters that are determined by the steady states
    parmsDf = transformParms(transformParms=paraStarList, useParms=allIniStatesList + doseParmsList + PKparms + paraSnowList, df=parmsDf, matrix=fR)

    # Get rParam from the df
    if len(paraStarList) > 0:
        rParam = np.array(parmsDf["value"][-len(paraStarList):])
    else:
        rParam = np.array([])

    # Check whether one of the star parameters is negative
    # if np.min(rParam)<0:
    #     print("Warning: At least one of the star parameters is negative:")
    #     print(rParam)

    return rParam

def rhs(t, wx, p, d, doseParmsList, pks, PKparms, dose, fautowrap, Ns):
    '''
    Function to create the right hand side of the ODE equations
    :param t: time
    :param wx: initial states
    :param p: parameters in linear scale
    :param dose: dose
    :param fautowrap: fautowrap
    :param Ns: number of state variables
    :return: f; the right hand side of differential equation
    '''
    # Make f
    # sText = makeTextRunAutowrap(nwvarName="f", autowrapName="fautowrap", statelistName="wx", parlistName="p", statelist=wx, parlist=p)
    sText = makeTextRunAutowrap(autowrapName="fautowrap", statelistName="wx", parlistName="p", statelist=wx, parlist=p)
    f = eval(sText)
    # print("Executed: f = %s" % sText)

    # Convert f matrix format to array format
    f = np.array(f.tolist(), dtype=float).reshape(Ns, )

    return f

def fRightHandSide(x, t, arguments):
    '''
    IMPORTANT: this function is utilized in the odeint function
    :param x: initial states in linear scale
    :param t: time
    :param arguments: List of additional arguments
    :return: fall; the right hand of the differential equations f, dotR and dotS, including the state and the sensitivity items
    IMPORTANT: this function is utilized in the pyx file for the sundial/cvode solver
    '''

    # Unravel the arguments
    [pplusdose, stateOList, doseParmsList, PKparms, Ns, Nse, Np, Nso, Nd, Npk, Nsnow, Nstar, rParam, autowrapList] = arguments
    # Unravel pplusdose, the parameters in linear scale
    [x0, doseParms, pkParms, thetaSnow, dose] = pplusdose

    # Unravel the autowrapList
    [fautowrap, gautowrap, fXMautowrap, \
     fStarMautowrap, fSnowMautowrap, fRXOMautowrap, fRSnowMautowrap, \
     gXMautowrap, gStarMautowrap, gSnowMautowrap, fUMautowrap, gUMautowrap] = autowrapList

    # Retrieve x0 from pplusdose
    # x0_e = np.exp(x0)

    # Assign parameter values in x0_e to the parameter names in stateOList
    # First make a string called sText that can be executed
    # sText = assignParametersText(stateOList, "x0")
    # exec(sText)
    # print("Executed: %s" % sText)

    # Initial conditions wx, r0, s0 as specified by x0_ode_sensitivities
    wx = x[0:Ns] # first Ns are the initial state
    r = x[Ns:Ns * Nse + Ns].reshape(Ns, Nse) # r: derivative x wrt xo
    s = x[Ns * Nse + Ns:(Ns * Nsnow) + (Ns * Nse + Ns)].reshape(Ns, Nsnow)  # s: derivative x wrt thetaSnow
    u = x[(Ns * Nsnow) + (Ns * Nse + Ns):].reshape(Ns, (Nd+Npk))

    # Doses
    d = doseParms

    # PK parameters
    pks = pkParms

    # Stack star, concentrations, pk and snow parameters
    p = np.hstack((rParam, d, pks, thetaSnow))

    # Make f
    f = rhs(t, wx, p, d, doseParmsList, pks, PKparms, dose, fautowrap, Ns)

    # Make fT
    sText = makeTextRunAutowrap(autowrapName="fUMautowrap", statelistName="wx", parlistName="p", statelist=wx, parlist=p)
    utUN = eval(sText)
    # print("Executed: utUN = %s" % sText)
    utUN = np.array(utUN)

    # Make fX
    sText = makeTextRunAutowrap(autowrapName="fXMautowrap", statelistName="wx", parlistName="p", statelist=wx, parlist=p)
    fXN = eval(sText)
    # print("Executed: fXN = %s" % sText)

    # Make fStar
    sText = makeTextRunAutowrap(autowrapName="fStarMautowrap", statelistName="wx", parlistName="p", statelist=wx, parlist=p)
    fStarN = eval(sText)
    # print("Executed: fStarN = %s" % sText)
    fStarN = np.array(fStarN)

    # Make fSnow
    sText = makeTextRunAutowrap(autowrapName="fSnowMautowrap", statelistName="wx", parlistName="p", statelist=wx, parlist=p)
    fSnowN = eval(sText)
    # print("Executed: fSnowN = %s" % sText)
    fSnowN = np.array(fSnowN)

    # Make fRXO
    sText = makeTextRunAutowrap(autowrapName="fRXOMautowrap", statelistName="x0", parlistName="p", statelist=x0, parlist=p)
    fRXON = eval(sText)
    # print("Executed: fRXON = %s" % sText)
    fRXON = np.array(fRXON)

    # Make fRSnow
    sText = makeTextRunAutowrap(autowrapName="fRSnowMautowrap", statelistName="x0", parlistName="p", statelist=x0, parlist=p)
    fRSnowN = eval(sText)
    # print("Executed: fRSnowN = %s" % sText)
    fRSnowN = np.array(fRSnowN)

    fXR = np.dot(np.array(fXN.tolist(), dtype=float), np.array(r))
    # Check whether there are free initial/steady states. If not, fRXON should be a empty array of shape (Nstar,0)
    if Nse > 0:
        fStarfRXO = np.dot(np.array(fStarN.tolist(), dtype=float), np.array(fRXON))
    else:
        fRXON = np.empty((Nstar,0))
        fStarfRXO = np.dot(np.array(fStarN.tolist(), dtype=float), np.array(fRXON))
    fXS = np.dot(np.array(fXN.tolist(), dtype=float), np.array(s))
    fStarfRSnow = np.dot(np.array(fStarN.tolist(), dtype=float), np.array(fRSnowN))
    fXU = np.dot(np.array(fXN.tolist(), dtype=float), np.array(u))

    # Check whether the star matrices are not empty
    if Nstar == 0:
        dotR = fXR
        dotS = fSnowN + fXS
    else:
        dotR = fXR + fStarfRXO
        dotS = fSnowN + fXS + fStarfRSnow
    # Check whether fXU is empty
    if fXU.size:
        dotU = fXU + utUN
    else:
        dotU = fXU

    fall = np.hstack((np.hstack((np.hstack((f, (dotR.flatten()))), (dotS.flatten()))), dotU.flatten()))  # the extended term with both state and sensitivities
    # print("Dimensions f: %s" % (f.shape,))
    # print("Dimensions dotR: %s" % (dotR.shape, ))
    # print("Dimensions dotS: %s"% (dotS.shape, ))
    return fall

def map2Jy(stateList, allIniStatesList, stateOKnownList, stateOKnownNameList, paraStarList, paraSnowList, stateOList, obsList, paraList, doseParmsList, PKparms, pLog, tspan, Nm, Nso, Nd, Ns, Nse, Np, Nstar, Nsnow, Npk, dose, x0_ode_sensitivities, rlen, slen, autowrapList, matrixList):
    #pLog is: [log(initial states, x0), concentration(s), log(PKparameters), log(snowParameters), log(scalingParameters), offsetParameters]

    # Unravel the matrixList
    [f, fR, g] = matrixList

    # Unravel the autowrapList
    [fautowrap, gautowrap, fXMautowrap, \
     fStarMautowrap, fSnowMautowrap, fRXOMautowrap, fRSnowMautowrap, \
     gXMautowrap, gStarMautowrap, gSnowMautowrap, fUMautowrap, gUMautowrap] = autowrapList

    # Get the initial states from x0_ode_sensitivities and the snow parameters
    x0Log = x0_ode_sensitivities[0:Nse]
    thetaSnowLog = pLog[Nse+Nd+Npk:]

    # Calculate the exponent of all log-transformed initial values and snow parameters
    x0 = np.exp(x0Log)
    doseParms = pLog[Nse:Nse+Nd]
    pkParms = np.exp(pLog[Nse+Nd:Nse+Nd+Npk])
    thetaSnow = np.exp(thetaSnowLog)
    if Nso > 0:
        thetaSnow[-int(Nso/2):] = thetaSnowLog[-int(Nso/2):]
    dpksnowParms = np.append(doseParms,np.append(pkParms,thetaSnow))
    # print('x0:')
    # print(x0)
    # print('thetaSnow:')
    # print(thetaSnow)

    # Extract the known initial states and the x0 parameters from the dictionary. They should now be in the order of allIniStatesList
    # Make a dictionary of all the initial states and their current values
    d_tmp = {}
    for i, state in enumerate(stateOList):
        d_tmp[state] = x0[i]
    for i, state in enumerate(stateOKnownNameList):
        d_tmp[state] = stateOKnownList[i]
    x0_withKnowns = [d_tmp.get(key) for key in allIniStatesList]

    # Calculate the star parameters
    rParam = xOParam2Rparam(x0_withKnowns, dpksnowParms, Nso, paraStarList, paraSnowList, doseParmsList, PKparms, fR, allIniStatesList)
    # print("Starparameters:")
    # print(rParam)

    # Bind x0Log (initial states), thetaSnowLog (snow parameters) and the current dose
    # pplusdose = [x0Log, thetaSnowLog, dose]
    pplusdose = [x0, doseParms, pkParms, thetaSnow, dose]

    # Make a arguments list
    argList = [pplusdose, stateOList, doseParmsList, PKparms, Ns, Nse, Np, Nso, Nd, Npk, Nsnow, Nstar, rParam, autowrapList]

    # Run odeint from scipy.integrate package. Given a set of parameters passed by argList and initial states, the
    # set of ODE equations can be solved. This is done for every time point t in tspan. The ODEs of the initial system
    # of equations is solved together with the additional ODEs dotR and dotS.
    # y0 should be the parameters on linear scale.
    sol = odeint(func=fRightHandSide, y0=np.hstack([x0_withKnowns, x0_ode_sensitivities[Nse:]]), t=tspan, args=(argList,))
    # print("Dimensions sol: %s" % (sol.shape, ))
    # print(sol)

    # Extract the dynamics of f_i, i.e. the solution of the original ODEs, from the matrix
    stateSimu = sol[:, 0:Ns]

    # Split sol into the r and s and u matrices
    dotRmatrix = sol[:, Ns:Ns + (Ns * Nse)]
    dotSmatrix = sol[:, Ns + (Ns * Nse):Ns + (Ns * Nse) + (Ns * Nsnow)] # sol[:, Ns + Ns * Ns:]
    dotUmatrix = sol[:, Ns + (Ns * Nse) + (Ns * Nsnow):]

    # Initialise a matrix obsSimu. We calculated the state simulation from the ODE system in f and with the
    # parameter values and initial states, we can calculate g0 and g1. We save the output in a obsSimu matrix.
    obsSimu = np.array([]).reshape(0, Nm)

    # Make parameter set p
    p = np.hstack((rParam, doseParms, pkParms, thetaSnow))

    # gautowrap needs time <t> and <dose>, even though these are not used in gautowrap, so it doesn't matter what the value is.
    t = 0

    # Loop over the index of the time points
    for tindex in range(0,len(tspan)):
        # Pass the parameter values as arguments to gautowrap to calculate g0 and g1 for every time point
        sText = makeTextRunAutowrap(autowrapName="gautowrap",
                                    statelistName="stateSimu[tindex,:]", parlistName="p",
                                    statelist=stateSimu[tindex,:], parlist=p, withTimeDose = True)
        #print(sText)
        gautowrapTemp = eval(sText)
        # print("Executed: gautowrapTemp = %s" % sText)

        # Save the output of gautowrap saved in gautowrapTemp in an array called obsSimu
        gtn = np.array(gautowrapTemp.tolist(), dtype=float).reshape(Nm, )
        obsSimu = np.vstack((obsSimu,gtn)) # dim: Nt * Nm

    # Initialise a gix numpy array
    d1 = {}
    for i in range(0, Nm):
        varName = 'g' + str(i + 1) + 'x'
        d1[varName] = np.array([]).reshape(0, Nse+Nd+Npk+Nsnow) # np.array([]).reshape(0, Ns+Nsnow)

    # Now we have g0 and g1 saved in obsSimu, but we also need the solutions for the other partial differential equations
    # gSnow, gStar, gSnowStar etc.
    # Loop over the rows of sol, i.e. over the time points (matrix.shape[0] = nrows and matrix.shape[1] = ncols)
    for tindex in range(0, sol.shape[0]):
        # Get the dotRmatrix at time point tindex
        tempR = dotRmatrix[tindex, ].reshape(Ns, Nse)
        tempS = dotSmatrix[tindex, ].reshape(Ns, Nsnow)
        tempU = dotUmatrix[tindex, ].reshape(Ns, Nd+Npk)
        wx = stateSimu[tindex,]

        # Pass the parameter values as arguments to gXMautowrap
        sText = makeTextRunAutowrap(autowrapName="gXMautowrap",
                                    statelistName="wx", \
                                    parlistName="p", statelist=wx, parlist=p, withTimeDose=True)
        gXN = eval(sText)
        #print("Executed: gXN = %s" % sText)

        # Pass the parameter values as arguments to gSnowMautowrap
        sText = makeTextRunAutowrap(autowrapName="gSnowMautowrap",
                                    statelistName="wx", \
                                    parlistName="p", statelist=wx, parlist=p, withTimeDose=False)
        gSnowN = eval(sText)
        #print("Executed: gSnowN = %s" % sText)

        # Pass the parameter values as arguments to gUMautowrap
        sText = makeTextRunAutowrap(autowrapName="gUMautowrap",
                                    statelistName="wx", \
                                    parlistName="p", statelist=wx, parlist=p, withTimeDose=True)
        gUN = eval(sText)
        #print("Executed: gSnowN = %s" % sText)

        # Pass the parameter values as arguments to gStarMautowrap
        sText = makeTextRunAutowrap(autowrapName="gStarMautowrap",
                                    statelistName="wx", \
                                    parlistName="p", statelist=wx, parlist=p, withTimeDose=False)
        gStarN = eval(sText)
        #print("Executed: gStarN = %s" % sText)

        # Pass the parameter values as arguments to fRXOMautowrap
        t = -100
        sText = makeTextRunAutowrap(autowrapName="fRXOMautowrap",
                                    statelistName="x0", \
                                    parlistName="p", statelist=x0, parlist=p, withTimeDose=True)
        fRXON = eval(sText)
        #print("Executed: fRXON = %s" % sText)

        # Pass the parameter values as arguments to fRSnowMautowrap
        t = -100
        sText = makeTextRunAutowrap(autowrapName="fRSnowMautowrap",
                                    statelistName="x0", \
                                    parlistName="p", statelist=x0, parlist=p, withTimeDose=True)
        fRSnowN = eval(sText)
        #print("Executed: fRSnowN = %s" % sText)

        # Check whether there are free initial/steady states. If not, fRXON should be a empty array of shape (Nstar,0)
        if Nse > 0:
            fRXON = fRXON
        else:
            fRXON = np.empty((Nstar, 0))

        # Check whether the star matrices are not empty
        if len(rParam) == 0:
            # Calculate g0 and gSnow and gU
            gON = np.dot(gXN, tempR)
            gSnowN = gSnowN + np.dot(gXN, tempS)
        else:
            # Calculate g0 and gSnow and gU
            gON = np.dot(gXN, tempR) + np.dot(gStarN, fRXON)
            gSnowN = gSnowN + np.dot(gXN, tempS) + np.dot(gStarN, fRSnowN)

        # Check whether fXU is empty
        if gUN.size:
            gUN = np.dot(gXN, tempU) + gUN
        else:
            gUN = np.dot(gXN, tempU)

        # Here the transformation plays a role.
        # Transformation of dg/log(dX0)
        TO = np.diag(x0)
        # Transformation of dg/log(dSnow)
        TSnowList = thetaSnow
        if Nso > 0:
            TSnowList[-int(Nso/2):] = 1.0
        TSnow = np.diag(TSnowList)
        # Transformation of dg/log(dU)
        TUList = np.append(doseParms,pkParms)
        if Nd > 0:
            TUList[0:Nd] = 1.0
        TU = np.diag(TUList)

        gON = np.dot(gON, TO)
        gSnowN = np.dot(gSnowN, TSnow)
        gUN = np.dot(gUN, TU)

        # Output yg are the solutions to the g ODEs.
        yg = np.hstack((gON, gSnowN, gUN))
        #print(yg)

        for i in range(0,Nm):
            varName = 'g' + str(i+1) + 'x'
            d1[varName] = eval('np.vstack(( d1["g' + str(i+1) + 'x"], yg[' + str(i) + ',:]))')

    endN = Nm * dotRmatrix.shape[0]
    gall = np.zeros([endN, Nse+Nd+Npk+Nsnow])
    for i in range(0, Nm):
        gall[range(i, endN, Nm),:] = eval('d1["g' + str(i + 1) + 'x"]')

    return stateSimu, obsSimu, gall

def t0Sensitivities(Ns, Nse, Nsnow, Nd, Npk, stateOList, stateList):
    '''
    This function initializes the sensitivity matrices for both the states (x) and parameters (p)
    :param Ns: number of states
    :param Nsnow: number of purely free parameters
    :return:
    Two matrices: rmatrix and smatrix. rmatrix is the matrix dx/dxo and smatrix is the matrix dx/dthetaSnow
    '''
    # Make the rmatrix
    # Initial matrix dimensions
    rmatrix = np.zeros([Ns, Nse])
    # Get the state names for which the initial values are not known
    unknownStates = [str.split(stateO, '_')[0] for stateO in stateOList]
    # Check which of the states is unknown (True if known, False if not known)
    bool_list = [state in unknownStates for state in stateList]
    # Put 1 at correct positions in the rmatrix
    row = 0
    i = 0
    for b in bool_list:
        if b:
            rmatrix[row, i] = 1
            i += 1
        row += 1

    # Make the s and u matrices
    smatrix = np.zeros([Ns, Nsnow])
    umatrix = np.zeros([Ns, (Nd+Npk)])

    return(np.array(rmatrix, dtype = float), np.array(smatrix, dtype = float),np.array(umatrix, dtype = float))

def p2JacobianPrediction(pLog, Nso, Nd, Nm, Ns, Nse, Np, Nstar, Nsnow, Npk, tspan, doseList, paraNamesList, autowrapList, matrixList):
    # pLog is: [log(initial states, x0), log(snowParameters), log(scalingParameters), offsetParameters]

    # Unravel the list of lists of parameter names
    allIniStatesList, stateOKnownList, stateOKnownNameList, stateList, stateOList, obsList, paraStarList, paraSO, paraSnowList, paraList, paraFreeList, doseParmsList, PKparms = paraNamesList

    # Check whether the paraFreeList
    if len(paraFreeList) != len(pLog):
        print("The number of parameters in paraFreeList are not same as specified in pLog. Please check again.")
        sys.exit()

    # Make x0_ode_sensitivities, i.e. an array of the steady state values for all partial differential equations:
    # x0, r0, s0; or x0, Identity matrix, zero-matrix; or x0, I, O
    # Get the initial values for the state variables
    x0Log = pLog[range(0, Nse)]
    # Make the Identity and Zero matrix
    rMatrix, sMatrix, uMatrix = t0Sensitivities(Ns, Nse, Nsnow, Nd, Npk, stateOList, stateList)
    # Reshape the matrices a array
    rReshape = rMatrix.reshape(Ns * Nse)
    rlen = len(rReshape)
    sReshape = sMatrix.reshape(Ns * Nsnow)
    slen = len(sReshape)
    uReshape = uMatrix.reshape(Ns * (Nd + Npk))
    ulen = len(uReshape)
    x0_ode_sensitivities = np.hstack((x0Log, rReshape, sReshape, uReshape))

    # Initialise the final arrays that are filled and serve as the output of the function
    gall_allDoses = np.array([]).reshape(0, Nse + Nd + Npk + Nsnow)
    stateSimu_allDoses = np.array([]).reshape(0, Ns)
    obsSimu_allDoses = np.array([]).reshape(0, Nm)

    # Run map2Jy for every dose separately. The output of map2Jy per dose will be combined in a matrix, where
    # every row is the output of one dose
    if len(doseList) > 0:
        for dosei in doseList:
            # gall is the jacobian matrix
            stateSimu, obsSimu, gall = map2Jy(stateList, allIniStatesList, stateOKnownList, stateOKnownNameList, paraStarList, paraSnowList, stateOList, obsList, paraList, doseParmsList, PKparms, pLog, \
                   tspan, Nm, Nso, Nd, Ns, Nse, Np, Nstar, Nsnow, Npk, dosei, x0_ode_sensitivities, rlen, slen, autowrapList, matrixList)

            # Add the partial differential equations excluded from the Sensitivity equation, i.e. the concentration parameters and pharmocokinetic parameters
            gall_states = gall[:,0:Nse]
            gall_snow = gall[:,Nse:Nse+Nsnow] # gall[:,Nse:]
            gall_doses = gall[:,Nse+Nsnow:Nse+Nsnow+Nd] # np.zeros((gall.shape[0],Nd))
            gall_pk = gall[:,Nse+Nsnow+Nd:Nse+Nsnow+Nd+Npk] # np.zeros((gall.shape[0],Npk))

            #print(gall_doses)

            # Bind all the gall parts
            gallNew = np.hstack((np.hstack((np.hstack((gall_states,gall_doses)),gall_pk)),gall_snow))

            # Combine the output per dose in a matrix, where every row is the output of one dose
            gall_allDoses = np.vstack((gall_allDoses, gallNew))
            stateSimu_allDoses = np.vstack((stateSimu_allDoses, stateSimu))
            obsSimu_allDoses = np.vstack((obsSimu_allDoses, obsSimu))
    else:
        # gall is the jacobian matrix
        stateSimu, obsSimu, gall = map2Jy(stateList, allIniStatesList, stateOKnownList, stateOKnownNameList,
                                          paraStarList, paraSnowList, stateOList, obsList, paraList, doseParmsList,
                                          PKparms, pLog, tspan, Nm, Nso, Nd, Ns, Nse, Np, Nstar, Nsnow, Npk, 0,
                                          x0_ode_sensitivities, rlen, slen, autowrapList, matrixList)

        # Add the partial differential equations excluded from the Sensitivity equation, i.e. the concentration parameters and pharmocokinetic parameters
        gall_states = gall[:, 0:Nse]
        gall_snow = gall[:, Nse:Nse + Nsnow]  # gall[:,Nse:]
        gall_doses = gall[:, Nse + Nsnow:Nse + Nsnow + Nd]  # np.zeros((gall.shape[0],Nd))
        gall_pk = gall[:, Nse + Nsnow + Nd:Nse + Nsnow + Nd + Npk]  # np.zeros((gall.shape[0],Npk))

        # print(gall_doses)

        # Bind all the gall parts
        gallNew = np.hstack((np.hstack((np.hstack((gall_states, gall_doses)), gall_pk)), gall_snow))

        # Combine the output per dose in a matrix, where every row is the output of one dose
        gall_allDoses = np.vstack((gall_allDoses, gallNew))
        stateSimu_allDoses = np.vstack((stateSimu_allDoses, stateSimu))
        obsSimu_allDoses = np.vstack((obsSimu_allDoses, obsSimu))

    # Create a symmetric matrix J
    J = np.dot(np.transpose(gall_allDoses), gall_allDoses)

    return stateSimu_allDoses, obsSimu_allDoses, gall_allDoses, J

def differencePreData(pLog, argP):
    # pLog is: [log(initial states, x0), log(snowParameters), log(scalingParameters), offsetParameters]
    # Make gall_allDoses and obsSimu_allDoses global variables, to use it in p2JacobianLsq
    global gall_allDoses, obsSimu_allDoses

    # Unravel the argP components
    [Nso, Nd, Nm, Ns, Nse, Np, Nstar, Nsnow, Npk, Nt, tspan, doseList, plateID_list, dataM, dataS, dataR, optim, paraNamesList, autowrapList, matrixList] = argP

    # Run p2JacobianPrediction
    stateSimu_allDoses, obsSimu_allDoses, gall_allDoses, J = p2JacobianPrediction(pLog, Nso, Nd, Nm, Ns, Nse, Np, Nstar, Nsnow, Npk, tspan, doseList, paraNamesList, autowrapList, matrixList)

    if len(doseList) == 0 and len(plateID_list) == 0:
        difference = ((obsSimu_allDoses.reshape(Nm * Nt, ) - dataM))
    elif len(doseList) == 0 and len(plateID_list) > 0:
        difference = ((obsSimu_allDoses.reshape(Nm * Nt, ) - dataM) / dataS)
    elif len(doseList) > 0 and len(plateID_list) == 0:
        difference = ((obsSimu_allDoses.reshape(Nm * Nt * len(doseList), ) - dataM))
    elif len(doseList) > 0 and len(plateID_list) > 0:
        difference = ((obsSimu_allDoses.reshape(Nm * Nt * len(doseList), ) - dataM) / dataS)

    # Get squared residuals
    sqdist = []
    for data_repl in dataR:
        if len(doseList) == 0:
            # Calculate squared difference
            diff = (obsSimu_allDoses.reshape(Nm * Nt, ) - data_repl) ** 2
            # Replace nans with zeros
            diff[np.isnan(diff)] = 0
            sqdist.append(diff)
        elif len(doseList) > 0:
            # Calculate squared difference
            diff = (obsSimu_allDoses.reshape(Nm * Nt * len(doseList), ) - data_repl) ** 2
            # Replace nans with zeros
            diff[np.isnan(diff)] = 0
            sqdist.append(diff)

    # Get sum
    sumsqdist = sum(sqdist)
    # Get MSE
    MSE = sumsqdist / len(sqdist)

    if optim == "MSE":
        out = MSE
    elif optim == "difference":
        out = difference

    return out # MSE # difference

def p2JacobianLsq(pLog, argP):
    # Unravel the argP components
    [Nso, Nd, Nm, Ns, Nse, Np, Nstar, Nsnow, Npk, Nt, tspan, doseList, plateID_list, dataM, dataS, dataR, optim, paraNamesList, autowrapList, matrixList] = argP

    if len(plateID_list) > 0:
        jaclsq = gall_allDoses * np.tile(((1.0 / (dataS)) ).reshape(len(dataS), 1), [1, len(pLog)])
    elif len(plateID_list) == 0:
        jaclsq = gall_allDoses

    return jaclsq

def runMultiesti(zippedArgs):
    # Run multiesti with unzipped arguments
    print("Running estimation")
    return multiesti(*zippedArgs)

def multiesti(seed, argPplus):
    [pInit, argP, testRun, startTime, timeMax, path, DATE, MODEL, bounds] = argPplus
    np.random.seed(seed)
    # The bounds for all parameters are unlimited: from -Inf to Inf
    #boundub = ([-np.inf] * len(pInit[seed]), [np.inf] * len(pInit[seed]))
    boundub = bounds

    [Nso, Nd, Nm, Ns, Nse, Np, Nstar, Nsnow, Npk, Nt, tspan, doseList, plateID_list, dataM, dataS, dataR, optim, paraNamesList, autowrapList, matrixList] = argP

    [f, fR, g] = matrixList

    allIniStatesList, stateOKnownList, stateOKnownNameList, stateList, stateOList, obsList, paraStarList, paraSO, paraSnowList, paraList, paraFreeList, doseParmsList, PKparms = paraNamesList

    p = pInit[seed]
    #print(p)

    # Log-transform all the parameter values, except the offset parameters and the doses
    # (where int(Nso/2) = the number of offset parameters)
    # pLog is: [log(initial states, x0), log(concentration(s)),log(snowParameters), log(scalingParameters), offsetParameters]
    pLog = p.copy()
    # Log-transform the initial states
    pLog[0:Nse] = np.log(p[0:Nse])
    # Log-transform the pk, snow and scaling parameters
    pLog[Nse+Nd:len(p)-int(Nso/2)] = np.log(p[Nse+Nd:len(p)-int(Nso/2)])
    print(pLog)

    # Parameter values at the start and end are the same before estimation
    pStart = pLog
    pEnd = pLog

    cost = 1.0e36

    if optim == "MSE":
        jacobian = '2-point'
    elif optim == "difference":
        jacobian = p2JacobianLsq

    #try:
    if testRun:
        print("Only running a testrun to check whether least_squares doesn't give an error... ")
        least_squares(differencePreData, pLog, jac=p2JacobianLsq, bounds=boundub, args=(argP,), xtol=1.0e-10,
                      ftol=1.0e-10, verbose=2)
    else:
        # Do parameter estimation the first time
        est = least_squares(differencePreData, pLog, jac=jacobian, bounds=boundub, args=(argP,), xtol=1.0e-5,
                        ftol=1.0e-5, verbose=2)

        # Update parameters
        # print("first Log update:")
        # print(est.x)
        pLog = est.x

        # Do parameter estimation the second time
        est = least_squares(differencePreData, pLog, jac=jacobian, bounds=boundub, args=(argP,), xtol=1.0e-5,
                            ftol=1.0e-5, verbose=2)

        # Retrieve the cost of cost function at estimation round 2
        cost = est.cost

        # Update parameters
        pLog = est.x

        # Update final parameter values
        pEnd = pLog

        # Do parameter estimation until convergence to minimum
        flagContinue = 1
        count = 0
        cycle = 0

        while flagContinue and count < 10:
            print("Cycle = %i" % count)
            # Do parameter estimation the count-th time
            est = least_squares(differencePreData, pLog, jac=jacobian, bounds=boundub, args=(argP,), xtol=1.0e-5,
                                ftol=1.0e-5, verbose=2)

            # Check whether the cost still decreases
            count += 1
            # If two subsequent costs are very similar, estimation stops
            if cost - est.cost < 1.0e-5 or (time.time() - startTime > timeMax):
                flagContinue = 0
            # If cost is still decreasing a lot, the count is reset to 0
            elif cost - est.cost > 1.0e-4:
                count = 0
            
            if est.cost > 1.0e3 and cycle > 3:
                flagContinue = 0

            # Update parameters
            pLog = est.x
            # Update cost
            cost = est.cost
            # Update final parameter values
            pEnd = pLog

            cycle += 1

    #except:
    #    sys.exit("least_squares gives a fatal error. Please check you input functions and debug. ")

    d = {'pStart': pStart, 'pEnd': pEnd, 'cost': cost}

    pathToCSVfile = path + DATE + "_MH_" + MODEL + "Model_FinishedRuns_costs.csv"
    with open(pathToCSVfile, "a") as g:
        #fcntl.flock(g, fcntl.LOCK_EX)
        g.write("parset_" + str(seed + 1) + "," + str(seed + 1) + "," + str(cost) + "\n")
        g.close()
        #fcntl.flock(g, fcntl.LOCK_UN)

    print("Estimation of parameter set %i done! Final cost was %f" % ((seed+1), cost))

    # Make a df of the initial parameters
    parmsDf = pd.DataFrame({'init_value': p}, index=paraFreeList)
    # Calculate the star parameters from the (initial or) steady states and other parameters
    pInitAllDf = transformParms(transformParms=paraStarList, useParms=stateOList + paraSnowList, df=parmsDf, matrix=fR)

    # Transform the output parameters back to linear scale
    outputParmsLog = pEnd
    outputParms = np.exp(outputParmsLog)
    outputParms[Nse:Nse+Nd] = outputParmsLog[Nse:Nse+Nd]
    if Nso > 0:
        outputParms[-int(Nso / 2)::] = outputParmsLog[-int(Nso / 2)::]

    # Make a data frame of the new parameters
    newParmsDf = pd.DataFrame({'est_value':outputParms}, index = paraFreeList)

    # Calculate the star parameters from the (initial or) steady states and other parameters
    outputParmsDf = transformParms(transformParms=paraStarList, useParms=stateOList + paraSnowList, df=newParmsDf, matrix=fR)
    outputDf = pd.concat([pInitAllDf,outputParmsDf], axis = 1)
    pAllNew = outputParmsDf["est_value"].values

    # Save the new data frame
    print("outputParmsDf:")
    print(outputParmsDf)

    try:
        outputDf.to_csv(path + DATE + "_MH_" + MODEL +"Model_parameterEstimates_parmset" + str(seed+1) + "_cost_" + "{0:.2f}".format(cost) + ".csv")
        print("Saved %s" % (path + DATE + "_MH_" + MODEL +"Model_parameterEstimates_parmset" + str(seed+1) + "_cost_" + "{0:.2f}".format(cost) + ".csv"))
    except:
        outputDf.to_csv(path + DATE + "_MH_" + MODEL +"Model_parameterEstimates_parmset" + str(seed+1) + "_cost_" + str(cost) + ".csv")
        print("Saved %s" % (path + DATE + "_MH_" + MODEL +"Model_parameterEstimates_parmset" + str(seed+1) + "_cost_" + str(cost) + ".csv"))
    return pStart, pEnd, cost, d
 
