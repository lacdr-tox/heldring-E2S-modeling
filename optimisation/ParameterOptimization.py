###################################################
# Part 0 - Import modules
###################################################

# Import all packages
import copy
import os
os.environ["MKL_CBWR"] = "SSE2"
os.environ["OMP_NUM_THREADS"] = '1'
import matplotlib
matplotlib.use('Agg')
from matplotlib import colors
from pyDOE import *
from odeFunctions import *
from multiprocessing import Pool, freeze_support
import itertools
import re
import sys
import pandas as pd
import platform
import warnings
from datetime import datetime

# Use sys arguments or arguments defined here
useSysArgv = True
# Choose from "MSE" or "difference"
optim = "difference"

if __name__ == '__main__':
    from sympy.utilities.autowrap import autowrap

    freeze_support()

    print("Script started...")

    print("All packages imported")
    print("Python version %s.%s.%s" %(sys.version_info[0], sys.version_info[1], sys.version_info[2]))
    computerName = platform.node()

    print("This is the name of the script: %s" % sys.argv[0])

    ###################################################
    # Part 0 - Set user variables
    ###################################################

    # Set the model name, the date and the run index, to prevent overwriting output from previous runs
    if not useSysArgv:
        # Model name
        MODEL = "M006"

        # Define the path to the input folder
        pathToInputFolder = "/Users/muriel/Documents/LACDR/Projects/PHH/Models/v2_firstRevision/"# "/Users/muriel/Documents/LACDR/Methods/00_OptimPackage/TestModel/"

        # Names of the input files
        # The data file MUST contain a name of 1-5 characters/digits containing info on the data origin
        # (e.g. the treatment) concatenated to 'data'. For example "DMSOdata" or "SIMUdata"
        dataFile = "20201006_MH_CDDPdata_wFUCCI_DDP.csv"
        odeFile = "M006_20220308_ODE_p53model.txt"

        # Paths to output folder
        pathToOutput = "/Users/muriel/Documents/LACDR/Methods/00_OptimPackage/TestModel/TestOutput/" + MODEL + "/"

    else:
        # Get dataDile and odeFile
        try:
            dataPathAndFile = sys.argv[1]
            odePathAndFile = sys.argv[2]

            # Create the path and the file
            m1 = re.search("/", dataPathAndFile)
            if m1:
                pathToInputFolder = "/".join(dataPathAndFile.split("/")[0:-1]) + "/"
                dataFile = dataPathAndFile.split("/")[-1]
            else:
                pathToInputFolder = ""
                dataFile = dataPathAndFile

            m2 = re.search("/", dataPathAndFile)
            if m2:
                pathToInputFolder = "/".join(odePathAndFile.split("/")[0:-1]) + "/"
                odeFile = odePathAndFile.split("/")[-1]
            else:
                pathToInputFolder = ""
                odeFile = odePathAndFile

            pathToOutput = sys.argv[3]
        except IndexError:
            sys.exit("Please enter the (path to the) dataFile and/or odeFile as sys.argv.")

        # Get the model name
        try:
            MODEL = sys.argv[4]
        except IndexError:
            warnings.warn("No model name passed. Getting the model name from the odeFile. ")
            MODEL = odeFile.split(".txt")[0]

        pathToOutput = pathToOutput + "/" + MODEL + "/"

    DATETIME = datetime.today().strftime('%Y%m%d_%H%M%S')

    # Colors for the plot (should have same length as len(plateID_list)). Names should be part of the Python colors.
    palette = ["steelblue", "steelblue", "steelblue", "steelblue", 'r', 'g']  # ,'grey',"r",'y',"g"]

    # Create a list of markers to use in the plots
    markers = ["o", "v", "s", "*", "d", "X", "^", "<", ">", "o", "v", "s", "*", "d", "X", "^", "<", ">"]

    # Do a testRun?
    testRun = False

    ###################################################
    # Part I - Create directories and save information
    # on input and output paths and files.
    ###################################################

    # Print some of the settings
    print("This is the current date and time: %s" % DATETIME)

    # Make paths to the files
    pathToData = pathToInputFolder + dataFile
    pathToODEs = pathToInputFolder + odeFile

    # Create directories if they do not exist
    if os.path.isdir(pathToOutput):
        print(pathToOutput + " does already exist. Overwriting previous output...")
    else:
        print(pathToOutput + " does not exist. Creating a new directory and writing the output to " + pathToOutput)
        os.mkdir(pathToOutput)

    path0 = pathToOutput + re.findall("[a-zA-Z0-9]{1,5}data", dataFile)[0].split("data")[0] + "/"
    if os.path.isdir(path0):
        print(path0 + " does already exist. Overwriting previous output...")
    else:
        print(path0 + " does not exist. Creating a new directory and writing the output to " + path0)
        os.mkdir(path0)

    path = path0 + DATETIME + "/"
    if os.path.isdir(path):
        print(path + " does already exist. Overwriting previous output...")
    else:
        print(path + " does not exist. Creating a new directory and writing the output to " + path)
        os.mkdir(path)

    ###################################################
    # Part II - Read in the parameters and user settings
    # check whether most sets are defined correctly.
    ###################################################

    # Read in the user settings and parameters
    try:
        fh = open(pathToODEs, 'r')
        print("Path to text file with ODEs and parms exists")
        with open(pathToODEs) as f:
            lines = f.readlines()
        startReading = False
        stopReading = False

        print("The user settings:")
        for i, line in enumerate(lines):

            match1 = re.search("####START_USER_SETTINGS####", line)
            match2 = re.search("####END_PARAMETERS####\n", line)
            if match1:
                startReading = True
            elif match2:
                stopReading = True

            if startReading and not stopReading:
                print(line)
                exec(line)

    except FileNotFoundError:
        sys.exit("Path to text file with user settings, ODEs and parameters does not exist.")

    print("Number of estimation runs: %i" % NiniP)
    print("Number of CPUs: %i" % Ncpu)
    print("Maximum estimation time: %f hours" % (timeMax/60/60))

    # Save the metadata
    Metadata = ["MODEL = \"" + MODEL + "\"",
                "DATETIME = \"" + DATETIME + "\"",
                "pathToInputFolder = \"" + pathToInputFolder + "\"",
                "dataFile = \"" + dataFile + "\"",
                "pathToData = \"" + pathToData + "\"",
                "odeFile = \"" + odeFile + "\"",
                "pathToODEs = \"" + pathToODEs + "\"",
                "pathToOutput = \"" + pathToOutput + "\"",
                "NiniP = " + str(NiniP),
                "Ncpu = " + str(Ncpu)]
    # Write all metadata to a file
    with open(path + "metadata.txt","w+") as h:
        for listitem in Metadata:
            h.write('%s\n' % listitem)
    h.close()
    print("Metadatafile saved in %smetadata.txt" % path)

    # Check whether obsList and obsStateVarList are of the same length:
    if not len(obsList) == len(obsStateVarList):
        sys.exit("Oops! obsList and obsStateVarList as defined in the %s file are not of the same length." % (odeFile))

    ###################################################
    # Part III - Check whether input in <user settings>
    # corresponds to input in <parameters> part.
    ###################################################

    # Check data input settings. StateVarName, doseColName and plateIDName must exist.
    # Check the state variable column name input
    error1 = "Error: your <%s> is not specified in file %s. Please indicate the name of the column in your data that contains the information on %s."
    error2 = "Error: <%s> is not defined in file %s. Define either the name, or, if there is no information on %s, define <%s> as None."

    if not 'stateVarName' in globals():
        sys.exit(error1 % ("stateVarName", odeFile,"state variables"))
    elif stateVarName is None:
        sys.exit(error1 % ("stateVarName", odeFile,"state variables"))

    if not 'doseColName' in globals():
        sys.exit(error2 % ("doseColName", odeFile,"concentrations", "doseColName"))

    if not 'plateIDName' in globals():
        sys.exit(error2 % ("plateIDName", odeFile,"plate/experiment/replicate IDs", "plateIDName"))

    # Check the dose input
    dosesFromData = False
    w = "Your <doseList> contains dose information, but the column name stored in <doseColName> is not specified. The <doseList> is set to an empty list. "
    if doseColName is None:
        warnings.warn("Your <doseColName> is set to None. The script continues with the assumption that you do not have multiple concentrations to fit to.")
        if not 'doseList' in globals():
            warnings.warn(w)
            doseList = []
        elif doseList is None:
            doseList = []
        elif len(doseList) > 0:
            warnings.warn(w)
            doseList = []
    elif not doseColName is None:
        if not 'doseList' in globals():
            warnings.warn("No <doseList> is specified. <doseList> will be retrieved from the data.")
            dosesFromData = True
        elif doseList is None or doseList == []:
            warnings.warn("No <doseList> is specified. <doseList> will be retrieved from the data.")
            dosesFromData = True

    # Check the plateID input
    plateIDsFromData = False
    w = "Your <plateID_list> contains plate information, but the column name stored in <plateIDName> is not specified. The <plateID_list> is set to an empty list."
    if plateIDName is None:
        warnings.warn("Your <plateIDName> is set to None. The script continues with the assumption that you do not have multiple experimental replicates to fit to.")
        if not 'plateID_list' in globals():
            plateID_list = []
        elif plateID_list is None:
            plateID_list = []
        elif len(plateID_list) > 0:
            warnings.warn(w)
            plateID_list = []
    elif not plateIDName is None:
        if not 'plateID_list' in globals():
            warnings.warn("No <plateID_list> is specified. <plateID_list> will be retrieved from the data.")
            plateIDsFromData = True
        elif plateID_list is None or plateID_list == []:
            warnings.warn("No <plateID_list> is specified. <plateID_list> will be retrieved from the data.")
            plateIDsFromData = True

    ###################################################
    # Part IV - Read the data
    ###################################################

    # Get the file name
    filename = pathToData

    # Check whether file exists
    try:
        # Read in the data
        data = pd.read_csv(filename)
        print("Path to the data file exists. Reading data...")
    except FileNotFoundError:
        sys.exit("Error: Path to the data file does not exist.")

    # Check if data column names and input correspond.
    message = "Error: Requested column name \"%s\" does not exist in the data. Adjust the data or the <doseColName> in file %s."
    # Check if stateVarName exists in the data
    try:
        data[stateVarName]
    except KeyError:
        sys.exit(message % (stateVarName, odeFile))
    # Check if doseColName exists in the data
    try:
        data[doseColName]
    except KeyError:
        if doseColName is None:
            pass
        else:
            sys.exit(message % (doseColName, odeFile))
   # Check if plateIDName exists in the data
    try:
        data[plateIDName]
    except KeyError:
        if plateIDName is None:
            pass
        else:
            sys.exit(message % (plateIDName, odeFile))
    # Check if plateIDName exists in the data
    try:
        data[realDataColName]
    except KeyError:
        sys.exit(message % (realDataColName, odeFile))
    # Check if plateIDName exists in the data
    try:
        data[realTimeColName]
    except KeyError:
        sys.exit(message % (realTimeColName, odeFile))
    # Check if plateIDName exists in the data
    try:
        data[interpolDataColName]
    except KeyError:
        sys.exit(message % (interpolDataColName, odeFile))
    # Check if plateIDName exists in the data
    try:
        data[interpolTimeColName]
    except KeyError:
        sys.exit(message % (interpolTimeColName, odeFile))

    # Get the doses as defined in the data
    if dosesFromData:
        # Get the doses from the data
        print("Retrieving concentrations from the data...")
        doseList = list(data[doseColName].unique())
        print("Found the following concentrations: {}".format(doseList))

    # Get the plateIDs as defined in the data
    if plateIDsFromData:
        print("Retrieving plateIDs from the data...")
        plateID_list = []
        for sv in obsStateVarList:
            tmp = list(data[plateIDName].unique())
            plateID_list.append(tmp)
        print("Found the following plateIDs: {}".format(plateID_list))

    # Get the file name
    filename = pathToData

    try:
        # Read in the data
        data = pd.read_csv(filename)
        print("Path to the data file exists. Reading data...")

    except FileNotFoundError:
        sys.exit("Error: Path to the data file does not exist.")

    # Make a separate dataTemp data frame per cell line, because every cell line can have a different number
    # of replicates and save those data frames in a list:
    dfListInterpol = []
    dfListReal = []
    # Loop over the reporters
    for cell_linei, celllinet in enumerate(obsStateVarList):
        # Create NrOfDoses matrices with NrOfReplicates rows and NrOfTimepoints columns
        # For example: 22 time points, 3 replicates, 4 reporter cell lines, 5 doses

        # If there are multiple doses and replicates, do the following:
        # The output is a list of the state variables with, per state variable,
        # a list that contains state variable data per dose as lists,
        # with a list of replicates per dose
        if len(plateID_list) > 0 and len(doseList) > 0:
            dataTempInterpol = np.zeros((len(doseList), len(plateID_list[cell_linei]), nrOfTps))
            dataTempReal = np.zeros((len(doseList), len(plateID_list[cell_linei]), nrOfTps))
            # Loop over the doses
            for dosei, stresslevel in enumerate(doseList):
                # Loop over the replicates
                for replicatei, replicateID in enumerate(plateID_list[cell_linei]):

                    # Check if celllinet, stresslevel and replicateID exist in the data frame
                    if data[data[stateVarName] == celllinet].empty:
                        sys.exit("Oops!  Requested state variable %s does not exist in the data. Adjust the data or the obsStateVarList in file %s." % (celllinet, pathToODEs))
                    elif data[data[doseColName] == stresslevel].empty:
                        sys.exit("Oops!  Requested dose %s does not exist in the data. Adjust the data or the doseList in file %s." % (stresslevel, pathToODEs))
                    elif data[data[plateIDName] == replicateID].empty:
                        sys.exit("Oops!  Requested plateID %s does not exist in the data. Adjust the data or the obsStateVarList in file %s." % (replicateID, pathToODEs))

                    # Get the time points for this part of the data
                    a = (data[doseColName] == stresslevel).tolist()
                    b = (data[stateVarName] == celllinet).tolist()
                    c = (data[plateIDName] == replicateID).tolist()
                    selection = [ai and bi and ci for ai, bi, ci in zip(a, b, c)]
                    tps = list(data[selection][interpolTimeColName])[0:nrOfTps]

                    # Loop over the time points
                    for tpi,tp in enumerate(tps):
                        datasample = data[selection].iloc[[tpi]]

                        # Check if the data that is selected is unique
                        if not datasample.shape[0] == 1:
                            sys.exit("Oops!  Multiple data rows are selected with state variable %s, replicate ID %s, dose %s and time point %s." % (
                                    celllinet, replicateID, stresslevel, tp))

                        dataTempReal[dosei, replicatei, tpi] = datasample[realDataColName]
                        dataTempInterpol[dosei, replicatei, tpi] = datasample[interpolDataColName]

            # Get data timepoints
            dataTimePoints = data[selection][interpolTimeColName][0:nrOfTps]
            realTimePoints = data[selection][realTimeColName][0:nrOfTps]

        # If there are multiple doses and but no replicates, do the following.
        # The output is a list of the state variables with, per state variable,
        # a list that contains state variable data per dose as lists
        elif len(plateID_list) == 0 and len(doseList) > 0:
            dataTempInterpol = np.zeros((len(doseList), nrOfTps))
            dataTempReal = np.zeros((len(doseList), nrOfTps))

            # Loop over the doses
            for dosei, stresslevel in enumerate(doseList):
                # Check if celllinet, stresslevel and replicateID exist in the data frame
                if data[data[stateVarName] == celllinet].empty:
                    sys.exit(
                        "Error: Requested state variable %s does not exist in the data. Adjust the data or the obsStateVarList in file %s." % (
                        celllinet, pathToODEs))
                elif data[data[doseColName] == stresslevel].empty:
                    sys.exit(
                        "Error: Requested dose %s does not exist in the data. Adjust the data or the doseList in file %s." % (
                        stresslevel, pathToODEs))

                # Get the time points for this part of the data
                a = (data[doseColName] == stresslevel).tolist()
                b = (data[stateVarName] == celllinet).tolist()
                selection = [ai and bi for ai, bi in zip(a, b)]
                tps = list(data[selection][interpolTimeColName])[0:nrOfTps]
                # Loop over the time points
                for tpi, tp in enumerate(tps):
                    datasample = data[selection].iloc[[tpi]]

                    # Check if the data that is selected is unique
                    if not datasample.shape[0] == 1:
                        sys.exit(
                            "Error: Multiple data rows are selected with state variable %s, dose %s and time point %s." % (
                                celllinet, stresslevel, tp))

                    dataTempReal[dosei, tpi] = datasample[realDataColName]
                    dataTempInterpol[dosei, tpi] = datasample[interpolDataColName]

            # Get data timepoints
            dataTimePoints = data[selection][interpolTimeColName][0:nrOfTps]
            realTimePoints = data[selection][realTimeColName][0:nrOfTps]

        # If there are multiple replicates and but no doses, do the following:
        # The output is a list of the state variables with, per state variable,
        # a list that contains state variable data per replicate as lists
        elif len(plateID_list) > 0 and len(doseList) == 0:
            dataTempInterpol = np.zeros((len(plateID_list[cell_linei]), nrOfTps))
            dataTempReal = np.zeros((len(plateID_list[cell_linei]), nrOfTps))

            # Loop over the replicates
            for replicatei, replicateID in enumerate(plateID_list[cell_linei]):
                # Check if celllinet, stresslevel and replicateID exist in the data frame
                if data[data[stateVarName] == celllinet].empty:
                    sys.exit(
                        "Error: Requested state variable %s does not exist in the data. Adjust the data or the obsStateVarList in file %s." % (
                            celllinet, pathToODEs))
                elif data[data[plateIDName] == replicateID].empty:
                    sys.exit(
                        "Error: Requested plateID %s does not exist in the data. Adjust the data or the obsStateVarList in file %s." % (
                        replicateID, pathToODEs))

                # Get the time points for this part of the data
                a = (data[plateIDName] == replicateID).tolist()
                b = (data[stateVarName] == celllinet).tolist()
                selection = [ai and bi for ai, bi in zip(a, b)]
                tps = list(data[selection][interpolTimeColName])[0:nrOfTps]
                # Loop over the time points
                for tpi, tp in enumerate(tps):
                    datasample = data[selection].iloc[[tpi]]

                    # Check if the data that is selected is unique
                    if not datasample.shape[0] == 1:
                        sys.exit(
                            "Error: Multiple data rows are selected with state variable %s, replicate ID %s and time point %s." % (
                                celllinet, replicateID, tp))

                    dataTempReal[replicatei, tpi] = datasample[realDataColName]
                    dataTempInterpol[replicatei, tpi] = datasample[interpolDataColName]

            # Get data timepoints
            dataTimePoints = data[selection][interpolTimeColName][0:nrOfTps]
            realTimePoints = data[selection][realTimeColName][0:nrOfTps]

        elif len(plateID_list) == 0 and len(doseList) == 0:
            dataTempInterpol = np.zeros((nrOfTps))
            dataTempReal = np.zeros((nrOfTps))

            # Check if celllinet, stresslevel and replicateID exist in the data frame
            if data[data[stateVarName] == celllinet].empty:
                sys.exit(
                    "Error: Requested state variable %s does not exist in the data. Adjust the data or the obsStateVarList in file %s." % (
                        celllinet, pathToODEs))

            # Get the time points for this part of the data
            selection = (data[stateVarName] == celllinet).tolist()
            tps = list(data[selection][interpolTimeColName])[0:nrOfTps]
            # Loop over the time points
            for tpi, tp in enumerate(tps):
                datasample = data[selection].iloc[[tpi]]

                # Check if the data that is selected is unique
                if not datasample.shape[0] == 1:
                    sys.exit(
                        "Error: Multiple data rows are selected with state variable %s and time point %s." % (celllinet, tp))

                dataTempReal[tpi] = datasample[realDataColName]
                dataTempInterpol[tpi] = datasample[interpolDataColName]

            # Get data timepoints
            dataTimePoints = data[selection][interpolTimeColName][0:nrOfTps]
            realTimePoints = data[selection][realTimeColName][0:nrOfTps]
        else:
            sys.exit("Error: Something wrong with the plateID_list or doseList...")


        dfListInterpol.append(dataTempInterpol)
        dfListReal.append(dataTempReal)

    dfListMeanReal = []
    dfListStdReal = []
    dfListMean = []
    dfListStd = []

    dfReps = []
    dataR = []

    # Loop over the dataframes, i.e. the state variables
    for cell_linei, celllinet in enumerate(obsStateVarList):
        df = dfListInterpol[cell_linei]
        dfReal = dfListReal[cell_linei]

        if len(plateID_list) > 0 and len(doseList) > 0:
            dfTemp = np.array([]).reshape(len(plateID_list[cell_linei]), 0)
            dfTempReal = np.array([]).reshape(len(plateID_list[cell_linei]), 0)
            for dosei in np.arange(0, len(doseList), 1):
                # Bind the doses horizontally to obtain a NrOfReplicates x (NrOfDoses x NrOfTimepoints) matrix
                dfTemp = np.hstack((dfTemp, df[dosei]))
                dfTempReal = np.hstack((dfTempReal, dfReal[dosei]))
        elif len(plateID_list) == 0 and len(doseList) > 0:
            dfTemp = np.array([]).reshape(0)
            dfTempReal = np.array([]).reshape(0)
            for dosei in np.arange(0, len(doseList), 1):
                # Bind the doses horizontally to obtain a NrOfReplicates x (NrOfDoses x NrOfTimepoints) matrix
                dfTemp = np.hstack((dfTemp, df[dosei]))
                dfTempReal = np.hstack((dfTempReal, dfReal[dosei]))
            dfTemp = np.array([dfTemp])
            dfTempReal = np.array([dfTempReal])
        elif len(plateID_list) > 0 and len(doseList) == 0:
            dfTemp = df
            dfTempReal = dfReal
        elif len(plateID_list) == 0 and len(doseList) == 0:
            dfTemp = np.array([df])
            dfTempReal = np.array([dfReal])
        else:
            sys.exit("Error: Something wrong with the plateID_list or doseList")

        # Save the replicates
        dfReps.append(dfTemp)

        # Calculate the mean and std per column, i.e. the mean and std of the replicates
        dfTempM = np.mean(dfTemp, axis=0)
        dfTempS = np.std(dfTemp, axis=0)
        dfTempMReal = np.mean(dfTempReal, axis=0)
        dfTempSReal = np.std(dfTempReal, axis=0)

        # Make a list of the means, where the length of the list is NrOfCelllines
        dfListMean.append(dfTempM)
        dfListStd.append(dfTempS)
        dfListMeanReal.append(dfTempMReal)
        dfListStdReal.append(dfTempSReal)

    # Bind the mean lists as rows to form a matrix with dimensions NrOfCellines x (NrOfDoses x NrOfTimepoints)
    dataM = np.vstack(dfListMean)
    dataS = np.vstack(dfListStd)
    dataMReal = np.vstack(dfListMeanReal)
    dataSReal = np.vstack(dfListStdReal)

    # Transpose the matrices to form a matrix with dimensions (NrOfDoses x NrOfTimepoints) x NrOfCelllines
    dataMT = np.transpose(dataM)
    dataST = np.transpose(dataS)

    # Append all the rows to form a matrix with dimensions 1 x (NrOfCelllines x NrOfDoses x NrOfTimepoints)
    if len(doseList) > 0:
        dataMTR = dataMT.reshape(len(obsList) * nrOfTps * len(doseList), )
        dataSTR = dataST.reshape(len(obsList) * nrOfTps * len(doseList), )
    elif len(doseList) == 0:
        dataMTR = dataMT.reshape(len(obsList) * nrOfTps, )
        dataSTR = dataST.reshape(len(obsList) * nrOfTps, )

    dataNoise = [dataMTR, dataSTR]

    # Add option for MSE
    dfRepsM = np.ones(
        (np.max([ps.shape[0] for ps in dfReps]), len(dfReps) * dfReps[0].shape[1])) * np.nan  # define empty array
    for i, mat in enumerate(dfReps):
        # populate columns
        dfRepsM[0: mat.shape[0], i * mat.shape[1]: i * mat.shape[1] + mat.shape[1]] = mat

    # Make data array per replicate
    for row in range(0, dfRepsM.shape[0]):
        # Get first replicate
        dfRepTemp = dfRepsM[row, :]
        # Reshape into same format as dataM
        dfRepTemp = dfRepTemp.reshape(len(obsList), nrOfTps * len(doseList))
        # Transpose
        dfRepTempT = np.transpose(dfRepTemp)
        # Append all the rows to form a matrix with dimensions 1 x (NrOfCelllines x NrOfDoses x NrOfTimepoints)
        if len(doseList) > 0:
            dfRepTempTR = dfRepTempT.reshape(len(obsList) * nrOfTps * len(doseList), )
        elif len(doseList) == 0:
            dfRepTempTR = dfRepTempT.reshape(len(obsList) * nrOfTps, )

        dataR.append(dfRepTempTR)


    ###################################################
    # Part V - Save the settings based on the ODE file,
    # to be used in other scripts.
    ###################################################

    # Check whether stressName is defined
    if not "stressName" in globals():
        sys.exit("Error: <stressName> parameter is not defined in %s. Choose a <stressName> or define as empty list." % odeFile)

    # If stressName is missing in the odeFile, then the data input cannot be used for fitting.
    # Checking whether, if there is perturbation, the stressName is defined.
    if not doseColName is None:
        if len(stressName) == 0:
            sys.exit("Error: <stressName> parameter is not defined, but since <doseColName> is defined, there should be a state variable indicated that represents stress.")

    # TODO: add option for non-explicit stress input.
    # Is the stress represented as ODE or as explicit numeric function?
    if stressName[0] in stateList:
        warnings.warn("<stressName> is the same as a state variable. Script continues with the assumption that you have no explicit stress function.")
        stressExplicit = False
    else:
        stressExplicit = True

    # Is the stress input fixed?
    if stressExplicit:
        if (len(PKparms) + len(doseParmsList)) == 0 and len(doseList) > 0:
            stressInput = "fixed"
            print("Stress input parameters are fixed.")
        elif (len(PKparms) + len(doseParmsList)) > 0 and len(doseList) > 0:
            stressInput = "free"
            print("Stress input parameters are not fixed.")
        elif len(doseList) == 0:
            stressInput = None
            print("There is no explicit stress input function, but only ODEs. ")
            if not len(stressName) == 0:
                warnings.warn("There is no information on concentrations, but <stressName> in %s is still defined. Setting <stressName> to an empty list." % odeFile)
                stressName = []
            if len(doseNames) > 0:
                warnings.warn("There is no information on concentrations, but <doseNames> in %s are still defined. Setting <doseNames> to an empty list." % odeFile)
                doseNames = []
    elif not stressExplicit:
        if len(doseList) > 0:
            stressInput = "multidose"
        elif len(doseList) == 0:
            # TODO check whether this is right!
            stressInput = None
            print("There is no stress input. ")
            if not len(stressName) == 0:
                warnings.warn("There is no information on concentrations, but <stressName> in %s is still defined. Setting <stressName> to an empty list." % odeFile)
                stressName = []
            if len(doseNames) > 0:
                warnings.warn("There is no information on concentrations, but <doseNames> in %s are still defined. Setting <doseNames> to an empty list." % odeFile)
                doseNames = []

    ###################################################
    # Part VI - Modify stress input to suit the stress
    # type and origin and create parameter lists
    ###################################################

    # Check whether the doseList and doseNames are of equal length
    if not len(doseNames) == len(doseList):
        sys.exit("List of the doses as denoted in the script or derived from the data (doseList) and doseNames (in ODE system file of %s) are not of equal length." % (odeFile))

    # Get the index of the ODE equation that belongs to the stress name
    bool_list = [state in stressName for state in stateList]
    idx_true = [i for i, x in enumerate(bool_list) if x]

    # Get original stateOList and doseParmsList
    stateOList_orig = stateOList
    stateList_orig = stateList
    doseParmsList_orig = doseParmsList
    Nse_orig = len(stateOList_orig)
    Ns_orig = len(stateList_orig)
    Nd_orig = len(doseParmsList_orig)

    if not stressExplicit:
        if stressInput == "multidose":
            # 1. Make new state variables list <stateList>
            stateListNew1 = [i for i in doseNames] # [stressName[0] + "_" + i for i in doseNames]
            stateListNew2 = [state for state in stateList if state not in set(stressName)]
            stateListNew = stateListNew1 + stateListNew2
            print("New state variable list stateList is {}".format(stateListNew))

            # 2. Make new list of free initial states <stateOList>
            stateOListNew = []
            stressinistate = [stressinistate for stressinistate in stateOList if stressName[0] in stressinistate]
            if len(stressinistate) == 0:
                sys.exit("Error: Original name for initial stress input cannot be found.")
            stressinistate = stressinistate[0]
            for inistate in stateOList:  # ['E2_init', 'ER_init']
                if stressinistate == inistate:
                    for i in doseParmsList:
                        # Replace E2 with E2_EC1
                        inistatenew = inistate.replace(stressName[0], i) # stressName[0] + "_" + i)
                        stateOListNew.append(inistatenew)
                else:
                    stateOListNew.append(inistate)
            print("New list of free state variables is {}".format(stateOListNew))

            # 3. Make new list with all initial states <allIniStatesList>
            allIniStatesListNew1 = [stressinistate.replace(stressName[0],i) for i in doseNames]
            allIniStatesListNew2 = [state for state in allIniStatesList if state not in set([stressinistate])]
            allIniStatesListNew = allIniStatesListNew1 + allIniStatesListNew2
            print("New state variable list allIniStatesListNew is {}".format(allIniStatesListNew))

            # 4. Make new list for known state names <stateOKnownNameList>
            # Get the known doseNames
            knownDoseParms = [doseparm for doseparm in doseNames if doseparm not in set(doseParmsList)]
            # Add E2_init and replace E2 with E2_EC1
            if not stressinistate in stateOKnownNameList:
                stateOKnownNameListNew = []
                [stateOKnownNameListNew.append(knownState) for knownState in stateOKnownNameList]
                # Replace E2 with E2_EC1
                for i in knownDoseParms:
                    stressinistateNew = stressinistate.replace(stressName[0], i) # stressName[0] + "_" + i)
                    stateOKnownNameListNew.append(stressinistateNew)
            else:
                stateOKnownNameListNew = []
                for knownstate in stateOKnownNameList:
                    if knownstate == stressinistate:
                        for i in knownDoseParms:
                            stressinistateNew = knownstate.replace(stressName[0], i) #stressName[0] + "_" + i)
                            stateOKnownNameListNew.append(stressinistateNew)
                    else:
                        stateOKnownNameListNew.append(knownstate)
            # Order the stateOKnownNameListNew
            indices = {c: i for i, c in enumerate(allIniStatesListNew)}
            stateOKnownNameListNew = sorted(stateOKnownNameListNew, key=indices.get)
            # Print output
            print("New list of known initial states is {}".format(stateOKnownNameListNew))

            # 5. Make new list for known states <stateOKnownList>
            stateOKnownListNewString = "stateOKnownListNew = " + str(stateOKnownNameListNew).replace("\'","").replace("\"","")
            print("New list of known initial states is {}".format(stateOKnownListNewString))
            # Create initial dosing regime variable
            # EC1 = 0 needs to become EC1_init = 0
            for knownstate in stateOKnownNameListNew:
                for knowndose in knownDoseParms:
                    if knowndose in knownstate:
                        sText = knownstate + " = " + knowndose
                        print(sText)
                        exec(sText)
            exec(stateOKnownListNewString)

            # Update lists
            stateList = stateListNew
            stateOList = stateOListNew
            stateOKnownList = stateOKnownListNew
            stateOKnownNameList = stateOKnownNameListNew
            allIniStatesList = allIniStatesListNew

    # # TODO original code
    # # Parameter list
    # paraList = paraStarList + doseParmsList + PKparms + paraSnowList
    # # List of free parameters
    # paraFreeList = stateOList + doseParmsList + PKparms + paraSnowList

    # Parameter list
    if stressExplicit:
        paraList = paraStarList + doseParmsList + PKparms + paraSnowList
        # List of free parameters
        paraFreeList = stateOList + doseParmsList + PKparms + paraSnowList
    else:
        # TODO check this!!
        paraList = paraStarList + paraSnowList
        # List of free parameters
        paraFreeList = stateOList + paraSnowList
        # Set doseParmaList to an empty list
        doseParmsList = []

    # Make a list of all the lists of parameter names
    paraNamesList = [allIniStatesList, stateOKnownList, stateOKnownNameList, stateList, stateOList, obsList,
                     paraStarList, paraSO, paraSnowList, paraList, paraFreeList, doseParmsList, PKparms]

    # Number of states
    Ns = len(stateList)
    # Number of initial states to estimate
    Nse = len(stateOList)
    # Number of doses to estimate
    Nd = len(doseParmsList)
    # Number of PK parameters to estimate
    Npk = len(PKparms)
    # Number of parameters and initial states to estimate
    Np = len(paraFreeList)
    # Number of star parameters
    Nstar = len(paraStarList)
    # Number of snow parameters
    Nsnow = len(paraSnowList)
    # Number of scaling/offset parameters
    Nso = len(paraSO)
    # Number of observations Nm
    Nm = len(obsList)
    # Number of time points
    Nt = nrOfTps
    # End time point
    tend = nrOfTps
    # Time points; np.linspace(start, stop, num)
    tspan = np.linspace(1, tend, Nt)

    # Check whether there are lower or upper bounds
    if not ("ubounds" in globals() and "lbounds" in globals()):
        bounds = ([-np.inf] * (Nse + Nd + Npk + Nsnow), [np.log(1000)] * (Nse + Nd + Npk + Nsnow))
    elif ("ubounds" in globals() and "lbounds" in globals()):
        # Make a base list with lower and upper bounds from -Inf to Inf
        baseList_lb = [-np.inf] * (Nse + Nd + Npk + Nsnow)
        baseList_ub = [np.log(1000)] * (Nse + Nd + Npk + Nsnow)

        if not stressExplicit:
            # Change the index names in the lbounds an ubounds
            # lbounds = ([1,np.log(10)],[2,np.log(0.5)],[3,-np.inf])
            # Get the original indices for the doses
            oldid = list(np.arange(Nse_orig, Nse_orig + Nd_orig))
            newid = []
            for i,dosename in enumerate(doseParmsList_orig):
                # Get the new index for the doses in doseParmsList_orig
                newid = newid + [j for j,stateO in enumerate(stateOList) if dosename in stateO]

            # Get the original indices for the states
            # TODO change old indices for states to not include the stress
            oldis = []
            for i, stateO_orig in enumerate(stateOList_orig):
                if stressName[0] in stateO_orig:
                    pass
                elif stateO_orig in stateOList:
                    oldis.append(i)
            newis = []
            for i,stateO_orig in enumerate(stateOList_orig):
                # Get the new index for the states in stateOList_orig
                newis = newis + [j for j, stateO in enumerate(stateOList) if stateO_orig in stateO]

            # Replace the indices in the lower bounds
            j = 0
            k = 0
            for ilb, lb in enumerate(lbounds):
                if lb[0] in oldid:
                    lbounds[ilb][0] = newid[j]
                    j += 1
                elif lb[0] in oldis:
                    lbounds[ilb][0] = newis[k]
                    k += 1
                else:
                    lbounds[ilb][0] = lbounds[ilb][0]-1

           # Replace the indices in the upper bounds
            j = 0
            k = 0
            for iub, ub in enumerate(ubounds):
                if ub[0] in oldid:
                    ubounds[iub][0] = newid[j]
                    j += 1
                elif ub[0] in oldis:
                    ubounds[iub][0] = newis[k]
                    k += 1
                else:
                    ubounds[iub][0] = ubounds[iub][0]-1

            # Also change the indices for the priors
            j = 0
            k = 0
            for ip, prior in enumerate(priors):
                if prior[0] in oldid:
                    priors[ip][0] = newid[j]
                    j += 1
                elif prior[0] in oldis:
                    priors[ip][0] = newis[k]
                    k += 1
                else:
                    priors[ip][0] = priors[ip][0] - 1
            # Change the order of priors
            priors = list(priors)
            priors.sort()
            priors = tuple(priors)

        # Change the base list according to the changes in the txt file
        for lb in lbounds:
            baseList_lb[lb[0]] = lb[1]
        for ub in ubounds:
            baseList_ub[ub[0]] = ub[1]

        # Make one bounds tuple with the lower and upper bounds
        bounds = (baseList_lb, baseList_ub)

        print("Bounds seem to be correctly defined.")
    else:
        sys.exit(
            "Oops! Upper bounds and lower bounds are not well defined. The variable names should be 'lbounds' and 'ubounds'")

    Settings = ["testRun = " + str(testRun),
                "stressExplicit = " + str(stressExplicit),
                "stressInput = \"" + stressInput + "\"",
                "plateID_list = " + str(plateID_list),
                "realDataColName = \"" + realDataColName + "\"",
                "interpolDataColName = \"" + interpolDataColName + "\"",
                "interpolTimeColName = \"" + interpolTimeColName + "\"",
                "realTimeColName = \"" + realTimeColName + "\"",
                "palette = " + str(palette),
                "doseList = " + str(doseList),
                "doseFunctionType = \"" + doseFunctionType + "\"",
                "markers = " + str(markers),
                "paraList = " + str(paraList),
                "paraFreeList = " + str(paraFreeList),
                "paraNamesList = " + str(paraNamesList),
                "Ns = " + str(Ns),
                "Nse = " + str(Nse),
                "Nd = " + str(Nd),
                "Npk = " + str(Npk),
                "Np = " + str(Np),
                "Nstar = " + str(Nstar),
                "Nsnow = " + str(Nsnow),
                "Nso = " + str(Nso),
                "Nm = " + str(Nm),
                "Nt = " + str(Nt),
                "tend = " + str(tend),
                "tspan = " + str(list(tspan)),
                "doseList = " + str(doseList),
                "stateOList_orig = " + str(stateOList_orig),
                "stateList_orig = " + str(stateList_orig),
                "doseParmsList_orig = " + str(doseParmsList_orig),
                "Nse_orig = " + str(Nse_orig),
                "Ns_orig = " + str(Ns_orig),
                "Nd_orig = " + str(Nd_orig)]

    # Make symbols of the variables
    symbs = ["tsym","dose"] + doseParmsList + stateList + stateOList + obsList + paraStarList + paraSnowList + paraList + PKparms + stressName
    sText = makeSymbols(symbs)
    exec(sText)
    print("Executed: %s" % sText)
    Settings = Settings + [sText]

    # Check whether priors are in globals
    if not "priors" in globals():
        priors = ()

    # Initialize the input function as the 'heaviside function' (with the hill-based formula)
    if stressExplicit:
        if stressInput == "fixed":
            sText = "exposurePiecewise = Piecewise("
            for i, doseName in enumerate(doseNames):
                sText = sText + "(" + doseName + " * exp(-tsym * tau1), (dose >= " + str(doseList[i] - 0.00001) + ") & (dose <= " + str(doseList[i] + 0.00001) + ")),"
            sText = sText + "(0, True))"
            exec(sText)
            print("Executed: %s" % sText)
            Settings = Settings + [sText]
            sText = "Stress = exposurePiecewise * Piecewise((1, (tsym > 0)), (0,True))"
            exec(sText)
            print("Executed: %s" % sText)
            Settings = Settings + [sText]
        elif stressInput is None:
            print("No stress function made")
        elif stressInput == "free":
            # Make ut and utautowrap
            sText = makeExposurePWText(doseNames, doseList, funcType = doseFunctionType)
            exec(sText)
            print("Executed: %s" % sText)
            Settings = Settings + [sText]
            sText = makeutText(PWfunctionName ="exposurePiecewise", tName = "tsym")
            exec(sText)
            print("Executed: %s" % sText)
            Settings = Settings + [sText]
    elif not stressExplicit:
        if stressInput == "multidose":
            # Insert Piecewise function in ODE equation.
            # Make Stress functions
            sText = "Stress = Piecewise("
            for i, doseName in enumerate(doseNames):
                sText = sText + "(" + doseName + ", (dose >= " + str(
                    doseList[i] - 0.00001) + ") & (dose <= " + str(doseList[i] + 0.00001) + ")),"
            sText = sText + "(0, True))"
            # 'Stress = Piecewise((EC1, (dose >= -1e-05) & (dose <= 1e-05)),(EC2, (dose >= 4.99999) & (dose <= 5.00001)),(EC3, (dose >= 9.99999) & (dose <= 10.00001)),(0, True))'
            exec(sText)
            print("Executed: %s" % sText)
            Settings = Settings + [sText]
            # sText = "Stress = exposurePiecewise * Piecewise((1, (tsym > 0)), (0,True))"
            # # 'Stress = exposurePiecewise * Piecewise((1, (tsym > 0)), (0,True))'
            # exec(sText)
            # print("Executed: %s" % sText)
            # Settings = Settings + [sText]


    ###################################################
    # Part VII - Make starting parameter sets with
    # Latin Hypercube Sampling
    ###################################################

    # Make initial parameter values
    pInitList = list()
    if NiniP >= 1:
        if stressExplicit:
            parmsList = lhs(Nse + Nd + Npk + Nsnow, samples=NiniP)
        else:
            parmsList = lhs(Nse + Npk + Nsnow, samples=NiniP)

        for parmset in range(0, len(parmsList)):
            pInit = parmsList[parmset]

            # Change parameter values according to the priors
            for prior in priors:
                pInit[prior[0]] = prior[1]

            # Append parameters to list
            pInitList.append(pInit)
        # print(pInitList)
    else:
        sys.exit(
            "Error: Number of initial parameter sets should be 1 or an integer larger than 1. Please correct the value of \'NiniP\'")

    ###################################################
    # Part VIII - Read the ODEs
    ###################################################

   # Read in the ODEs from a file
    # Define the system of ODEs either by reading a file or explicitly in this script
    # Read in the ODEs from a file
    try:
        fh = open(pathToODEs, 'r')
        print("Path to text file with ODEs exists.")
        with open(pathToODEs) as f:
            lines = f.readlines()
        startReading = False
        stopReading = False
        # Make a dictionary that will contain the ODEs as specified in the file
        odeDict = {}
        obsDict = {}
        print("The system of ODES:")
        k = 0
        for i, line in enumerate(lines):
            if "####START_ODES####\n" == lines[i - 1]:
                startReading = True
            elif "####END_ODES####\n" == line:
                stopReading = True

            if startReading and not stopReading:
                line_orig = line
                # exec(line)

                # Replace state variable for stress with exposure piecewise with stress names
                if not stressExplicit:
                    print("Replacing state variable %s for Stress..." % stressName[0])
                    # Replace all E2s with Stress
                    line = line.replace(stressName[0], "Stress")

                line_repl = line

                # Make a dictionary that contains the ODEs
                # It is a match if the first chars of line are f0, f1, ..., fn
                match = re.search('^f[0-9]+', line)
                if match:
                    parts = line.split(" = ")
                    # print(parts)
                    # print(odeDict)
                    # print("---")
                    # Multiply the f ODE function that contains the stress
                    if not stressExplicit and (len(doseNames) > 1) and (int(parts[0][1]) in idx_true):
                        for j in range(0,len(doseNames)):
                            # Change the odeDict
                            odeDict[parts[0][0]+str(j)] = parts[1][:-1]
                            # Change the line
                            line = line_repl.replace(parts[0], parts[0][0]+str(j))
                            exec(line)
                            print("Executed: %s" % line)
                            k += 1
                    else: #if not stressExplicit and (len(doseNames) > 1) and not (int(parts[0][1]) in idx_true):
                        # Change the dir
                        odeDict[parts[0][0] + str(k)] = parts[1][:-1]
                        # Change the line
                        line = line_repl.replace(parts[0], parts[0][0] + str(k))
                        exec(line)
                        print("Executed: %s" % line)
                        k += 1
                    # else:
                    #     odeDict[parts[0][0]+str(k)] = parts[1][:-1]
                    #     k += 1

                # Make a dictionary that contains the observables
                # It is a match if the first chars of line are g0, g1, ..., gn
                match = re.search('^g[0-9]+', line)
                if match:
                    parts = line.split(" = ")
                    obsDict[parts[0]] = parts[1][:-1]

                # Make new lines
                print("Original line: %s" % (line_orig))
                print("New line: %s " % line)
                exec(line)
                print("Executed: %s" % line)

    except FileNotFoundError:
        sys.exit("Error: Path to text file with ODEs does not exist.")


    ###################################################
    # Part IX - Make a README file and initiate a file
    # containing finished runs
    ###################################################

    # Make a README.txt file
    f = open(path + "README.txt", "w+")
    f.write("This folder contains the output of the parameter estimation for model %s on %s. \n" % (MODEL, DATETIME))
    f.write("All output is saved in %s. \n" % path)
    f.write("%s: %s \n\n" % (MODEL, OPTIONAL_TEXT))
    f.write("Stress input was %s. \n" % (str(stressInput)))
    f.write("The data used for parameter estimation was in file %s. \n" % dataFile)
    f.write("The ode functions used for parameter estimation were in file %s. \n" % odeFile)
    f.write("The number of initial parameter sets was %s run on %s cores. \n" % (NiniP, Ncpu))
    f.write("The function describing the stress over time was of type %s .\n" % doseFunctionType)
    f.write("The concentrations used as starting point were %s. \n" % (', '.join([str(i) for i in doseList])))
    f.write("The free parameters that were estimated are %s. \n" % (', '.join(paraFreeList)))
    f.write("The star parameters that were calculated from the estimates are %s. \n" % (', '.join(paraStarList)))
    f.close()

    # Initialise the file containing a cost matrix
    pathToCSVfile = path + DATETIME + "_MH_" + MODEL + "Model_FinishedRuns_costs.csv"
    with open(pathToCSVfile, "a") as g:
        # fcntl.flock(g, fcntl.LOCK_EX)
        g.write("," + "Parmset" + "," + "Cost" + "\n")
        # fcntl.flock(g, fcntl.LOCK_UN)
        g.close()

    ###################################################
    # Part X - Create matrices and autowrap functions
    ###################################################

    if stressExplicit:
        if not stressInput is None:
            # Make an autowrap function for Stress
            sText = makeAutowrapText(attach="Stress", backend="cython", arguments=paraList + ["tsym", "dose"])
            exec(sText)
            print("Executed: %s" % sText)
            Settings = Settings + [sText]

            # Replace stressName in f functions with Stress
            sTextList = makeReplaceSText(stressName[0], attach="f", iRange=range(0, Ns))
            for sText in sTextList:
                exec(sText)
                print("Executed: %s" % sText)
        else:
            Stressautowrap = None
    else:
        # TODO check this!
        Stressautowrap = None

    # Make f and fautowrap
    # f a matrix of f0, ..., fn
    sText = makeMatrixText(attach="f", iRange=range(0, Ns))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = makeAutowrapText(attach="f", backend="cython", arguments=stateList + paraList + ["tsym", "dose"])
    exec(sText)
    print("Executed: %s" % sText)

    # Make fR
    # fR a matrix of fR1, ..., fn
    sText = makeMatrixText(attach="fR", iRange=range(0, Nstar), add1=True)
    exec(sText)
    print("Executed: %s" % sText)

    # Make g and gautowrap
    # g a matrix of g0, ..., gn
    sText = makeMatrixText(attach="g", iRange=range(0, Nm), add1=False)
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = makeAutowrapText(attach="g", backend="cython", arguments=stateList + paraList + ["tsym", "dose"])
    exec(sText)
    print("Executed: %s" % sText)

    # Make fRSnowi
    # Partial derivative of fR to Snow parameters
    txtList = makeDerivText(derivName="Snow", attach="fR", iRange=range(0, Nsnow), List=paraSnowList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make fRXOi
    # Partial derivative of fR to the unknown initial states
    txtList = makeDerivText(derivName="XO", attach="fR", iRange=range(0, Nse), List=stateOList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make fXi
    # Partial derivative of f to the state variables
    txtList = makeDerivText(derivName="X", attach="f", iRange=range(0, Ns), List=stateList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make fSnowi
    # Partial derivative of f to Snow parameters
    txtList = makeDerivText(derivName="Snow", attach="f", iRange=range(0, Nsnow), List=paraSnowList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make fUi
    # Partial derivative of f to dose and pk parameters
    txtList = makeDerivText(derivName="U", attach="f", iRange=range(0, Nd+Npk), List=doseParmsList+PKparms)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make fStari
    # Partial derivative of f to Star parameters
    txtList = makeDerivText(derivName="Star", attach="f", iRange=range(0, Nstar), List=paraStarList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make gXi
    # Partial derivative of g to state variables
    txtList = makeDerivText(derivName="X", attach="g", iRange=range(0, Ns), List=stateList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make gSnowi
    # Partial derivative of g to Snow parameters
    txtList = makeDerivText(derivName="Snow", attach="g", iRange=range(0, Nsnow), List=paraSnowList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make gUi
    # Partial derivative of g to dose and pk parameters
    txtList = makeDerivText(derivName="U", attach="g", iRange=range(0, Nd+Npk), List=doseParmsList+PKparms)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make gStari
    # Partial derivative of g to Star parameters
    txtList = makeDerivText(derivName="Star", attach="g", iRange=range(0, Nstar), List=paraStarList)
    for txt in txtList:
        exec(txt)
        print("Executed: %s" % txt)

    # Make matrices by binding rows
    sText = makeMatrixText2(attach="fX", iRange=range(0, Ns))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText=makeAutowrapText(attach="fXM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"])
    exec(sText)
    print("Executed: %s" % sText)

    sText = makeMatrixText2(attach="fSnow", iRange=range(0, Nsnow))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = makeAutowrapText(attach="fSnowM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"])
    exec(sText)
    print("Executed: %s" % sText)

    if len(range(0, Nd+Npk)) > 0:
        sText = makeMatrixText2(attach="fU", iRange=range(0, Nd+Npk))
        exec(sText)
        print("Executed: %s" % sText)
        # Do autowrap
        sText = makeAutowrapText(attach="fUM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"])
        exec(sText)
        print("Executed: %s" % sText)
    elif len(range(0, Nd+Npk)) == 0:
        fUM = Matrix([])
        # Do autowrap
        sText = makeAutowrapText(attach="fUM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"])
        exec(sText)
        print("Executed: %s" % sText)
    #else:
    #    print("No fU and fUM matrix made, because there is no explicitly defined stress input function.")

    sText = makeMatrixText2(attach="fStar", iRange=range(0, Nstar))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = (makeAutowrapText(attach="fStarM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"]))
    exec(sText)
    print("Executed: %s" % sText)

    sText = makeMatrixText2(attach="fRSnow", iRange=range(0, Nsnow))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = (makeAutowrapText(attach="fRSnowM", backend="cython", arguments=stateOList + paraList + ["tsym", "dose"]))
    exec(sText)
    print("Executed: %s" % sText)

    sText = makeMatrixText2(attach="fRXO", iRange=range(0, Nse))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = (makeAutowrapText(attach="fRXOM", backend="cython", arguments=stateOList + paraList + ["tsym", "dose"]))
    exec(sText)
    print("Executed: %s" % sText)

    sText = makeMatrixText2(attach="gX", iRange=range(0, Ns))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = (makeAutowrapText(attach="gXM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"]))
    exec(sText)
    print("Executed: %s" % sText)

    sText = makeMatrixText2(attach="gSnow", iRange=range(0, Nsnow))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = (makeAutowrapText(attach="gSnowM", backend="cython", arguments=stateList + paraList))
    exec(sText)
    print("Executed: %s" % sText)

    if len(range(0, Nd+Npk)) > 0:
        sText = makeMatrixText2(attach="gU", iRange=range(0, Nd+Npk))
        exec(sText)
        print("Executed: %s" % sText)
        # Do autowrap
        sText = (makeAutowrapText(attach="gUM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"]))
        exec(sText)
        print("Executed: %s" % sText)
    elif len(range(0, Nd + Npk)) == 0:
        gUM = Matrix([])
        # Do autowrap
        sText = makeAutowrapText(attach="gUM", backend="cython", arguments=stateList + paraList + ["tsym", "dose"])
        exec(sText)
        print("Executed: %s" % sText)
    # else:
    #     print("No gU and gUM matrix made, because there is no explicitly defined stress input function.")

    sText = makeMatrixText2(attach="gStar", iRange=range(0, Nstar))
    exec(sText)
    print("Executed: %s" % sText)
    # Do autowrap
    sText = (makeAutowrapText(attach="gStarM", backend="cython", arguments=stateList + paraList))
    exec(sText)
    print("Executed: %s" % sText)

    # Make a list of all autowrap functions and matrices
    autowrapList = [fautowrap, gautowrap, fXMautowrap, \
                    fStarMautowrap, fSnowMautowrap, fRXOMautowrap, fRSnowMautowrap, \
                    gXMautowrap, gStarMautowrap, gSnowMautowrap, fUMautowrap, gUMautowrap]
    matrixList = [f, fR, g]

    # Write all settings to a file
    with open(path + "settings.txt", "w+") as f:
        for listitem in Settings:
            f.write('%s\n' % listitem)
    print("Settings saved in %ssettings.txt" % path)

    ###################################################
    # Part XI - Perform parameter estimation
    ###################################################

    # Make a data frame of the parameters
    parmsDf = pd.DataFrame({'init_value':pInit}, index = paraFreeList)

    # Calculate the star parameters from the (initial or) steady states and other parameters
    pInitAllDf = transformParms(transformParms=paraStarList, useParms=stateOList + paraSnowList, df=parmsDf, matrix=fR)
    pInitAll = pInitAllDf["init_value"].values

    # Print an example
    print("Example of pInitAllDf:")
    print(pInitAllDf)

    # Save the current time
    t0 = time.time()

    # Arguments used to make the system of ODEs
    argP = [Nso, Nd, Nm, Ns, Nse, Np, Nstar, Nsnow, Npk, Nt, tspan, doseList, plateID_list, dataMTR, dataSTR, dataR, optim, paraNamesList, autowrapList, matrixList]

    # Create a numpy array of 1 matrix of dimensions pInit.ndim
    pInitList = np.array(pInitList)

    # Make a list of pInit and argP, in which argP is also a list of arguments
    # argPplus = [pInit, argP, testRun]
    argPplus = [pInitList, argP, testRun, t0, timeMax, path, DATETIME, MODEL, bounds]

    # Make a range of 0 to the number of ?
    index_est_all = range(0, NiniP) # i = {0}

    # Parallelize the execution of code by the specified number of CPU's
    pool = Pool(Ncpu)

    # 'itertools.izip' works in Python 2.7, but when run in Python 3 'itertools.izip' can be changed to 'zip'
    if sys.version_info[0] < 3:
        zipArgs = itertools.izip(index_est_all, itertools.repeat(argPplus, len(index_est_all)))
    elif sys.version_info[0] == 3:
        zipArgs = zip(index_est_all, itertools.repeat(argPplus, len(index_est_all)))


    # Do parameter estimation
    result_multicost = pool.map(runMultiesti, zipArgs)

    # Close the pooled task
    pool.close()
    pool.join()

    # Make an output dictionary
    dict1 = dict(estimationR=result_multicost)

    # Print the computation time
    t1 = time.time()
    estTime = t1-t0
    print("Minutes used for estimation: %.2f minutes" % (estTime/60.0))

    ###################################################
    # Part XII - Plot the output
    ###################################################

    # Loop over the output to save a list of ordered costs and the parameter sets
    parameterSets = []
    costs = []

    # Start plotting the output of the parameter estimation
    print("Finished parameter estimation. Busy plotting the output...")

    for pIniti in range(0,len(pInitList)):
        parameterSets.append(int(pIniti+1))
        print("This is parmset %i" % (pIniti+1))

        # Make a df of the initial parameters
        parmsDf = pd.DataFrame({'init_value': pInitList[pIniti]}, index=paraFreeList)
        # Calculate the star parameters from the (initial or) steady states and other parameters
        pInitAllDf = transformParms(transformParms=paraStarList, useParms=stateOList + paraSnowList, df=parmsDf, matrix=fR)

        # Append the cost
        costs.append(result_multicost[pIniti][2])

        # Comment out the following chunk to test plotting
        # Transform the output parameters back to linear scale
        outputParmsLog = result_multicost[pIniti][1]

        # TODO new code
        # Revert Log-transform of the initial states
        outputParmsLog[0:Nse] = np.exp(outputParmsLog[0:Nse])
        # Log-transform the pk, snow and scaling parameters
        outputParmsLog[Nse + Nd:len(outputParmsLog) - int(Nso / 2)] = np.exp(outputParmsLog[Nse + Nd:len(outputParmsLog) - int(Nso / 2)])
        outputParms = outputParmsLog
        # TODO old code
        # outputParms = np.exp(outputParmsLog)
        # outputParms[Ns:Ns+Nd] = outputParmsLog[Ns:Ns+Nd]
        # if Nso > 0:
        #     outputParms[-int(Nso / 2)::] = outputParmsLog[-int(Nso / 2)::]

        # # Uncomment and use to test plotting
        # outputParms = pInitList[pIniti]

        # Make a data frame of the new parameters
        newParmsDf = pd.DataFrame({'est_value':outputParms}, index = paraFreeList)

        # Calculate the star parameters from the (initial or) steady states and other parameters
        outputParmsDf = transformParms(transformParms=paraStarList, useParms=stateOList + paraSnowList, df=newParmsDf, matrix=fR)
        outputDf = pd.concat([pInitAllDf,outputParmsDf], axis = 1)
        pAllNew = outputParmsDf["est_value"].values

        # Make a plot of the data together with the simulation
        t = tspan

        # Make the new parameters global variables
        for row in outputParmsDf.iterrows():
            sText = str(row[0]) + " = " + str(float(row[1]))
            # print(sText)
            exec(sText)

        def replaceStress(d, dose, evaluate = False):
            dictExec = d
            dose = dose
            # Replace Stress with Piecewise function
            if not evaluate:
                for j, key in enumerate(dictExec.keys()):
                    dictExec[key] = dictExec[key].replace("Stress",str(Stress))
            else:
                for j, key in enumerate(dictExec.keys()):
                    dictExec[key] = dictExec[key].replace("Stress", str(eval(str(Stress))))
            return dictExec

        def replaceStateVars(stateList, d, string=""):
            dictExec = d

            # Replace the state variables by z[0], z[1], ..., z[n]
            for i, stateVar in enumerate(stateList):
                for j, key in enumerate(dictExec.keys()):
                    splitres = re.split(r"([^a-zA-Z]" + stateVar + "[^a-zA-Z])", dictExec[key])
                    subres = []
                    for res in splitres:
                        if re.search(r"((^|[^a-zA-Z])" + stateVar + "([^a-zA-Z]|$))", res):
                            subres = subres + [re.sub(stateVar, string + str(i) + "]", res)]
                        else:
                            subres = subres + [res]
                    joinedres = "".join(subres)
                    dictExec[key] = joinedres
            return dictExec

        def modelSimulation(z, t, args):
            [odeDictExec, Stressautowrap, p] = args

            if not Stressautowrap is None:
                # Stress level
                sText = "Stressautowrap("
                for i, var in enumerate(p):
                    if i < len(p) - 1:
                        sText = sText + str(var) + ","
                    else:
                        sText = sText + str(var) + ")"

                # Assign Stressautowrap to Stress parameter
                sText = stressName[0] + " = " + sText
                exec(sText)

            # Get the correct order of f equations:
            ordr = []
            for i in range(0, len(stateList)):
                ordr.append('f' + str(i))

            # Evaluate the ODEs
            output = []
            for i, key in enumerate(ordr):
                output.append(eval(odeDictExec[key]))

            # Output
            dzdt = output
            return dzdt

        # initial condition
        # Make a global variable z that contains the initial conditions, for example:
        # z = [Uo, Vo]
        sText = "z = ["
        allIniStates = allIniStatesList
        for stateVari in np.arange(0,len(allIniStates),1):
            if stateVari < len(allIniStates) - 1:
                sText = sText + allIniStates[stateVari] + ","
            else:
                sText = sText + allIniStates[stateVari] + "]"
        exec(sText)
        print("Executed: %s" % sText)

        # Replace Stress with Piecewise
        if not stressExplicit:
            odeDict = replaceStress(odeDict, dose, evaluate = False)
            # odeDict = {'f0': '- s_c1 * Stress * ER / (1 + (Stress / K_e2 ) + (ER / K_er))', 'f1': '- s_c1 * Stress * ER / (1 + (Stress / K_e2 ) + (ER / K_er)) + s_er - ER * d_er', 'f2': 's_c1 * Stress * ER / (1 + (Stress / K_e2 ) + (ER / K_er)) - d_c1 * C1'}

        # Replace the state variables in odeDict by z[0], z[1], ..., z[n]
        odeDictExec = replaceStateVars(stateList, odeDict, string = "z[")

        # Simulate the model with the estimated parameters (simulations) and transform the output
        # by using the mapping g (mappedSimu)
        if len(doseList) > 0:
            simulations = np.zeros(( len(doseList),len(t),len(stateList) ))
            mappedSimu = np.zeros(( len(doseList),len(t),len(obsList) ))
            # TODO if doseList == [], this does not work?
            for i,dose in enumerate(doseList):
                if Nstar == 0:
                    p = pAllNew[Nse:].tolist() + ["t"] + [dose]
                else:
                     p = pAllNew[-Nstar:].tolist() + pAllNew[Nse:-Nstar].tolist() + ["t"] + [dose]

                args = [odeDictExec, Stressautowrap, p]
                simu = odeint(modelSimulation, z, t, args=(args,))
                simulations[i,:,:] = simu

                # Replace the state variables in obsDict by simu[0], simu[1], ..., simu[n]
                # obsDict = {'g0': 'Stress', 'g1': 'ER', 'g2': 'C1'}
                if not stressExplicit:
                    obsDict_orig = copy.deepcopy(obsDict)
                    obsDictExecTmp = replaceStress(obsDict_orig, dose, evaluate=True)
                    obsDictExec = replaceStateVars(stateList, obsDictExecTmp, string="simu[:,")
                else:
                    obsDictExec = replaceStateVars(stateList, obsDict, string="simu[:,")

                if len(obsList) > 1:
                    sText = "mappedSimu[i,:,:] = np.transpose(np.array(["
                    for j, key in enumerate(sorted(obsDictExec.keys())):
                        if j < len(obsList) - 1:
                            sText = sText + obsDictExec[key] + ","
                        else:
                            sText = sText + obsDictExec[key] + "]))"
                    exec(sText)
                    print("Executed: %s" % sText)
                else:
                    sText = "mappedSimu[i,:,0] = " + list(obsDictExec.values())[0]
                    exec(sText)
                    print("Executed: %s" % sText)

        else:
            simulations = np.zeros((len(t), len(stateList)))
            mappedSimu = np.zeros((len(t), len(obsList)))
            # TODO if doseList == [], this does not work?
            if Nstar == 0:
                p = pAllNew[Nse:].tolist() + ["t"] + [0]
            else:
                p = pAllNew[-Nstar:].tolist() + pAllNew[Nse:-Nstar].tolist() + ["t"] + [0]

            args = [odeDictExec, Stressautowrap, p]
            simu = odeint(modelSimulation, z, t, args=(args,))
            simulations[ :, :] = simu

            if len(obsList) > 1:
                sText = "mappedSimu[:,:] = np.transpose(np.array(["
                for j, key in enumerate(sorted(obsDictExec.keys())):
                    if j < len(obsList) - 1:
                        sText = sText + obsDictExec[key] + ","
                    else:
                        sText = sText + obsDictExec[key] + "]))"
                exec(sText)
                print("Executed: %s" % sText)
            else:
                sText = "mappedSimu[:,0] = " + list(obsDictExec.values())[0]
                exec(sText)
                print("Executed: %s" % sText)

        # Check output path availability
        figPath = path + "Figures/"
        if os.path.isdir(figPath):
            print(figPath + " does already exist. Overwriting/adding output...")
        else:
            print(figPath + " does not exist. Creating a new directory and writing the output to " + figPath)
            os.mkdir(figPath)

        # Make plot titles
        if len(doseList) > 0:
            plot_titles = [(str(i) + " uM") for i in doseList]
        else:
            plot_titles = [""]

        try:
            #############################################################################################################
            # Plot 1: the simulations of the state variables together with the data mean
            #############################################################################################################

            if len(doseList) > 0:
                fig, ax = plt.subplots(Nm, len(doseList), sharex='col', sharey='row')
            else:
                fig, ax = plt.subplots(Nm, 1, sharex='col', sharey='row')

            begin = 0
            # This one is fully functional
            if (len(doseList) > 1) & (len(obsStateVarList) > 1):
                # Plotting of multiple doses (i) and state variables (j)
                for i, dose in enumerate(doseList):
                    for j, stateVar in enumerate(obsStateVarList):
                        # Get the base color for the plot of this state variable
                        base_col = colors.to_rgb(palette[j])
                        # Plot the data
                        ax[j, i].plot(realTimePoints, dataMReal[j, begin:begin + Nt], 'o', label=stateVar,
                                      color=base_col)
                        # Plot the errorbars around the data
                        ax[j, i].errorbar(realTimePoints, dataMReal[j, begin:begin + Nt], yerr=dataSReal[j, begin:begin + Nt],
                                          capsize=4, fmt="none",
                                          color=base_col)
                        # Plot the model simulation
                        ax[j, i].plot(dataTimePoints, mappedSimu[i, :, j], label=stateVar, color="#000000", linewidth=3)
                        ax[j, len(doseList) - 1].set_ylabel(stateVar, rotation=-90, fontsize=24)
                        ax[j, len(doseList) - 1].yaxis.set_label_coords(1.03, 0.5)
                        ax[j, i].tick_params(labelsize=16)
                    begin += Nt
                    ax[0, i].set_title(plot_titles[i],
                                       fontdict={'fontsize': 24})


            # TO DO: This one might not be fully functional
            elif (len(doseList) in [0,1]) & (len(obsStateVarList) > 1):
                # Plotting of one dose and multiple state variables (j)
                for j, stateVar in enumerate(obsStateVarList):
                    # Get the base color for the plot of this state variable
                    base_col = colors.to_rgb(palette[j])
                    # Plot the data
                    ax[j].plot(realTimePoints, dataMReal[j, begin:begin + Nt], 'o', label=stateVar,
                               color=base_col)
                    # Plot the errorbars around the data
                    ax[j].errorbar(realTimePoints, dataMReal[j, begin:begin + Nt], yerr=dataSReal[j, begin:begin + Nt],
                                   capsize=4, fmt="none",
                                   color=base_col)
                    # Plot the model simulation
                    if len(doseList) > 0:
                        ax[j].plot(dataTimePoints, mappedSimu[0, :, j], label=stateVar, color="#000000", linewidth=3)
                    elif len(doseList) == 0:
                        ax[j].plot(dataTimePoints, mappedSimu[:, j], label=stateVar, color="#000000", linewidth=3)
                    ax[j].set_ylabel(stateVar, rotation=-90, fontsize=24)
                    ax[j].yaxis.set_label_coords(1.03, 0.5)
                    ax[j].tick_params(labelsize=16)
                # begin += Nt
                ax[0].set_title(plot_titles[0],
                                fontdict={'fontsize': 24})


            # TO DO: This one might not be fully functional
            elif (len(doseList) > 1) & (len(obsStateVarList) == 1):
                # Plotting of multiple doses (i) and one state variable
                # Get the base color for the plot of this state variable
                base_col = colors.to_rgb(palette[0])
                for i, dose in enumerate(doseList):
                    # Plot the data
                    ax[i].plot(realTimePoints, dataMReal[0, begin:begin + Nt], 'o', label=obsStateVarList[0],
                               color=base_col)
                    # Plot the errorbars around the data
                    ax[i].errorbar(realTimePoints, dataMReal[0, begin:begin + Nt], yerr=dataSReal[0, begin:begin + Nt],
                                   capsize=4, fmt="none",
                                   color=base_col)
                    # Plot the model simulation
                    ax[i].plot(dataTimePoints, mappedSimu[i, :, 0], label=obsStateVarList[0], color="#000000", linewidth=3)
                    ax[i].tick_params(labelsize=16)
                    begin += Nt
                    ax[i].set_title(plot_titles[i],
                                    fontdict={'fontsize': 24})
                ax[len(doseList)-1].set_ylabel(obsStateVarList[0], rotation=-90, fontsize=24)
                ax[len(doseList)-1].yaxis.set_label_coords(1.03, 0.5)


            # TO DO: This one might not be fully functional
            else:
                # Plotting of one dose and one state variable
                # Get the base color for the plot of this state variable
                base_col = colors.to_rgb(palette[0])
                # Plot the data
                ax.plot(realTimePoints, dataMReal[0, begin:begin + Nt], 'o', label=obsStateVarList[0],
                           color=base_col)
                # Plot the errorbars around the data
                ax.errorbar(realTimePoints, dataMReal[0, begin:begin + Nt], yerr=dataSReal[0, begin:begin + Nt],
                               capsize=4, fmt="none",
                               color=base_col)
                # Plot the model simulation
                ax.plot(dataTimePoints, mappedSimu[0, :, 0], label=obsStateVarList[0], color="#000000", linewidth=3)
                ax.set_ylabel(obsStateVarList[0], rotation=-90, fontsize=24)
                ax.yaxis.set_label_coords(1.03, 0.5)
                ax.tick_params(labelsize=16)
                ax.set_title(plot_titles[0],
                                fontdict={'fontsize': 24})


            # add a big axes, hide frame
            fig.add_subplot(111, frameon=False)
            fig.subplots_adjust(hspace=0.2)
            fig.set_figheight(max(18,((Nm*2.5/1.5)+3)))
            fig.set_figwidth(max(10, ((len(doseList) * 4 / 2) + 2)))
            # fig.suptitle(x = 0.5, y = 0.9, t = 'Model vs. Data', fontsize=28)
            fig.text(0.02, 0.5, 'Intensity (a.u.)', va='center', rotation='vertical', fontsize=28)

            # hide tick and tick label of the big axes
            plt.grid(False)
            plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
            plt.xlabel("Time (h)", fontsize=28)
            plt.savefig(figPath + DATETIME + "_MH_" + MODEL + "Model_parmset" + str(pIniti + 1) + "_modelSimulationMeanData.pdf")
            plt.close()

            print("Plot 1 created!")

            #############################################################################################################
            # Plot 2: the simulations of the state variables together with the data per replicate
            #############################################################################################################

            if len(doseList) > 0:
                fig, ax = plt.subplots(Nm, len(doseList), sharex='col', sharey='row')
            else:
                fig, ax = plt.subplots(Nm, 1, sharex='col', sharey='row')

            # This one is fully functional.
            # 1. Case with multiple doses and multiple observables.
            # 1A. with replicates
            # 1B. without replicates
            if (len(doseList) > 1) & (len(obsStateVarList) > 1):
                # Plotting of multiple doses (i) and state variables (j)
                for i, dose in enumerate(doseList):
                    # 1. Case with multiple doses and multiple observables.
                    # 1A. with replicates
                    if len(plateID_list) > 0:
                        try:
                            print("Plot 2; Block 1A")
                            for j, stateVar in enumerate(obsStateVarList):
                                # Get the base color for the plot of this state variable
                                base_col = colors.to_rgb(palette[j])
                                # Plot the replicates
                                for k, rep in enumerate(plateID_list[j]):
                                    numOfRepl = len(plateID_list[j])

                                    # Add some factor to change the colors a bit
                                    scales = []
                                    for l in range(1, numOfRepl + 1):
                                        scales.append(l / numOfRepl)
                                    tuples = []
                                    for scale in scales:
                                        tuples.append((scale, scale, scale))

                                    # Check whether there are enough markers
                                    if k > (len(markers) - 1):
                                        warnings.warn(
                                            "There are not enough markers defined to plot the data. Restarting with choosing markers at the beginning of the marker list.")
                                        markers = markers + markers

                                    a = (data[stateVarName] == stateVar).tolist()
                                    b = (data[doseColName] == dose).tolist()
                                    c = (data[plateIDName] == rep).tolist()
                                    selection = [ai and bi and ci for ai, bi, ci in zip(a, b, c)]
                                    realDataCol = data[selection].iloc[0:Nt,].loc[:,realDataColName]

                                    ax[j, i].plot(realTimePoints, realDataCol, markers[k], label=stateVar + ' rep ' + str(rep),
                                                  color=tuple(l * r for l, r in zip(base_col, tuples[
                                                      k])))
                                    ax[j, i].tick_params(labelsize=16)

                                # Plot the model simulations
                                ax[j, i].plot(dataTimePoints, mappedSimu[i, :, j], label=stateVar, color="#000000", linewidth=3)
                                ax[j, len(doseList) - 1].set_ylabel(stateVar, rotation=-90, fontsize=24)
                                ax[j, len(doseList) - 1].yaxis.set_label_coords(1.03, 0.5)
                            ax[0, i].set_title(plot_titles[i],
                                               fontdict={'fontsize': 24})
                        except:
                            sys.exit("Error: plotting went wrong in plot 2, block 1A.")
                    # 1. Case with multiple doses and multiple observables.
                    # 1B. without replicates
                    else:
                        try:
                            print("Plot 2; Block 1B")
                            for j, stateVar in enumerate(obsStateVarList):
                                # Get the base color for the plot of this state variable
                                base_col = colors.to_rgb(palette[j])

                                # Plot the data per replicate
                                a = (data[stateVarName] == stateVar).tolist()
                                b = (data[doseColName] == dose).tolist()
                                selection = [ai and bi for ai, bi in zip(a, b)]
                                realDataCol = data[selection].iloc[0:Nt,].loc[:,realDataColName]

                                ax[j, i].plot(realTimePoints, realDataCol, markers[0],
                                              label=stateVar,
                                              color=base_col)
                                ax[j, i].tick_params(labelsize=16)

                                # Plot the model simulations
                                ax[j, i].plot(dataTimePoints, mappedSimu[i, :, j], label=stateVar, color="#000000", linewidth=3)
                                ax[j, len(doseList) - 1].set_ylabel(stateVar, rotation=-90, fontsize=24)
                                ax[j, len(doseList) - 1].yaxis.set_label_coords(1.03, 0.5)
                            ax[0, i].set_title(plot_titles[i],
                                               fontdict={'fontsize': 24})
                        except:
                            sys.exit("Error: plotting went wrong in plot 2, block 1B.")

            # 2. Case with one dose and multiple observables.
            # 2A. with replicates
            # 2B. without replicates
            elif (len(doseList) in [0,1]) & (len(obsStateVarList) > 1):
                # 2. Case with one dose and multiple observables.
                # 2A. with replicates
                # Plotting of one dose and multiple state variables (j)
                if len(plateID_list) > 0:
                    try:
                        print("Plot 2; Block 2A")
                        for j, stateVar in enumerate(obsStateVarList):
                            # Get the base color for the plot of this state variable
                            base_col = colors.to_rgb(palette[j])
                            for k, rep in enumerate(plateID_list[j]):
                                numOfRepl = len(plateID_list[j])

                                # Add some factor to change the colors a bit
                                scales = []
                                for l in range(1, numOfRepl + 1):
                                    scales.append(l / numOfRepl)
                                tuples = []
                                for scale in scales:
                                    tuples.append((scale, scale, scale))

                                # Check whether there are enough markers
                                if k > (len(markers) - 1):
                                    warnings.warn(
                                        "There are not enough markers defined to plot the data. Restarting with choosing markers at the beginning of the marker list.")
                                    markers = markers + markers

                                # Plot the data per replicate
                                if len(doseList) == 1:
                                    a = (data[stateVarName] == stateVar).tolist()
                                    b = (data[doseColName] == doseList[0]).tolist()
                                    c = (data[plateIDName] == rep).tolist()
                                    selection = [ai and bi and ci for ai, bi, ci in zip(a, b, c)]
                                    realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]
                                elif len(doseList) == 0:
                                    a = (data[stateVarName] == stateVar).tolist()
                                    b = (data[plateIDName] == rep).tolist()
                                    selection = [ai and bi for ai, bi in zip(a, b)]
                                    realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]

                                ax[j].plot(realTimePoints, realDataCol, markers[0],
                                           label=stateVar + ' rep ' + str(rep),
                                           color=tuple(l * r for l, r in zip(base_col, tuples[
                                               k])))
                                ax[j].tick_params(labelsize=16)

                            # Plot the model simulations
                            if len(doseList) == 1:
                                ax[j].plot(dataTimePoints, mappedSimu[0,:, j], label=stateVar, color="#000000", linewidth=3)
                            # TODO: I think this is correct, but needs a check
                            elif len(doseList) == 0:
                                ax[j].plot(dataTimePoints, mappedSimu[ :, j], label=stateVar, color="#000000",
                                           linewidth=3)
                            ax[j].set_ylabel(stateVar, rotation=-90, fontsize=24)
                            ax[j].yaxis.set_label_coords(1.03, 0.5)
                        ax[0].set_title(plot_titles[0],
                                           fontdict={'fontsize': 24})
                    except:
                        sys.exit("Error: plotting went wrong in plot 2, block 2A.")

                # 2. Case with one dose and multiple observables.
                # 2B. without replicates
                else:
                    try:
                        print("Plot 2; Block 2B")
                        for j, stateVar in enumerate(obsStateVarList):
                            # Get the base color for the plot of this state variable
                            base_col = colors.to_rgb(palette[j])

                            # Plot the data per replicate
                            if len(doseList) == 1:
                                a = (data[stateVarName] == stateVar).tolist()
                                b = (data[doseColName] == doseList[0]).tolist()
                                selection = [ai and bi for ai, bi in zip(a, b)]
                                realDataCol = data[selection].iloc[0:Nt,].loc[:,realDataColName]
                            elif len(doseList) == 0:
                                selection = (data[stateVarName] == stateVar).tolist()
                                realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]

                            ax[j].plot(realTimePoints, realDataCol, markers[0],
                                          label=stateVar,
                                          color=base_col)
                            ax[j].tick_params(labelsize=16)

                            # Plot the model simulations
                            ax[j].plot(dataTimePoints, mappedSimu[:, j], label=stateVar, color="#000000", linewidth=3)
                            ax[j].set_ylabel(stateVar, rotation=-90, fontsize=24)
                            ax[j].yaxis.set_label_coords(1.03, 0.5)
                        ax[0].set_title(plot_titles[0],
                                           fontdict={'fontsize': 24})
                    except:
                        sys.exit("Error: plotting went wrong in plot 2, block 2B.")

            # 3. Case with multiple doses and one observable.
            # 3A. with replicates
            # 3B. without replicates
            elif (len(doseList) > 1) & (len(obsStateVarList) == 1):
                # 3. Case with multiple doses and one observable.
                # 3A. with replicates
                if len(plateID_list) > 0:
                    try:
                        print("Plot 2; Block 3A")
                        # Get the base color for the plot of this state variable
                        base_col = colors.to_rgb(palette[0])
                        # Plotting of multiple doses (i) and one state variable
                        for i, dose in enumerate(doseList):
                            for k, rep in enumerate(plateID_list[0]):
                                numOfRepl = len(plateID_list[0])

                                # Add some factor to change the colors a bit
                                scales = []
                                for l in range(1, numOfRepl + 1):
                                    scales.append(l / numOfRepl)
                                tuples = []
                                for scale in scales:
                                    tuples.append((scale, scale, scale))

                                # Check whether there are enough markers
                                if k > (len(markers) - 1):
                                    warnings.warn(
                                        "There are not enough markers defined to plot the data. Restarting with choosing markers at the beginning of the marker list.")
                                    markers = markers + markers

                                if len(doseList) == 1:
                                    a = (data[stateVarName] == obsStateVarList[0]).tolist()
                                    b = (data[doseColName] == dose).tolist()
                                    c = (data[plateIDName] == rep).tolist()
                                    selection = [ai and bi and ci for ai, bi, ci in zip(a, b, c)]
                                    realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]
                                elif len(doseList) == 0:
                                    a = (data[stateVarName] == stateVar).tolist()
                                    b = (data[plateIDName] == rep).tolist()
                                    selection = [ai and bi for ai, bi in zip(a, b)]
                                    realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]

                                ax[i].plot(realTimePoints, realDataCol, markers[0],
                                           label=stateVar,
                                           color=base_col)
                                ax[i].tick_params(labelsize=16)

                            # Plot the model simulations
                            ax[i].plot(dataTimePoints, mappedSimu[i, :, 0], label=obsStateVarList[0], color="#000000", linewidth=3)
                            ax[i].set_title(plot_titles[i],
                                            fontdict={'fontsize': 24})
                        ax[len(doseList) - 1].set_ylabel(obsStateVarList[0], rotation=-90, fontsize=24)
                        ax[len(doseList) - 1].yaxis.set_label_coords(1.03, 0.5)
                    except:
                        sys.exit("Error: plotting went wrong in plot 2, block 3A.")

                # 3. Case with multiple doses and one observable.
                # 3B. without replicates
                else:
                    try:
                        print("Plot 2; Block 3B")
                        # Get the base color for the plot of this state variable
                        base_col = colors.to_rgb(palette[0])
                        # Plotting of multiple doses (i) and one state variable
                        for i, dose in enumerate(doseList):

                            a = (data[stateVarName] == stateVar).tolist()
                            b = (data[doseColName] == dose).tolist()
                            selection = [ai and bi for ai, bi in zip(a, b)]
                            realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]

                            ax[i].plot(realTimePoints, realDataCol, markers[0],
                                       label=stateVar,
                                       color=base_col)
                            ax[i].tick_params(labelsize=16)

                            # Plot the model simulations
                            ax[i].plot(dataTimePoints, mappedSimu[i, :, 0], label=stateVar, color="#000000", linewidth=3)
                            ax[i].set_title(plot_titles[i],
                                            fontdict={'fontsize': 24})
                        ax[len(doseList) - 1].set_ylabel(obsStateVarList[0], rotation=-90, fontsize=24)
                        ax[len(doseList) - 1].yaxis.set_label_coords(1.03, 0.5)
                    except:
                        sys.exit("Error: plotting went wrong in plot 2, block 3B.")

            # 4. Case with one dose and one observable.
            # 4A. with replicates
            # 4B. without replicates
            else:
                # Get the base color for the plot of this state variable
                base_col = colors.to_rgb(palette[0])
                # 4. Case with one dose and one observable.
                # 4A. with replicates
                if len(plateID_list) > 0:
                    try:
                        print("Plot 2; Block 4A")
                        # Plotting of one dose and one state variable
                        for k, rep in enumerate(plateID_list[0]):
                            numOfRepl = len(plateID_list[0])

                            # Add some factor to change the colors a bit
                            scales = []
                            for l in range(1, numOfRepl + 1):
                                scales.append(l / numOfRepl)
                            tuples = []
                            for scale in scales:
                                tuples.append((scale, scale, scale))

                            # Check whether there are enough markers
                            if k > (len(markers) - 1):
                                warnings.warn(
                                    "There are not enough markers defined to plot the data. Restarting with choosing markers at the beginning of the marker list.")
                                markers = markers + markers

                            if len(doseList) == 1:
                                a = (data[stateVarName] == obsStateVarList[0]).tolist()
                                b = (data[doseColName] == doseList[0]).tolist()
                                c = (data[plateIDName] == rep).tolist()
                                selection = [ai and bi and ci for ai, bi, ci in zip(a, b, c)]
                                realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]
                            elif len(doseList) == 0:
                                a = (data[stateVarName] == obsStateVarList[0]).tolist()
                                b = (data[plateIDName] == rep).tolist()
                                selection = [ai and bi for ai, bi in zip(a, b)]
                                realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]

                            ax.plot(realTimePoints, realDataCol, markers[0],
                                       label=obsStateVarList[0],
                                       color=tuple(l * r for l, r in
                                                zip(base_col, tuples[k])))
                            ax.tick_params(labelsize=16)

                        # Plot the model simulations
                        ax.plot(dataTimePoints, mappedSimu[ :, 0], label=obsStateVarList[0], color="#000000",
                                   linewidth=3)
                        ax.set_title(plot_titles[0],
                                        fontdict={'fontsize': 24})
                        ax.set_ylabel(obsStateVarList[0], rotation=-90, fontsize=24)
                        ax.yaxis.set_label_coords(1.03, 0.5)

                    except:
                        sys.exit("Error: plotting went wrong in plot 2, block 4A.")

                # 4. Case with one dose and one observable.
                # 4B. without replicates
                else:
                    try:
                        print("Plot 2; Block 4B")
                        if len(doseList) == 1:
                            a = (data[stateVarName] == obsStateVarList[0]).tolist()
                            b = (data[doseColName] == doseList[0]).tolist()
                            selection = [ai and bi for ai, bi in zip(a, b)]
                            realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]
                        elif len(doseList) == 0:
                            selection = (data[stateVarName] == obsStateVarList[0]).tolist()
                            realDataCol = data[selection].iloc[0:Nt, ].loc[:, realDataColName]

                        ax.plot(realTimePoints, realDataCol, markers[0],
                                   label=obsStateVarList[0],
                                   color=base_col)
                        ax.tick_params(labelsize=16)

                        # Plot the model simulations
                        ax.plot(dataTimePoints, mappedSimu[ :, 0], label=obsStateVarList[0], color="#000000",
                                   linewidth=3)
                        ax.set_title(plot_titles[0],
                                        fontdict={'fontsize': 24})
                        ax.set_ylabel(obsStateVarList[0], rotation=-90, fontsize=24)
                        ax.yaxis.set_label_coords(1.03, 0.5)

                    except:
                        sys.exit("Error: plotting went wrong in plot 2, block 4B.")


            # add a big axes, hide frame
            fig.add_subplot(111, frameon=False)
            fig.subplots_adjust(hspace=0.4)
            fig.set_figheight(max(18,((Nm*2.5/1.5)+3)))
            fig.set_figwidth(max(10, ((len(doseList) * 4 / 2) + 2)))
            # fig.suptitle('Model vs. Data', fontsize=28)
            fig.text(0.02, 0.5, 'Intensity (a.u.)', va='center', rotation='vertical', fontsize=28)

            # hide tick and tick label of the big axes
            plt.grid(False)
            plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
            plt.xlabel("Time (h)", fontsize=28)
            plt.savefig(figPath + DATETIME + "_MH_" + MODEL + "Model_parmset" + str(pIniti+1) + "_modelSimulationReplicateData.pdf")
            plt.close()

            print("Plot 2 created!")

            #############################################################################################################
            # Plot 3: the simulations for all state variables
            #############################################################################################################

            if not stressExplicit:
                Ns = Ns_orig
                # Get the indices for the concentrations
                bool_list = [not state in stateList_orig for state in stateList]
                idx_true = [i for i, x in enumerate(bool_list) if x]
                idx_false = [i for i, x in enumerate(bool_list) if not x]
                simunew = np.zeros([simulations.shape[0],simulations.shape[1],len(stateList_orig)])
                for i in range(0,simulations.shape[0]):
                    pos = 0
                    for j in range(0,simulations.shape[2]):
                        if j == i or j in idx_false:
                            simunew[i,:,pos] = simulations[i,:,j]
                            pos += 1
                simulations = simunew

            if len(doseList) > 1:
                fig, ax = plt.subplots(nrows = Ns, ncols = len(doseList), sharex='col', sharey='row')
            else:
                fig, ax = plt.subplots(nrows = Ns, ncols = 1, sharex='col', sharey='row')

            if (len(doseList) > 1) & (len(stateList) > 1):
                # Plotting of multiple doses (i) and state variables (j)
                for i in range(0, len(doseList)):
                    for j in range(0, Ns):
                        ax[j, i].plot(dataTimePoints, simulations[i, :, j], label=stateList_orig[j], color='#000000', linewidth=3)
                        ax[j, len(doseList) - 1].set_ylabel(stateList_orig[j], rotation=-90, fontsize=24)
                        ax[j, len(doseList) - 1].yaxis.set_label_coords(1.03, 0.5)
                        ax[j, i].tick_params(labelsize=16)
                    ax[0, i].set_title(plot_titles[i],
                                       fontdict={'fontsize': 24})

            elif (len(doseList) in [0,1]) & (len(stateList) > 1):
                if len(doseList) > 0:
                    # Plotting of one dose and multiple state variables (j)
                    for j in range(0, Ns):
                        ax[j].plot(dataTimePoints, simulations[0, :, j], label=stateList_orig[j], color='#000000', linewidth=3)
                else:
                    # Plotting of one dose and multiple state variables (j)
                    for j in range(0, Ns):
                        ax[j].plot(dataTimePoints, simulations[ :, j], label=stateList_orig[j], color='#000000', linewidth=3)
                    ax[j].set_ylabel(stateList_orig[j], rotation=-90, fontsize=24)
                    ax[j].yaxis.set_label_coords(1.03, 0.5)
                    ax[j].tick_params(labelsize=16)
                ax[0].set_title(plot_titles[0],
                                fontdict={'fontsize': 24})

            elif (len(doseList) > 1) & (len(stateList) == 1):
                # Plotting of multiple doses (i) and one state variable
                for i, dose in enumerate(doseList):
                    ax[i].plot(dataTimePoints, simulations[i, :, 0], label=stateList_orig[0], color='#000000', linewidth=3)
                    ax[i].tick_params(labelsize=16)
                    ax[i].set_title(plot_titles[i],
                                    fontdict={'fontsize': 24})
                ax[len(doseList)-1].set_ylabel(stateList[0], rotation=-90, fontsize=24)
                ax[len(doseList)-1].yaxis.set_label_coords(1.03, 0.5)

            else:
                if len(doseList) > 0:
                    # Plotting of one dose and one state variable
                    ax.plot(dataTimePoints, simulations[0, :, 0], label=stateList_orig[0], color='#000000', linewidth=3)
                else:
                    # Plotting of one dose and one state variable
                    ax.plot(dataTimePoints, simulations[ :, 0], label=stateList_orig[0], color='#000000', linewidth=3)
                ax.set_ylabel(stateList_orig[0], rotation=-90, fontsize=24)
                ax.yaxis.set_label_coords(1.03, 0.5)
                ax.tick_params(labelsize=16)
                ax.set_title(plot_titles[0],
                             fontdict={'fontsize': 24})

            plt.tight_layout()

            # add a big axes, hide frame
            fig.add_subplot(111, frameon=False)
            fig.subplots_adjust(hspace=0.2)
            fig.set_figheight(max(20,((Nm*2.5/1.5)+3)))
            if len(doseList) > 0:
                fig.set_figwidth(max(10, ((len(doseList) * 4 / 2) + 2)))
            else:
                fig.set_figwidth(max(10, ((1 * 4 / 2) + 2)))
            fig.suptitle('Model', fontsize = 28)
            fig.text(0.02, 0.5, 'Intensity (a.u.)', va='center', rotation='vertical', fontsize=28)

            # hide tick and tick label of the big axes
            plt.grid(False)
            plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
            plt.xlabel("Time (h)", fontsize=28)
            plt.savefig(figPath + DATETIME + "_MH_" + MODEL +"Model_parmset" + str(pIniti+1) + "_InnerStates.pdf")
            plt.close()

            print("Plot 3 created!")

            print("Saved the figures")
        except:
            print("Plotting failed for parmset %s" % (str(pIniti+1)))

    costdf = pd.DataFrame(np.array([[int(i) for i in parameterSets],costs]),index = ['Parameter_set','Cost'])
    costdf = costdf.transpose()
    costdf = costdf.sort_values(by = "Cost",axis = 0)
    costdf.to_csv(path + DATETIME + "_MH_" + MODEL +"Model_costs.csv")

    print("Script ended!")


