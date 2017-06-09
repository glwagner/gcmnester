import os
import numpy as np

def replace_namelist_param(namepath, param, saveorig=True):
    
    """ Find and replace a parameter in the MITgcm namelist
    specified by datapath. The input param is a dictionary whose
    key gives the parameter and whose value gives the value.
    This program automatically determines and formats the value
    according to MITgcm's namelist requirements."""

    # Characters that identify the end of a parameter specification
    endchars = ['#', '&', '/'] 
    validtypes = [bool, str, int, float, list]
    linewidth = 8

    paramname = param.keys()[0]
    paramval = param.values()[0]

    # Convert numpy array to list 
    if type(paramval) is np.ndarray:
        if len(paramval.shape) > 1:
            raise ValueError("The numpy array that specifies the "
                    "parameter value can only have one dimension.")
        paramval = paramval.tolist()

    if not os.path.isfile(namepath):
        raise ValueError("The file {} does not exist!".format(namepath))

    if type(paramval) not in validtypes:
        raise ValueError("The parameter type is not valid")

    # Read the file into a list of strings
    with open(namepath, 'r') as namefile:
        namelist = namefile.readlines()

    # Find the first line where the parameter is specified
    linenum = 0
    for line in namelist:
        line = line.lstrip()
        if len(line) > 0:
            if line[0] not in endchars and paramname in line and '=' in line:
                if 'paramline' in locals():
                    raise RuntimeError( 
                        "Two lines have been found that contain the parameter name, \n"
                        "the equals sign '=', and do not appear to be a comment. \n"
                        "Check the input data file and remove the duplicate or "
                        "otherwise offensive line.")
                else:
                    paramline = linenum

        linenum += 1 

    # For multiline parameters, find the last line where the parameter is specified
    (linenum, endfound) = (paramline, False)
    while not endfound:
        linenum += 1
        # Determine if the end has been found
        line = namelist[linenum].lstrip()
        if len(line) == 0: 
            endfound = True
        elif line[0] in endchars or '=' in line:
            endfound = True
        elif linenum == len(namelist):
            raise RuntimeError("The end of the apparently multi-line parameter\n"
                    "specification could not be found.")
    paramendline = linenum

    # Save original file
    with open(namepath + '_orig', 'w') as namefile:
        for line in range(len(namelist)):
            namefile.write(namelist[line])

    # Write new file
    with open(namepath, 'w') as namefile:

        # Write the unmodified beginning lines
        for line in range(paramline):
            namefile.write(namelist[line])

        # Generate the new parameter lines
        if type(paramval) is str:
            newtext = " {} = '{}',\n".format(paramname, paramval)

        elif type(paramval) is bool:
            if paramval is False: paramtext = '.FALSE.'
            else: paramtext = '.TRUE.'
            newtext = " {} = {},".format(paramname, paramtext)

        elif type(paramval) in [int, float]:
            newtext = " {} = {},\n".format(paramname, paramval)

        elif type(paramval) is list:
            lineindent = 4 + len(paramname)
            newtext = " {} = ".format(paramname)
            for i in range(len(paramval)):
                newtext = newtext + "{:6}, ".format(paramval[i])
                if (i+1) % linewidth == 0 and (i+1) != len(paramval):
                    newtext = newtext + "\n{}".format(' '.ljust(lineindent))
            newtext = newtext + "\n"    

        # Write the new text into file
        print("The text replacing the parameter {} is".format(paramname)
            + "\n\n{}".format(newtext)
        )
        namefile.write(newtext)

        # Write the rest of the file
        for line in range(paramendline, len(namelist)):
            namefile.write(namelist[line])


def write_data_diagnostics(fields, freqs, levels, savedir='.', overwrite=False):

    """ Write a diagnostics file for MITgcm diagnostics.
    Each diagnostic is defined by a list of field names, frequency
    of output, vertical levels in the case of 3D diagnostics, and a 
    filename. The mandatory inputs fields, freqs, and levels are
    dictionaries whose keys correspond to the filename of the diagnostic
    and whose values are:

        fields[diagname] : List of fields in the diagnostic.

        freqs[diagname]  : Number (float) giving the frequency at which 
                           the diagnostic is outputted.

        levels[diagname] : For 3D diagnostics, the vertical levels from which
                           the diag will be extracted.
                
    The optional parameters savedir and overwrite give the directory to which
    the data.diagnostics file will be saved, and specify whether an existing
    file should be overwritten or not, respectively."""

    # Parameters
    levelLinewidth = 10
    fieldLinewidth = 3
    diagFilename = "{}/data.diagnostics".format(savedir)
    nDiags = len(fields.keys())

    # Check to see if file exists
    if os.path.isfile(diagFilename):
        if overwrite:
            os.remove(diagFilename)
        else:
            raise ValueError("File {} exists! "
                    "Either delete it or set overwrite=True.".format(diagFilename))

    (i, body) = (0, str())
    for diag in fields.keys():

        i += 1
        nFields = len(fields[diag])

        fileText  = " filename({:d}) = 'diags/{}',\n".format(i, diag)
        freqText  = " frequency({:d}) = {:.2f}\n".format(i, freqs[diag])

        # Construct text block for field data
        fieldText = " fields(1:{:d},{:d}) = ".format(nFields, i)
        for j in range(nFields):
            fieldText = fieldText + "'{}', ".format(fields[diag][j])
            if (j+1) % fieldLinewidth == 0 and (j+1) != nFields:
                fieldText = fieldText + "\n{:16}".format(' ')

        fieldText = fieldText + "\n"

        # Construct text block for level data
        if levels[diag] is None: 
            levelText = ''
        else:
            nLevels = len(levels[diag])
            levelText = " levels(1:{:d},{:d}) = ".format(nLevels, i)
            for j in range(nLevels):
                levelText = levelText + "{:3d}., ".format(levels[diag][j])
                if (j+1) % levelLinewidth == 0 and (j+1) != nLevels:
                    levelText = levelText + "\n{:18}".format(' ')
            levelText = levelText + "\n"

        body = body + fileText + freqText + fieldText + levelText

        if i != nDiags:
            body = body + "#\n#\n"

    header = (
        "# Diagnostic Package Choices\n"
        "#-----------------\n"
        "# For each output-stream:\n"
        "#  filename(n)      : prefix of the output file name (only 8.c long) for outp.stream n\n"
        "#  frequency(n):< 0 : write snap-shot output every multiple of |frequency| (iter)\n"
        "#               > 0 : write time-average output every multiple of frequency (iter)\n"
        "#  levels(:,n)      : list of levels to write to file (Notes: declared as REAL)\n"
        "#                     when this entry is missing, select all common levels of this list\n"
        "#  fields(:,n)      : list of diagnostics fields (8.c) (see 'available_diagnostics' file \n"
        "#                     for the list of all available diag. in this particular config)\n"
        "#--------------------------------------------------------------------\n"
        "#\n"
        " &diagnostics_list\n"
        "#\n"
        " dumpatlast = .TRUE.,\n"
        "#\n"
    )

    # Write the file
    diagFile = open("{}/data.diagnostics".format(savedir), 'w')
    diagFile.write(header + body)
    diagFile.close()
