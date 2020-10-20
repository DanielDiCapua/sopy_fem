import json
import os
import sys
import pkg_resources

exampleTypesSet = {"dynamics_TRUSS02",
                    "electrical_BR02",
                    "mechanics_BR02",
                    "mechanics_QU04",
                    "mechanics_TR03",
                    "structural_TRUSS02",
                    "thermal_BR02"
                }

def sopy_fem_help(exampleType="", basePath="", outputFile=""):
    if (exampleType in exampleTypesSet):
        if(basePath != ""):
            fileName = basePath + exampleType + "/data.json"
        else:
            resource_name = "Examples/" + exampleType + "/data.json"
            fileName = pkg_resources.resource_filename('sopy_fem', resource_name)
        with open(fileName, "r") as exampleFile:
            jsonText = json.load(exampleFile)
            jsonOuput = json.dumps(jsonText, indent=4)

        if(outputFile != ""):
            with open(outputFile, 'w') as f:
                print(jsonOuput, file=f)
        else:
            print(jsonOuput)
    else:
        print("Please choose one of the following examples types:\n")
        for exampleType in exampleTypesSet:
            print("  - ", exampleType)

        if (basePath == ""):
            print("\n")
            print("The example file can be sent to the console or to an output file:\n")
            print("  -For a console output use sopy_fem_help(<exampleType>) e.g.=> sopy_fem_help('mechanics_BR02') \n")
            print("  -For a file output use sopy_fem_help(<exampleType>,outputfile=<jasonfileName>) e.g.=>  sopy_fem_help('mechanics_BR02',outputFile='mechanics_bars.json') \n")


if __name__ == '__main__':
    args = sys.argv[1:]
    exampleType = ""
    basePath= "Examples/"
    if (len(args)  != 0):
        exampleType = args[0]
    sopy_fem_help(exampleType, basePath)
