import json
import os
import globalvars

def read_data():
    
    print("Input the name of the json data file with its relative path (e.g. project/example.json).")
    fileName = input("Press enter when the file is called data.json and is located in the current folder: ")
    if(fileName == ""):
        base_dir = os.getcwd()
        fileName = os.path.join(base_dir, "data.json")
    globalvars.dataFileName = fileName
    with open(fileName, "r") as read_file:
        globalvars.data = json.load(read_file)

