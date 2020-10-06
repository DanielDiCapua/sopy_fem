import json

def sopy_fem_help(exampleType):
    if(exampleType == "structural_TRUSS02"):
        jsonText = structural_TRUSS02_jsonData()

    print(json.dumps(jsonText, indent=4))



def structural_TRUSS02_jsonData():
    return {
        "ProblemType": "Structural_Mechanics",
        "AnalysisType": "StaticAnalysis",
        "Mesh": {
            "ElemType": "TRUSS02",
            "Nodes": [
                {
                    "x": 0.0,
                    "y": 0.0
                },
                {
                    "x": 1.0,
                    "y": 1.0
                },
                {
                    "x": 1.0,
                    "y": 0.0
                }
            ],
            "Elements": [
                {
                    "Connectivities": [1, 2],
                    "MaterialId": 1
                },
                {
                    "Connectivities": [2, 3],
                    "MaterialId": 1
                }
            ]
        },
        "Materials": [
            {
                "Young": 20000,
                "Density": 7800,
                "Area": 0.01
            }
        ],
        "Constraints": [
            {
                "Node": 1,
                "Activation": [1, 1],
                "Values": [0.0, 0.0]
            },
            {
                "Node": 3,
                "Activation": [1, 1],
                "Values": [0.0, 0.0]
            }

        ],
        "Loads": {
            "Point_Loads": [
                {
                    "Node": 2,
                    "Values": [10.0, 0]
                }
            ]
        },
        "Postprocess": {
            "Show_displacements": True,
            "Show_deformed": True,
            "Deformed_scale": 1.2,
            "Show_forces": True,
            "Show_reactions": True
        }
}

# if __name__ == '__sopy_fem_help__':
#     sopy_fem_help()