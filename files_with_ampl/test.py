import json

# Read the JSON file and load data into a dictionary
with open("data.json", "r") as file:
    data = json.load(file)

# Access and print individual elements from the loaded dictionary
print(data["solutionvect_stat"][0])

