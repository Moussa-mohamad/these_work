import subprocess
import os

def run_external_script(script_path, python_executable):
    try:
        result = subprocess.run([python_executable, script_path], capture_output=True, text=True, check=True)
        print("Output from the external script:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error occurred while running the external script:")
        print(e.stderr)

# Specify the absolute path to your external script
script_path = r"C:\Users\mmoussa\Desktop\rhino_test\2D_implementation\external_script.py"

# Specify the path to the external Python interpreter
external_python_executable = r"C:\Users\mmoussa\AppData\Local\Programs\Python\Python38\python.exe"  # Update this path to your external Python executable

# Check if the script file exists
if os.path.isfile(script_path):
    # Run the external script with the specified external Python interpreter
    run_external_script(script_path, external_python_executable)
else:
    print(f"The script file does not exist: {script_path}")
