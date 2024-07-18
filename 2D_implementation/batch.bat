@echo off
REM Specify the path to the specific Python version
set PYTHON_PATH=C:\Users\mmoussa\AppData\Local\Programs\Python\Python38\python.exe

REM Specify the path to the Python script
set SCRIPT_PATH=C:\Users\mmoussa\Desktop\rhino_test\2D_implementation\external_script.py

REM Run the Python script using the specified Python version
"%PYTHON_PATH%" "%SCRIPT_PATH%"

REM Open a new command prompt window, run the Python script, and keep the command line open
start cmd /k "%PYTHON_PATH% %SCRIPT_PATH%"



