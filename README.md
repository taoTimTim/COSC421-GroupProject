# COSC421-GroupProject

## Virtual Environment
This project uses a Python virtual environment to manage dependencies.  
Creating the virtual environment **should only be done once**, and activating it is required **each time you work on the project**.


The following are commands for windows

### Creating a virtual enviroment
python -m venv env
### Activate it
source env/Scripts/activate  # Git Bash
### OR
.\env\Scripts\activate      # PowerShell

## Installing Packages
All required packages are listed in requirements.txt
Install them with:
pip install -r requirements.txt

If you import a library that requires a new package, add it to requirements.txt

## Running Code
run in terminal using:

python src/module_name.py


## Data 
Add bdf data to the data folder but should not be commited to gihub. Notes the data folder is in the .gitignore file so don't store data outside of it.