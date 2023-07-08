# General
Python/R scripts made and data used for bachelor thesis aimed at develloping a diagnostic statistical tool for the purpose of interpreting aspects of bioremediation. 

## Structure
The project follows the structure below:

- `Sample_data/` -> Map containing the data used in the study to generate ordination plots and a template file.
  + `Figure4a_input_data.xlsx` -> Excel file containing the bioremediation data from `Raw_input_data.xlsx` used to make figure 4a from the thesis.
  + `Raw_input_data.xlsx` -> Excel file containing bioremediation data received for this project. No alterations were made.
  + `Restructured_input_data.xlsx` -> Excel file containing restructured bioremediation data. Restructuring follows the form of `Template_sheet.xlsx` to allow for it to be red by the Python tool.
  + `Template_sheet.xlsx` -> Excel file that shows the format that bioremediation data needs to have in order to be properly red by the Python tool.
- `Sample_figures` -> Map containing Python scripts used for generating specific plots from the thesis. Each Main_script here uses `Source_code/Main_code/analysis.py` as module. Serves as examples on how different parameters should be entered in the Main_script.
  + `Main_script_Figure4a` -> Script used to make figure 4a from the thesis. Input data is `Sample_data/Figure4a_input_data.xlsx`.
  + `Main_script_Figure5a` -> Script used to make figure 5a from the thesis. Input data is `Sample_data/Restructured_input_data.xlsx`.
  + `Main_script_Figure6f` -> Script used to make figure 6f from the thesis. Input data is `Sample_data/Restructured_input_data.xlsx`.
  + `Main_script_Figure7c` -> Script used to make figure 7c from the thesis. Input data is `Sample_data/Restructured_input_data.xlsx`.
- `Source_code/` -> Map containing all code used in the study.
  + `Main_code/` -> Map containing the Python main script and module used for generating bioremediation ordination plots.
    + `Main_script_template.py` -> Main script that calls on `analysis.py` to perform PCA, CCA or RDA ordiantion methods.
    + `analysis.py` -> Python module containing various functions needed to perform ordination.
  + `Scikit_bio_modules/` -> Map containing Python scripts that together perform RDA as implemented in the scikit-bio Python package. This provides an alternative method for when the scikit-bio package fails to install on the used Python IDE. Contains a seperate requirements file for the additional packages needed to get the RDA class running.
  + `R_validation` -> Python and R files used to validate the ordination methods implemented in the Python tool (in devellopment)
- `README.md` -> The file you are currently reading. Provides overview over the contents of this repository.
- `requirements` -> A .txt file denoting the required Python packages, see the 'requirements' section below.

## Requirements
The `requirements` file lists all the packages needed to be installed to get the main code running. The Scikit-bio module has a seperate requirements file in its map. 

It can be installed with the following method:
```bash
pip install -r requirements.txt
```
