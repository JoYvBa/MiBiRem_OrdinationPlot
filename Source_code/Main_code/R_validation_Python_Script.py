#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 17:19:24 2023

@author: jorrit
"""

import sys
sys.path.append('../Main_code')
import analysis as ana
import numpy as np
import pandas as pd

#%%

########################
### INPUT PARAMETERS ###
########################

# Ordination method parameters
Method = 'PCA' #Choose from 'PCA', 'CCA' and 'RDA'

#Set verbose to True to get messages about the script status in the console when running the code.
verbose = True

# Input file parameters
File_path = "../../Sample_data/Restructured_input_data.xlsx"
Units = False #Set to true if the row beneath the column names contains units.
Row_name = "well" #Put the name of the observation points here 
Species_sheet = ["DNA field"] #List with names of the sheets containing the species (dependant) variables.
Species_vars = ["Total bacteria 16SRrna", "Benzene carboxylase", "NirS", "NarG", "BssA_SRB", "BssA nitraat", "Peptococcus"]  #List of names from the species sheet to use for the analysis.
Environment_sheet = ["environ_variables"] #List with names of the sheets containing the Environmental variables.
Environment_vars = ["Sum GC", "B", "T", "E", "pm-X", "Indene", "N", "Redox", "Mn(III)", "Fe(II)", "NO3-", "SO42-", "B/T*100", "T/B*100"] #List of names from the environment sheet to use for the analysis.
Drop_row = [] #List of names of the observation points to disregard for the analysis. (Will be removed from both the Species and Environment)
Set_Null = 0 #Ordination methods require all cells to be filled. Set to False to remove rows with missing values. Set to "average" to replace the missing values with the average of their corresponding variable. Set to "median" to replace the missing values with the median of their corresponding variable. Set to any numeric value to replace all empty cells with that value.

# Transform parameters
Center = False #When set to true, centers the data; giving it an average of 0 
Standardize = True #When set to true, centers data and divides by standard deviation
LogScale_Species = False #Set to false to disable log transformation, set to True when wanting to log transform all data. Otherwise provide a list of all variable column names you want to transform
LogScale_Environment = False # Same as line above but for environment variables.
#Log transform is in the form of log10(Ax + B). Below, provide the desired A and B
LogScale_A = 1
LogScale_B = 1

# Plotting parameters
Plot_loadings = True # When set to True, the loadings will be plotted in the resulting ordination plot. At least one between Plot_loadings and Plot_scores needs to be true.
Plot_scores = False # When set to True, the scores will be plotted in the resulting ordination plot.
Full_extent = False # When set to True, the loadings and scores will be scaled to always have a loading close to 1. 
Scale_focus = "Loadings" # Can be set to either 'Loadings', 'Scores' or 'none'. If either Loadings or Scores, the other will be scaled to the extent of the other. If none, no such scaling will take place.

# The plotting parameters below are optional to change plot lay-out to your liking. They do not influence the results of the ordination.
Figsize = [6, 7] # Sets the dimensions of the resulting plot
Adjust_text = True # When set to True, scores and loading labels will be adjusted to not overlap with eachother. This may increase the execution time of the code, especially when both Plot_scores and Plot_loadings are set to True.
Subfigure = "a" # Is False by default. Can set to a string containing a letter to include it in the top left corner of the figure.
Axis_font = 16 # Font size of the axis titles
Label_font = 16 # Font size of the axis labels
Scoretext_font = 9 # Font size of the text displayed next to the site scores.
Loadingtext_font = 13 # Font size of the text displayed next to the loadings.
Arrow_width = 0.002 # Width of the vector tail 
Arrowhead_width = 0.02 # Width of the vector head
Arrowhead_length = 0.02 # length of the vector head
# Code will adjust the plot axis extent automatically to the loadings and scores. Only set these parameters to a value if you want specific axis extents.
Axis1_min = False # Forcibly sets the minimum value of the first ordiantion axis. Values between 0 and 1 recommended.
Axis1_max = False # Forcibly sets the maximum value of the first ordiantion axis. Values between 0 and 1 recommended.
Axis2_min = False # Forcibly sets the minimum value of the second ordiantion axis. Values between 0 and 1 recommended.
Axis2_max = False # Forcibly sets the maximum value of the second ordiantion axis. Values between 0 and 1 recommended.

#########################
### START MAIN SCRIPT ###
#########################

# Testing for faulty entries in parameters
MethodOptions = ['PCA', 'CCA', 'RDA']

if Method not in MethodOptions:
    raise Exception("Please select a valid method option, this can be one of the following: %s" % MethodOptions)

if not Species_sheet and not Environment_sheet:
    raise Exception("Both the Species_sheet and Environment_sheet parameter are empty. Please provide at least one sheet in one of the parameters.")

if (Method == "CCA" or Method == "RDA") and not Environment_sheet:
    raise Exception("A contrained ordination method is selected, but no Environmental sheet is given. Please enter the sheet name of the Environmental sheet in the Environment_sheet parameter or use an unconstrained ordination method.")

if (Method == "CCA" or Method == "RDA") and not Species_sheet:
    raise Exception("A contrained ordination method is selected, but no Species sheet is given. Please enter the sheet name of the Species sheet in the Species_sheet parameter or use an unconstrained ordination method.")

if not Plot_loadings and not Plot_scores:
    raise Exception("Right now, the loadings nor the scores will be plotted. Please set at least one of those to true.")

#Executing the analysis
InputFile = pd.read_excel(File_path, sheet_name = None)

# Preparing the Species dataframe if a sheet is given
if Species_sheet:
    Species = ana.reader(InputFile = InputFile, Units = Units, Select_sheet = Species_sheet, Row_name = Row_name, Keep_variables = Species_vars, verbose = verbose)
# Otherwise, Species will be an empty list for the rest of the script.
else:
    Species = []

# Preparing the Environment dataframe if a sheet is given
if Environment_sheet:
    Environment = ana.reader(InputFile = InputFile, Units = Units, Select_sheet = Environment_sheet, Row_name = Row_name, Keep_variables = Environment_vars, verbose = verbose)
# Otherwise, Environment will be an empty list for the rest of the script.
else:
    Environment = []

# The Species and Environment datasheets are put into a dictionary so they can be passed to and from functions more easily.
data = {"Species": Species, "Environment": Environment}

# Use the null function to remove any empty cells in the data.
data = ana.null(data = data, Drop_row = Drop_row, Set_Null = Set_Null, verbose = verbose)

# If there is a Species dataframe, transform the data.
if isinstance(data["Species"], pd.DataFrame):
    data["Species"] = ana.transform(dataframe = data["Species"], Center = Center, Standardize = Standardize, LogScale = LogScale_Species, LogScale_A = LogScale_A, LogScale_B = LogScale_B, verbose = verbose)

# If there is an Environment dataframe, transform the data.
if isinstance(data["Environment"], pd.DataFrame):
    data["Environment"] = ana.transform(dataframe = data["Environment"], Center = Center, Standardize = Standardize, LogScale = LogScale_Environment, LogScale_A = LogScale_A, LogScale_B = LogScale_B, verbose = verbose)

Species = data["Species"]
Environment = data["Environment"]

# Aquiring the row and column names
Species_head = list(Species.columns.values)
Environment_head = list(Environment.columns.values)
heads = Species_head + Environment_head
wells = list(Species.index)

# Saving the data to different csv files depending on which ordination method is used.
if Method == "PCA":
    pcavars = pd.merge(Species, Environment, on=Row_name)
    pcavars.to_csv('../../Sample_data/PCA_validation.csv')
elif Method == "CCA": 
    Species.to_csv('../../Sample_data/CCA_validation_species.csv')
    Environment.to_csv('../../Sample_data/CCA_validation_environment.csv')
elif Method == "RDA":
    Species.to_csv('../../Sample_data/RDA_validation_species.csv')
    Environment.to_csv('../../Sample_data/RDA_validation_environment.csv')

#%% Reading in the output of the R scripts

if Method == "PCA":
    loadings = pd.read_excel('../../Sample_data/PCA_validation_loadings.xlsx')
    loadings = loadings.to_numpy()
else:
    if Method == "CCA":
        loadings_spec = pd.read_excel('../../Sample_data/CCA_validation_species_loadings.xlsx')
        loadings_env = pd.read_excel('../../Sample_data/CCA_validation_environment_loadings.xlsx')
    elif Method == "RDA":
        loadings_spec = pd.read_excel('../../Sample_data/RDA_validation_species_loadings.xlsx')
        loadings_env = pd.read_excel('../../Sample_data/RDA_validation_environment_loadings.xlsx')
    
    # As the variable "N" is removed from canonical ordination in R due to co-linearity, it has to be removed from the list with variable names too.
    Environment_head.remove("N")
    heads.remove("N")
    
    loadings_spec = loadings_spec.to_numpy()
    loadings_env = loadings_env.to_numpy()
    loadings = np.append(loadings_spec, loadings_env, axis=0)
    
results = {"loadings": loadings, "scores": []}
names = {"heads": heads, "Species_head": Species_head, "Environment_head": Environment_head, "wells": wells}

# Plotting the data 
ana.plot(results = results, names = names, Method = Method, Plot_loadings = Plot_loadings, Plot_scores = Plot_scores, Adjust_text = Adjust_text, Full_extent = Full_extent, Scale_focus = Scale_focus, Figsize = Figsize, Subfigure = Subfigure, Axis_font = Axis_font, Label_font = Label_font, Scoretext_font = Scoretext_font, Loadingtext_font = Loadingtext_font, Arrow_width = Arrow_width, Arrowhead_width = Arrowhead_width, Arrowhead_length = Arrowhead_length, Axis1_min = Axis1_min , Axis1_max = Axis1_max, Axis2_min = Axis2_min, Axis2_max = Axis2_max, verbose = verbose)
print("Done! :)")
