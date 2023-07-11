#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
@author: Jorrit Bakker
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn import decomposition
from adjustText import adjust_text
import sys
sys.path.append('../Scikit_bio_modules')
import scikitbio_CCA_RDA as bio
import warnings

'''
Python module that contains functions for the purpose of performing ordination methods on input data and plotting the ordination scores and loadings.

reader : Takes a collection of dataframes and extracts the data for analysis.
null : Takes in dataframes and removes empty cells in the given manner.
transform : Takes in a dataframe and transforms the data.
PCA : Performs PCA on an input dataframe using scikit-learn package.
CCA : Performs CCA on an input dataframe using scikit-bio CCA module files
RDA : Performs RDA on an input dataframe using scikit-bio RDA module files
plot : Plots scores and loadings from an ordination method.
'''

def reader(InputFile, Units, Select_sheet, Row_name, Keep_variables, verbose = False):
    '''
    Function that takes in a dataframe sourced from an Excel file following a specific format and selects the desired columns/rows based on the input parameters

    Parameters
    ----------
    InputFile : Dataframe
        The input Excel file in the form of a dataframe. This can be the whole of the Excel file, as the code distinguishes the Environment and Species based on the rest of the input parameters.
    Units : Boolean
        Set to true when the second row contains units, those will have to be removed for the analysis. Otherwise set to False; when the second row contains the first of the observations.
    Select_sheet : List
        A list containing which sheet name(s) should be used. Multiple entries in the list means those sheets will be merged together as one dataframe.
    Row_name : String
        The name of the observation column (column containing the different observation points). In the source Excel file, this should be put somewhere into the first row.
    Keep_variables : List
        A list containing the names of the variables(columns) that should be used for the ordination.
    verbose : Boolean, optional
        Set to True to get messages in the Console about the status of the run code. The default is False.

    Returns
    -------
    Variables : Dataframe
        A dataframe containing the editted dataframe, representing either the Species or Environmental variables

    '''
    
    # Checking the amount of sheets that need to be merged together. The values is used to evaluete how much sheets remain.
    len_sheets = len(Select_sheet)
    
    #Checking if there are sheets named that are not contained in the input file. 
    ErrorName = [Sheet for Sheet in Select_sheet if Sheet not in InputFile]
    if ErrorName:
        raise Exception("It appears that the sheets %s are not contained in the given input file. Please make sure the right sheet names are entered." % ErrorName)
    
    # Preparing the data
    if verbose: print("Preparing the data...")
    # If statement checking if sheets need to be merged (if there is more than 1 name in Select_sheet),
    if len_sheets > 1:
        # For each sheet, remove the units row if present. Then set the observation column as the row names.
        for i in range(len_sheets):
            if Units:
                InputFile[Select_sheet[i]] = InputFile[Select_sheet[i]].drop('1')
            InputFile[Select_sheet[i]] = InputFile[Select_sheet[i]].set_index(Row_name)
        # Merge the first two entries in the Select_sheet list together on their row names.
        AllVariables = pd.merge(InputFile[Select_sheet[0]], InputFile[Select_sheet[1]], on = Row_name)
        # While there are still unmerged sheets remaining, take the next sheet name on the list and merge it to the dataframe.
        while len_sheets >= 3:
            i = 2
            AllVariables = pd.merge(AllVariables, InputFile[Select_sheet[i]], on=Row_name)
            i += 1
            len_sheets -= 1
    else:
        # WHen there is only one sheet, take that select sheet, drop the unit row if present and then set the observation colunn as row names.
        if Units:
            InputFile[Select_sheet[0]] = InputFile[Select_sheet[0]].drop('1')
        AllVariables = InputFile[Select_sheet[0]].set_index(Row_name)
    
    #Checking if there are variables named to use for ordination that are not contained in the selected sheet
    ErrorVars = [Var for Var in Keep_variables if Var not in AllVariables]
    if ErrorVars:
        raise Exception("It appears that the variables %s are not contained in the given input sheets. Please make sure the right sheet and variable names are entered" % ErrorVars)
    
    # Select only the variables named in the Keep_variables list.
    Variables = AllVariables[Keep_variables]
    
    # Sort the dataframe on the row names.
    Variables = Variables.sort_index()

    return(Variables)

def null(data, Drop_row = [], Set_Null = False, verbose = False):
    '''
    Function that takes in a dictionairy containing dataframes. Then removes select rows and mutates the cells containing NULL values based on the input parameters.

    Parameters
    ----------
    data : Dictionary
        A dictionary containing one or more dataframes. If the analysis method is PCA, there is one dataframe containing data and one empty list. Otherwise, it contains two dataframes.
    Drop_row : List, optional
        Parameter. List of rows that should be removed from (both) dataframes. Leave as default (empty list) if none should be removed. The default is [].
    Set_Null : String or Boolean, optional
        Parameter. Determines what to do with cells containing NULL values. Set to False to remove rows with missing values. Set to "average" to replace the missing values with the average of their corresponding variable. Set to "median" to replace the missing values with the median of their corresponding variable. Set to any numeric value to replace all empty cells with that value. The default is False.
    verbose : Boolean, optional
        Set to True to get messages in the Console about the status of the run code. The default is False.

    Returns
    -------
    data : Dictionary
        A dictionary containg one or more dataframes. If the analysis method is PCA, there is one dataframe containing data and one empty list. Otherwise, it contains two dataframes

    '''
    
    # Getting the Species dataframe from the data dictionary and removing the rows named in the Drop_row list.
    Species = data['Species']
    Nan_rows = []
    # If there is a Species dataframe:
    if isinstance(Species, pd.DataFrame):
        Species = Species.drop(Drop_row)
        # Identifying which rows and columns contain any amount of NULL cells and putting them in a list.
        Nan_rows = Species[Species.isna().any(axis=1)].index.tolist()
        Nan_cols = Species.columns[Species.isna().any()].tolist()
    
    #Getting the Environment dataframe from the data dictionary. If the method is PCA, this is instead an empty list.
    Environment = data['Environment']
    #If there is an Environment dataframe, remove the rows named in the Drop_row list. Then identify which rows and columns contain any amount of NULL cells.
    if isinstance(Environment, pd.DataFrame):
        Environment = Environment.drop(Drop_row)
        #Extend the list with rows containing NULL cells from the Species data with the rows containing NULL cells from the Environment data. 
        Nan_rows.extend(Environment[Environment.isna().any(axis=1)].index.tolist())
        Nan_cols_env = Environment.columns[Environment.isna().any()].tolist()
    
    # If there are any amount of rows containing NULL cells, the NULL values will be removed in the manner given by the Set_Null parameter.
    if Nan_rows:
        # If the entry of Set_Null is a float of interger value, fill every empty cell with that number.
        if type(Set_Null) == float or type(Set_Null) == int:
            if isinstance(Species, pd.DataFrame):
                Species = Species.fillna(Set_Null)
            if isinstance(Environment, pd.DataFrame):
                Environment = Environment.fillna(Set_Null)
            print('The values of the empty cells have been changed to the set value %s.' % Set_Null)
        
        # If the entry of Set_Null is "average", set every NULL cell will be set to the average of its corresponding variable
        elif Set_Null == "average":
            if isinstance(Species, pd.DataFrame):
                for var in Nan_cols:
                    Species[var] = Species[var].fillna(Species[var].mean(skipna = True))
            if isinstance(Environment, pd.DataFrame):
                for var in Nan_cols_env:
                    Environment[var] = Environment[var].fillna(Environment[var].mean(skipna = True))
            print('The values of the empty cells have been changed to the average of their corresponding variables.')
        
        # If the entry of Set_Null is "median", set every NULL cell will be set to the median of its corresponding variable
        elif Set_Null == "median":
            if isinstance(Species, pd.DataFrame):
                for var in Nan_cols:
                    Species[var] = Species[var].fillna(Species[var].median(skipna = True))
            if isinstance(Environment, pd.DataFrame):
                for var in Nan_cols_env:
                    Environment[var] = Environment[var].fillna(Environment[var].median(skipna = True))
            print('The values of the empty cells have been changed to the median of their corresponding variables.')
        
        # If the entry of Set_Null is none of the above, every row containing NULL cells will be removed from both the Species and Environment data. 
        else:
            if isinstance(Species, pd.DataFrame):
                Species = Species.drop(Nan_rows)
            if isinstance(Environment, pd.DataFrame):
                Environment = Environment.drop(Nan_rows)
                print('Environment has also been removed')
            print('The rows with the names %s have been removed since they contain NULL cells.' % Nan_rows)
    
    data = {'Species': Species, 'Environment': Environment}
     
    return(data)
    
def transform(dataframe, Center = False, Standardize = True, LogScale = False, LogScale_A = 1, LogScale_B = 1, verbose = False):
    '''
    A function that performs certain transformations on the input data based on the input parameters.

    Parameters
    ----------
    dataframe : Dataframe
        A dataframe containing data to be transformed.
    Center : Boolean, optional
        Parameter. When set to true, the data will be centered per variable, giving it an average of 0. The default is False.
    Standardize : Boolean, optional
        Parameter. When set to true, the data will be standardized per variable, giving it an average of 0 and a range between 0 and 1. The default is False.. The default is True.
    LogScale : Boolean or List, optional
        Parameter. Set to false to disable log transformation, set to True when wanting to log transform every variable. Otherwise provide a list of all variable column names you want to transform The default is False.
    LogScale_A : Integer or float, optional
        Parameter. Determines the method of Log transformation, where this is A in log10(Ax+B). The default is 1.
    LogScale_B : Interger or float, optional
       Parameter. Determines the method of Log transformation, where this is B in log10(Ax+B). The default is 1.
    verbose : Boolean, optional
         Set to True to get messages in the Console about the status of the run code. The default is False.

    Returns
    -------
    dataframe : Dataframe
        A dataframe containing the transformed input data.
    '''
    
    # If the LogScale parameter is in the form of a list, log transform every column of which the name is in the list.
    if type(LogScale) == list:
        if verbose: print("A selection of the data is being log scaled...")
        for Variable in LogScale:
            dataframe[Variable] = np.log10(LogScale_A * dataframe[Variable] + LogScale_B)
    # Otherwise, if the LogScale parameter is set to true, log scale every value in the dataframe.
    elif LogScale:
        if verbose: print("Data is being log scaled...")
        dataframe = np.log10(LogScale_A * dataframe + LogScale_B)
    
    # If the Center parameter is true, but not the Standardize parameter, center the data by subtracting the mean of every column from its values.
    if Center and not Standardize:
        if verbose: print("Data is being centered...")
        dataframe = dataframe.apply(lambda col: col-col.mean())
    # If the Standardize parameter is true, standardize all columns.
    elif Standardize:
        if verbose: print("Data is being standardized...")
        dataframe = stats.zscore(dataframe, axis=0)
    
    return(dataframe)

def PCA(data, Row_name, verbose = False):
    '''
    Function that performs Principal Component Analysis using sklearn.decomposition.PCA on the input data and gives the site scores and loadings.

    Parameters
    ----------
    data : Dictionary
       A dictionary containing the Species and Environment dataframes. The dataframes represents the data on which the PCA takes place. Either Environment or Species can be an empty list as long as the other is a dataframe.
    Row_name : String
       The name of the observation column (column containing the different observation points). In the source Excel file, this should be put somewhere into the first row.
    verbose : Boolean, optional
       Set to True to get messages in the Console about the status of the run code. The default is False.

    Returns
    -------
    results : Dictionary
        A dictionary containing the scores and loadings of the PCA, and the percentage of the variation explained by the first principal components.
    names : Dictionary
        A dictionary containing the names of the samples and the Environmental and Species variables
    '''
    
    Species = data["Species"]
    Environment = data["Environment"]
    Species_head = []
    Environment_head = []
    
    # If there are only Species variables, those become the data for the PCA
    if isinstance(Species, pd.DataFrame) and not isinstance(Environment, pd.DataFrame):
        Species_head = list(Species.columns.values)
        heads = Species_head
        Variables = Species
        
    # If there are only Environmental variables, those become the data for the PCA
    elif isinstance(Environment, pd.DataFrame) and not isinstance(Species, pd.DataFrame):
        Environment_head = list(Environment.columns.values)
        heads = Environment_head
        Variables = Environment
        
    # If there are both Species and Environmental variables, combine the dataframes into one for the purpose of the PCA
    elif isinstance(Environment, pd.DataFrame) and isinstance(Species, pd.DataFrame):
        Species_head = list(Species.columns.values)
        Environment_head = list(Environment.columns.values)
        heads = Species_head + Environment_head
        Variables = pd.merge(Species, Environment, on=Row_name)
    
    wells = list(Variables.index)
    
    # Checking if the dimensions of the dataframe allow for PCA
    length, width = Variables.shape
    if length < width:
        raise Exception("There are more variables than there are samples. PCA is only possible when there are at least as much samples as there are variables.")
    
    if verbose: print("Performing PCA...")
    
    # Using scikit.decomposoition.PCA with an amount of components equal to the amount of variables, then getting the loadings, scores and explained variance ratio.
    pca = decomposition.PCA(n_components=len(Variables.columns)) 
    pca.fit(Variables)
    loadings = pca.components_.T
    PCAscores = pca.transform(Variables)
    variances = pca.explained_variance_ratio_
    
    # Taking the first two PC for plotting
    loadings = loadings[:, [0,1]]
    scores = PCAscores[:, [0,1]]
    
    # Printing information about the succes of the PCA
    percent_explained = np.around(100*variances/np.sum(variances), decimals=2)
    
    # If verbose is True, will report on the correlation between the first two principal components and the explained variance by each principal component.
    if verbose:
        for i in range(len(percent_explained)):
            print('PC%s explains' % (i+1), '%s percent of the total variance' % percent_explained[i])
        print()
        coef = np.corrcoef(scores[:,0], scores[:,1])
        print('The correlation coefficient between PC1 and PC2 = %s.' % coef[0,1])
    
    results = {"loadings": loadings, "scores": scores, "percent_explained": percent_explained}
    names = {"heads": heads, "Species_head": Species_head, "Environment_head": Environment_head, "wells": wells}

    return(results, names)

def CCA(data, verbose = False):
    '''
    Function that performs Canonical Correspondence Analysis using skbio.stats.ordination.CCA on the input data and gives the site scores and loadings.

    Parameters
    ----------
    data : Dictionar
        A dictionary containing the Species and Environment dataframes. The dataframes represents the data on which the CCA takes place.
    verbose : Boolean, optional
        Set to True to get messages in the Console about the status of the run code. The default is False.

    Returns
    -------
    results : Dictionary
        A dictionary containing the scores and loadings of the RDA. The Environment variables are positionned first, then the Species variables.
    names : Dictionary
        A dictionary containing the names of the samples and the Environmental and Species variables

    '''
    
    Species = data["Species"]
    Environment = data["Environment"]
    Species_head = list(Species.columns.values)
    Environment_head = list(Environment.columns.values)
    heads = Environment_head + Species_head
    wells = list(Species.index)
    
    # Checking if the dimensions of the dataframe allow for CCA
    length_s, width_s = Species.shape
    length_e, width_e = Environment.shape
    if length_s < width_s or length_e < width_e:
        raise Exception("There are more variables than there are samples. CCA is only possible when there are at least as much samples as there are Species or Environmental variables")
    
    if verbose: print("Performing CCA...")
    # Performing CCA using the CCA function from scikit-bio.
    sci_cca = bio.cca(Species, Environment, scaling = 2)
    loadings_Species = sci_cca.features
    scores = sci_cca.samples
    loadings_Environment = sci_cca.biplot_scores
    
    # Transforming the loadings and scores to arrays and taking the first two canonical axes.
    loadings_Species = loadings_Species.to_numpy()
    loadings_Environment = loadings_Environment.to_numpy()
    loadings_Species = loadings_Species[:, [0,1]]
    loadings_Environment = loadings_Environment[:, [0,1]]
    loadings = np.append(loadings_Environment, loadings_Species, axis=0)
    scores = scores.to_numpy()
    scores = scores[:, [0,1]]
    
    results = {"loadings": loadings, "scores": scores}
    names = {"heads": heads, "Species_head": Species_head, "Environment_head": Environment_head, "wells": wells}
    
    return(results, names)

def RDA(data, verbose = False):
    '''
    Function that performs Redundancy Analysis using skbio.stats.ordination.RDA on the input data and gives the site scores and loadings.

    Parameters
    ----------
    data : Dictionary
        A dictionary containing the Species and Environment dataframes. The dataframes represents the data on which the CCA takes place.
    verbose : Boolean, optional
        Set to True to get messages in the Console about the status of the run code. The default is False.

    Returns
    -------
    results : Dictionary
        A dictionary containing the scores and loadings of the RDA. The Environment variables are positionned first, then the Species variables.
    names : Dictionary
        A dictionary containing the names of the samples and the Environmental and Species variables

    '''
    
    
    Species = data["Species"]
    Environment = data["Environment"]
    Species_head = list(Species.columns.values)
    Environment_head = list(Environment.columns.values)
    heads = Environment_head + Species_head
    wells = list(Species.index)
    
    # Checking if the dimensions of the dataframe allow for CCA
    length_s, width_s = Species.shape
    length_e, width_e = Environment.shape
    if length_s < width_s or length_e < width_e:
        raise Exception("There are more variables than there are samples. CCA is only possible when there are at least as much samples as there are Species or Environmental variables")
    
    if verbose: print("Performing RDA...")
    # Performing RDA using the RDA function from scikit-bio.
    
    sci_rda = bio.rda(Species, Environment, scaling = 2)
    loadings_Species = sci_rda.features
    scores = sci_rda.samples
    loadings_Environment = sci_rda.biplot_scores
    
    # Transforming the loadings and scores to arrays and taking the first two canonical axes.
    loadings_Species = loadings_Species.to_numpy()
    loadings_Environment = loadings_Environment.to_numpy()
    loadings_Species = loadings_Species[:, [0,1]]
    loadings_Environment = loadings_Environment[:, [0,1]]
    loadings = np.append(loadings_Environment, loadings_Species, axis=0)
    scores = scores.to_numpy()
    scores = scores[:, [0,1]]
    
    results = {"loadings": loadings, "scores": scores}
    names = {"heads": heads, "Species_head": Species_head, "Environment_head": Environment_head, "wells": wells}
    
    return(results, names)

def plot(results, names, Method, Plot_loadings = True, Plot_scores = False, Adjust_text = False, Full_extent = False, Scale_focus = "Loadings", Figsize = False, Subfigure = False, Axis_font = 12, Label_font = 12, Scoretext_font = 6, Loadingtext_font = 8, Arrow_width = 0.002, Arrowhead_width = 0.02, Arrowhead_length = 0.02, Axis1_min = False, Axis1_max = False, Axis2_min = False, Axis2_max = False, verbose = False):
    '''
    Function that takes in ordination loadings and scores and plots them in an ordination plot.

    Parameters
    ----------
    results : Dictionary
        Dictionary that contains ordination results as numpy arrays; ordination and scores.
    names : Dictionary
        A dictionary containing the names of the samples and the Environmental and Species variables
    Method : String
        The ordination method used in the analysis.
    Plot_loadings : Boolean, optional
        If set to True, will plot the ordination loadings. The default is True.
    Plot_scores : Boolean, optional
        If set to True, will plot the ordiantion scores. The default is False.
    Adjust_text : Boolean, optional
        When set to True, scores and loading labels will be adjusted to not overlap with eachother.
    Full_extent : Boolean, optional
        When set to True, the loadings and scores will be scaled to always have a loading close to 1. The default is False.
    Scale_focus : String, optional
        Can be set to either 'Loadings', 'Scores' or 'none'. If either Loadings or Scores, the other will be scaled to the extent of the other. If none, no such scaling will take place. The default is "Loadings".
    Figsize : List, optional
        Sets the dimensions of the resulting plot. Input as a list, where the first value is the size in the horizontal direction and the second value the size in the vertical direciton. The default is False.
    Subfigure : Boolean or string, optional
        When set to false, nothing happens. When set to a string, the string will be put in the top-left of the figure. The default is False.
    Axis_font : Integer, optional
        Sets the font size of the axis titles. Default is 12.
    Label_font : Integer, optional
        Sets the font size of the axis labels. Default is 12.
    Scoretext_font : Integer, optional
        Sets the font size of the text displayed next to the site scores. Default is 6.
    Loadingtext_font : Integer, optional
        Sets the font size of the text displayed next to the loadings. Default is 8.
    Arrow_width : Float, optional
        Sets the width of the vector tail. Default is 0.002.
    Arrowhead_width : Float, optional
        Sets the width of the vector head. Default is 0.02.
    Arrowhead_length : Float, optional
        Sets the width of the vector head. Default is 0.02.
    Axis1_min : Boolean, Interger or Float, optional
        Forcibly sets the minimum value of the first ordiantion axis. Values between 0 and 1 recommended. The default is False.
    Axis1_max : Boolean, Interger or Float, optional
        Forcibly sets the maximum value of the first ordiantion axis. Values between 0 and 1 recommended. The default is False.
    Axis2_min : Boolean, Interger or Float, optional
        Forcibly sets the minimum value of the second ordiantion axis. Values between 0 and 1 recommended. The default is False.
    Axis2_max : Boolean, Interger or Float, optional
        Forcibly sets the maximum value of the second ordiantion axis. Values between 0 and 1 recommended. The default is False.
    verbose : Boolean, optional
        Set to True to get messages in the Console about the status of the run code. The default is False.
    
    Returns
    -------
    None.
    
    Running the function will display the ordination plot.
    '''

    loadings = results["loadings"]
    scores = results["scores"]
    heads = names["heads"]
    Species_head = names["Species_head"]
    Environment_head = names["Environment_head"]
    wells = names["wells"]

    texts = []
    
    if verbose: print("Determining axes and scaling...")
    # A list of possible maximum and minimum axis extents.
    axis_numbers = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0]


    # Determing the largest values in the PCA scores.
    max_load = np.max(np.abs(loadings))
    if len(scores):
        max_score = np.max(np.abs(scores))
    
    # To make sure scores and loadings are between -1 and 1, they are all scaled down when one of them falls out of thay range.
    if max_load >= 1 or Full_extent:
        loadings_scaled = loadings / (max_load * 1.05)
    else:
        loadings_scaled = loadings
    
    if len(scores):
        if max_score >= 1 or Full_extent:
            scores_scaled = scores / (max_score * 1.05)
        else:
            scores_scaled = scores
    
    # When plotting both scores and loadings, scores or loadings are scaled depending on the extent of the other. Depends on the input of Scale_focus.
    if Plot_scores and Plot_loadings:
        if Scale_focus == "Scores":
            max_scaled_score = np.max(np.abs(scores_scaled))
            loadings_scaled = loadings_scaled * max_scaled_score
        elif Scale_focus == "None":
            pass
        else:
            max_scaled_load = np.max(np.abs(loadings_scaled))
            scores_scaled = scores_scaled * max_scaled_load
        
        # Append loadings and scores for the purpose of finding the axis extent neccecariy to fit all inside the plot.
        full_coords = np.append(loadings_scaled, scores_scaled, axis=0)
    
    elif Plot_loadings:
        full_coords = loadings_scaled
    else:
        full_coords = scores_scaled
    
    # Takes the manually set value of the maximum of first ordination axis if given. Otherwise, determines it based on the data.
    if Axis1_max:
        x_lim_pos = Axis1_max
    else:
        max_x = np.max(full_coords[:,0])
        x_lim_pos = np.min(np.where((axis_numbers - max_x) > 0, axis_numbers, 2))
    
    # Takes the manually set value of the minimum of first ordination axis if given. Otherwise, determines it based on the data.
    if Axis1_min:
        x_lim_neg = Axis1_min
    else:
        min_x = np.min(full_coords[:,0])
        x_lim_neg = np.max(np.where((axis_numbers - min_x) < 0, axis_numbers, -2))
    
    # Takes the manually set value of the maximum of second ordination axis if given. Otherwise, determines it based on the data.
    if Axis2_max:
        y_lim_pos = Axis2_max
    else:
        max_y = np.max(full_coords[:,1])
        y_lim_pos = np.min(np.where((axis_numbers - max_y) > 0, axis_numbers, 2))
    
    # Takes the manually set value of the minimum of second ordination axis if given. Otherwise, determines it based on the data.
    if Axis2_min:
        y_lim_neg = Axis2_min
    else:
        min_y = np.min(full_coords[:,1])
        y_lim_neg = np.max(np.where((axis_numbers - min_y) < 0, axis_numbers, -2))
    
    # Manually sets the figure dimensions if given, and sets de dpi to 300 for good figure resolution.
    if Figsize:
        plt.figure(figsize = (Figsize[0],Figsize[1]), dpi=300)
    else:
        plt.figure(dpi=300)

    # Setting the x and y axis limits with the previously determined values
    plt.xlim(x_lim_neg, x_lim_pos)
    plt.ylim(y_lim_neg, y_lim_pos)
    # Plotting lines that indicate the origin
    plt.plot([-1, 1], [0, 0], color='grey', linewidth=0.75, linestyle='--')
    plt.plot([0, 0], [-1, 1], color='grey', linewidth=0.75, linestyle='--')

    # Plotting the ordination scores by iterating over every coordinate in the scores array, if the Plot_scores parameter is set to true.
    if Plot_scores:
        if verbose: print("Plotting the scores...")
        for i, (x, y) in enumerate(scores_scaled):
            plt.scatter(x, y, color='grey', facecolor='none', edgecolor='grey')
            # Plotting the name of the scores and storing it in a list for the purpose of adjusting the position later
            tex = plt.text(x, y, wells[i], color='black', fontsize = Scoretext_font)
            texts.append(tex)

    # Plotting the ordination loadings by iterating over every coordinate in the loadings array, if the Plot_scores parameter is set to true.
    if Plot_loadings:
        if verbose: print("Plotting the loadings...")
        for i, (x, y) in enumerate(loadings_scaled):
            # Plots Environmental and Species variables with different colours and text formatting.
            if heads[i] in Environment_head:
                plt.arrow(0, 0, x, y, color='blue', width = Arrow_width, head_length = Arrowhead_length, head_width = Arrowhead_width)
                #Plotting the name of the loading and storing it in a list for the purpose of adjusting the position later
                tex = plt.text(x, y, heads[i], color='black', fontstyle='italic', fontsize = Loadingtext_font)
            else:
                plt.arrow(0, 0, x, y, color='red', width = Arrow_width, head_length = Arrowhead_length, head_width = Arrowhead_width)
                #Plotting the name of the loading and storing it in a list for the purpose of adjusting the position later
                tex = plt.text(x, y, heads[i], color='black', weight="bold", fontsize = Loadingtext_font)
            texts.append(tex)
    
    # Plots the axis title with the percentage of variation explained if the method was PCA
    if Method == "PCA" and len(results) > 2:
        percent_explained = results["percent_explained"]
        plt.xlabel('PC1 (%s%%)' % percent_explained[0], fontsize = Axis_font)
        plt.ylabel('PC2 (%s%%)' % percent_explained[1], fontsize = Axis_font)
    # Otherwise just plots a general axis title.
    else:
        plt.xlabel('Ordination axis 1', fontsize = Axis_font)
        plt.ylabel('Ordination axis 2', fontsize = Axis_font)
   
    # Changing the axis ticks and removing all but the first and last tick labels
    ax = plt.gca()
    x_ticks = np.around(ax.get_xticks(), decimals=1)
    y_ticks = np.around(ax.get_yticks(), decimals=1)

    x_tick_labels = [str(x_ticks[0])] + [''] * (len(x_ticks) - 2) + [str(x_ticks[-1])]
    y_tick_labels = [str(y_ticks[0])] + [''] * (len(y_ticks) - 2) + [str(y_ticks[-1])]

    ax.set_xticklabels(x_tick_labels, fontsize = Label_font)
    ax.set_yticklabels(y_tick_labels, fontsize = Label_font)
    if Subfigure:
        plt.text(-0.175,1, Subfigure, transform = ax.transAxes, fontsize = 16)
    
    if Adjust_text:
        if verbose: print("Adjusting text position...")
        #adjust_text gives an irrelevant warning for the purpose of this scipt. The statement below will make sure that warning does not show up.
        warnings.filterwarnings("ignore")
        adjust_text(texts)
