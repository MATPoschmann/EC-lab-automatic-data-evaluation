# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 16:59:33 2022

@author: poschmann
"""
import pandas as pd
import eclabfiles as ecf
from galvani import BioLogic
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob
import os
import itertools
import sys

#entry of filepath
#os.chdir(input('Enter full path of data folder: '))
# function to get user input in commandline. 
# Valueinputtext is the text printed for the user to ask for entry. 
# Valuetype is the data typa this function returns to the programm. 
# The while-loop and Try-function checks for the input beeing a number (int or float)
# If the input is no number there are two cases:
# Input is "e": the program terminates
# everything else: the program prints 'wrong input' and reruns the question for input
def getvalue(valueinputtext, valuetype): 
    valueinput = ''
    while (type(valueinput) != int) and (type(valueinput) != float):
        print("'e' for exit")
        valueinput = input(valueinputtext)
        try: 
            valueinput = valuetype(valueinput)
        except:
            if valueinput == 'e':
                sys.exit()
            else:
                print('wrong input')
    return valuetype(valueinput)

# plotandfit function is used here to plot single tafelplots, their tafelfit and annotated tafelslope value as text
# the function takes the values for x-axis and y-axis as lists or dataframes of same length as X and Y, respectively.
# fit_a is the slope of the fitted function
# fit_b is the y-value at x = 0
# fitparam is the found tafelslope
# subfolder is the name of subfolder the plot should be stored in
# title is the filename ending (incl. ".png") which is placed after the name of the datafile 
def plotandfit(X, Y, fit_a, fit_b, fitparam, subfolder, title):
    plt.figure(figsize=(10,8))
    plt.scatter(X, Y*1000)
    plt.plot(X, (fit_a*X+fit_b)*1000, 'r')
    #axis labels
    plt.xlabel('log(j)')
    plt.ylabel(r'$eta$ /mV')
    #legend
    plt.text(X.min(),Y.max()*1000,fitparam)
    #plotsaving
    new_filename = str(filename) + str(title)
    new_file_path = os.path.join(path, subfolder, new_filename)
    plt.savefig(new_file_path, dpi = 300)
    #if show != 0:
    #    plt.show()
    plt.close()
    return
# plotsave function is used to plot any x-y-scatterplots
# the function takes the values for x-axis and y-axis as lists or dataframes of same length as X and Y, respectively.
# subfolder is the name of subfolder the plot should be stored in
# title is the filename ending (incl. ".png") which is placed after the name of the datafile 
# xlabel is the x-axis-label as string
# ylabel is the y-axislabel as string
def plotsave(X, Y, subfolder, title, xlabel, ylabel):
    plt.figure(figsize=(10,8))
    plt.scatter(X, Y)
    #axis labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plotsaving
    new_filename = str(filename) + str(title)
    new_file_path = os.path.join(path, subfolder, new_filename)
    plt.savefig(new_file_path, dpi = 300)
    plt.close()
    return

# plotofloops function is used to plot any values of loop as scatter vs. applied potential
# datasheet is a DataFrame containing a column called 'Loop' and having columns called 'x_columnname' and 'y_columnname'
# the function plos the y_columnname vs. the x_columnname for each unique value in the loop column 
# subfolder is the name of subfolder the plot should be stored in
# title is the filename ending (incl. ".png") which is placed after the name of the datafile 
# xlabel is the x-axis-label as string
# ylabel is the y-axislabel as string
def plotofloops(datasheet, title, subfolder, x_columnname, y_columnname, x_label, y_label):
    plt.figure(figsize=(10,8))
    legend = []
    Ncolors = len(datasheet['Loop'].unique())
    colormap = plt.cm.viridis
    Ncolors = min(colormap.N,Ncolors)
    mapcolor = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]
    for loop, color in zip(datasheet['Loop'].unique(), itertools.product(mapcolor)):   
        legend.append(loop)
        plt.scatter(datasheet.loc[datasheet['Loop']==loop].loc[:, x_columnname], datasheet.loc[datasheet['Loop']==loop].loc[:,y_columnname]*1000, color = color)
    #axis labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    #legend
    norm = mpl.colors.Normalize(vmin=1, vmax=len(legend))
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colormap), label='Loop', ax=plt.gca())
    #plotsaving
    new_filename = str(filename) + str(title)
    new_file_path = os.path.join(path, subfolder, new_filename)     
    plt.savefig(new_file_path, dpi = 300)
    plt.close()
    return

# plotofloopsCV function is used to plot any potential values as scatter vs. the loopnumber
# datasheet is a DataFrame containing a column called 'Loop' and having columns called 'x_columnname' and 'y_columnname'
# the function plos the y_columnname vs. the x_columnname for each unique value in the loop column 
# subfolder is the name of subfolder the plot should be stored in
# title is the filename ending (incl. ".png") which is placed after the name of the datafile 
# xlabel is the x-axis-label as string
# ylabel is the y-axislabel as string
def plotofloopsCV(datasheet, title, subfolder, x_columnname, y_columnname, x_label, y_label):
    plt.figure(figsize=(10,8))
    legend = []
    Ncolors = len(datasheet['Loop'].unique())
    colormap = plt.cm.viridis
    Ncolors = min(colormap.N,Ncolors)
    mapcolor = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]
    for loop, color in zip(datasheet['Loop'].unique(), itertools.product(mapcolor)):   
        legend.append(loop)
        plt.scatter(datasheet.loc[datasheet['Loop']==loop].loc[:, x_columnname]*1000, datasheet.loc[datasheet['Loop']==loop].loc[:,y_columnname], color = color)
    #axis labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    #legend
    norm = mpl.colors.Normalize(vmin=1, vmax=len(legend))
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colormap), label='Loop', ax=plt.gca())
    #plotsaving
    new_filename = str(filename) + str(title)
    new_file_path = os.path.join(path, subfolder, new_filename)     
    plt.savefig(new_file_path, dpi = 300)
    plt.close()
    return


# Eatoverloop is a special plot-function that plots the EiR_at_I values determined from Modular Potentio vs. the loop
# first the funtion checks of all values are present by trying to make a plot.
# if datasets are missing it prints out 'No Values to plot for {Dataset name}'
# datasheet is a DataFrame containing a column called 'Loop' and columns called 'EiR_at_5mA/cm2', 'EiR_at_10mA/cm2', 'EiR_at_50mA/cm2', 'EiR_at_100mA/cm2'
# subfolder is the name of subfolder the plot should be stored in
# title is the filename ending (incl. ".png") which is placed after the name of the datafile 
def Eatoverloop(datasheet, subfolder, title):
    plt.figure(figsize=(10,8))
    for EAT in [r'$eta$_at_5mA/cm2', r'$eta$_at_10mA/cm2', r'$eta$_at_50mA/cm2', r'$eta$_at_100mA/cm2']:
        try:
            plt.scatter(datasheet['Loop'], (datasheet[EAT])*1000)
        except:
            print('No Values to plot for ' + str(EAT))
    #axis labels
    plt.xlabel('Loop')
    plt.ylabel(r'$\eta$ /mV')
    #legend
    plt.legend([r'$eta$_at_5mA/cm2', r'$eta$_at_10mA/cm2', r'$eta$_at_50mA/cm2', r'$eta$_at_100mA/cm2'])
    #plotsaving
    new_filename = str(filename) + str(title)
    new_file_path = os.path.join(path, subfolder, new_filename)     
    plt.savefig(new_file_path, dpi = 300)
    datasheet.to_csv(str(new_file_path) + '.txt')
    plt.close()


# getpointsnexttoj function finds the next two values neighbouring the j_value 
# tafelplottable is sorted on the 'j/mA/cm^2' column by absolute value of difference to j_value   
# after finding these values, which are ideally left and right of the j_value it makes a linear fit and returns the y_axis_value at the x_position of j_value (Eatvalue)   
# tafelplottable needs to be a DataFrame containing the column 'j/mA/cm^2' at second last position and 'EiR' at third-last position
def getpointsnexttoj(tafelplottable, j_value):
    tablesort = tafelplottable.iloc[(tafelplottable['j/mA/cm^2']- j_value).abs().argsort(), :]
    if j_value == 5:
        if (tablesort.iat[0,-2] <= j_value) & (tablesort.iat[1,-2] <= j_value):
            tablesort = tafelplottable.loc[(tafelplottable['j/mA/cm^2'] >= tablesort.iat[0,-2])]
    linear_model=np.polyfit(tablesort.iloc[0:2,-1], tablesort.iloc[0:2,-3],1)    
    #if j_value == 5:
    #    print(linear_model)
    Eatvalue = abs(linear_model[0])*np.log10(j_value)+linear_model[1]
    return Eatvalue

#class is used to process all imported datafiles resulting from PEIS-technique
class PEIS():
    def __init__(self):
        # these first lines find out at which lines in the dataframe correspond to each loop, by looking into the metadata (meta), because the data of every loop is appended to the datafile 
        self.loops = meta['loops']['indexes']
        self.Nrpointsloop = len(dataframe)-self.loops[-2]
        # following lines generate a 'Loop' column in the DataFrame using the metadata
        counter = 1
        df = []        
        for entry in self.loops[0:-1]:
            singleloop = dataframe.iloc[entry:entry+self.Nrpointsloop].copy()
            singleloop.loc[:,'Loop'] = counter
            counter += 1
            df.extend(singleloop.values.tolist())
        self.PEISdf = pd.DataFrame(df, columns = ['freq','Re(Z)','-Im(Z)','|Z|','Phase(Z)','time','<Ewe>','<I>','Cs','Cp','cycle number', 'I Range', '|Ewe|', '|I|', 'Ns', '(Q-Qo)', '<Ece>', '|Ece|','Phase(Zce)','|Zce|','Re(Zce)','-Im(Zce)','Phase(Zwe-ce)','|Zwe-ce|','Re(Zwe-ce)','-Im(Zwe-ce)','uts','Loop'])
        # REIMdf is just PEISdf reduced to 'freq', 'Loop', 'Re(Z)' and '-Im(Z)'
        self.REIMdf = self.PEISdf.iloc[:,[0, -1, 1, 2]]
        #following lines get the Re(Z)-value of the '-Im(Z)' with lowest absolute value if '-Im(Z)' crosses 0
        #if all '-Im(Z)-values are positive, the lines return the Re(Z) value of the datapoint with lowest '-Im(Z)
        df = []
        for entry in self.REIMdf.Loop.unique():  
            Im = self.REIMdf.loc[self.REIMdf['-Im(Z)'] == self.REIMdf[self.REIMdf['Loop'] == entry].iloc[:,-1].min()].reset_index(drop=True)
            if Im.iloc[-1, -1] >= 0:
                ReZ = self.REIMdf.loc[self.REIMdf['-Im(Z)'] == self.REIMdf[self.REIMdf['Loop'] == entry].iloc[:,-1].min()].reset_index(drop=True).iloc[-1, -2]
            else:
                ReZ = self.REIMdf.loc[self.REIMdf['-Im(Z)'].abs() == self.REIMdf[self.REIMdf['Loop'] == entry].iloc[:,-1].abs().min()].reset_index(drop=True).iloc[-1, -2]
            df.append([entry, ReZ])
        # Loop-numbers and respective Re(Z) values are put into ReZlist-DataFrame 
        self.ReZlist = pd.DataFrame(df, columns = ['Loop', 'Re(Z)'])
        self.ReZlist.Loop = self.ReZlist.Loop.astype(int)
        return

#class is used to process all imported datafiles resulting from OCV-technique
class OCV():
    def __init__(self):
        # these first lines find out at which lines in the dataframe correspond to each loop, by looking into the metadata (meta), because the data of every loop is appended to the datafile
        self.loops = meta['loops']['indexes']
        self.Nrpointsloop = len(dataframe)-self.loops[-2]
        # following lines generate a 'Loop' column in the DataFrame using the metadata
        counter = 1
        df = []
        for entry in self.loops[0:-1]:
            singleloop = dataframe.iloc[entry:entry+self.Nrpointsloop].copy()
            singleloop.loc[:,'Loop'] = counter
            counter += 1
            df.extend(singleloop.values.tolist())
        #generating new dataframe with selcted columns
        self.OCVdf = pd.DataFrame(df, columns = ['time','Ewe','Ece','mode','error', 'uts', 'Loop'])
        # reducing the number of columns to relevant
        self.OCVcutdf = self.OCVdf.iloc[:,[1, 2, 3, -1]]
        # if Re(Z) list from PEIS is present then Re(Z) values are added as additional column 
        try:
            correction = ReZlist.set_index('Loop').to_dict()
            for value in self.OCVdf['Loop']:
                self.OCVdf.loc[self.OCVdf['Loop'] == value, 'ReZ'] = correction[value]
        except:
            return
        return

#class is used to process all imported datafiles resulting from ZIR-technique
class ZIR():
    def __init__(self):
        # these first lines find out at which lines in the dataframe correspond to each loop, by looking into the metadata (meta), because the data of every loop is appended to the datafile
        self.loops = meta['loops']['indexes']
        self.Nrpointsloop = len(dataframe)-self.loops[-2]
        # following lines generate a 'Loop' column in the DataFrame using the metadata
        counter = 1
        df = []
        for entry in self.loops[0:-1]:
            singleloop = dataframe.iloc[entry:entry+self.Nrpointsloop].copy()
            singleloop.loc[:,'Loop'] = counter
            counter += 1
            df.extend(singleloop.values.tolist())
        # generation of new dataframe with selected columns
        self.ZIRdf = pd.DataFrame(df, columns = ['freq', 'Re(Z)', '-Im(Z)', '|Z|', 'Phase(Z)', 'time', '<Ewe>', '<I>', 'I Range', '|Ewe|', '|I|', '<Ece>', '|Ece|', 'Phase(Zce)', '|Zce|', 'Re(Zce)', '-Im(Zce)', 'Phase(Zwe-ce)', '|Zwe-ce|', 'Re(Zwe-ce)', '-Im(Zwe-ce)', 'uts', 'Loop'])
        return

#class is used to process all imported datafiles resulting from LSV-technique
class LSV():
    def __init__(self):
        # these first lines find out at which lines in the dataframe correspond to each loop, by looking into the metadata (meta), because the data of every loop is appended to the datafile
        self.loops = meta['loops']['indexes']
        self.Nrpointsloop = len(dataframe)-self.loops[-2]
        # following lines generate a 'Loop' column in the DataFrame using the metadata
        counter = 1
        df = []
        for entry in self.loops[0:-1]:
            singleloop = dataframe.iloc[entry:entry+self.Nrpointsloop].copy()
            singleloop.loc[:,'Loop'] = counter
            counter += 1
            df.extend(singleloop.values.tolist())
        #generating new dataframe with selcted columns
        self.LSVdf = pd.DataFrame(df, columns = ['time', 'control_V', 'Ewe', 'I', '(Q-Qo)', 'I Range', 'Ece', 'mode', 'ox/red', 'error', 'control changes', 'uts', 'Loop'])
        self.LSVcutdf = self.LSVdf.iloc[:,[2, 3, -1]]
        return

#class is used to process all imported datafiles resulting from CV-technique    
class CV():
    def __init__(self):
        # try making subfolder to save later results in
        try:
            os.makedirs('CVevaluation')
        except: print("path 'CVevaluation' already exists")
        # these following lines find out at which lines in the dataframe correspond to each loop, by looking into the metadata (meta), because the data of every loop is appended to the datafile
        self.loops = meta['loops']['indexes']
        self.Nrpointsloop = len(dataframe)-self.loops[-2]
        # following lines generate a 'Loop' column in the DataFrame using the metadata
        counter = 1
        df = []
        for entry in self.loops[0:-1]:
            singleloop = dataframe.iloc[entry:entry+self.Nrpointsloop].copy()
            singleloop.loc[:,'Loop'] = counter
            counter += 1
            df.extend(singleloop.values.tolist())
        # generating dataframe with selected columns
        self.CVdf = pd.DataFrame(df, columns = ['time', 'control_V', 'Ewe', '<I>', 'cycle number', '(Q-Qo)', 'I Range', '<Ece>', 'mode', 'ox/red', 'error', 'control changes', 'counter inc.', 'uts', 'Loop'])
        # reducing amount of columns in dataframe to relevant 
        self.CVcutdf = self.CVdf.iloc[:,[2, 3, -1]]
        # adding column with Re(Z) values of every Loop
        correction = ReZlist.set_index('Loop').to_dict()['Re(Z)']
        self.CVcorrdf = self.CVdf
        for value in ReZlist['Loop']:            
            self.CVdf.loc[self.CVdf['Loop'] == value, 'Re(Z)'] = correction[int(value)]
        # correcting the voltage of  Ewe by Reference elektrode, OER potential and resistivity measured with PEIS and putting the values in new column called Eta
        self.CVcorrdf['Eta'] = self.CVdf['Ewe'] - RefElek -1.23 - (self.CVdf['Re(Z)'] * self.CVdf['<I>']/1000)
        #reducing the number of columns to relevant  
        self.CVcorrcutdf = self.CVcorrdf.iloc[:,[-3, 2, -2, -1, 3]]
        # following lines generate dataframe with every loop as columnpair to be saved 
        safefile = pd.DataFrame()
        columnlist = []
        new_filename = str(filename) + "_data.csv"
        new_file_path = os.path.join(path,'CVevaluation', new_filename)
        for loop in self.CVcorrcutdf['Loop'].unique():
            columnlist.append('Loop_' +str(int(loop))+'_Eta/V')
            columnlist.append('Loop_' +str(int(loop))+'_I/mA')
            cut = self.CVcorrcutdf[self.CVcorrcutdf['Loop'] == loop].iloc[:, 3:5].reset_index(drop=True)
            safefile = pd.concat([safefile, cut], names = columnlist, axis = 1, join = 'outer', ignore_index= True).round(5)
        safefile.columns = columnlist
        safefile.to_csv(new_file_path, index = None, sep= ',')
        # makign graphic of data that will be saved as file
        plotofloopsCV(self.CVcorrcutdf, '_CVplots', 'CVevaluation', 'Eta', '<I>', r'$Eta$ /mV', 'I/mA')
        return

#class is used to process all imported datafiles resulting from Modular Potentio-technique    
class ModularPotentio():
    def __init__(self):
        # try making subfolder to save later results in
        try:
            os.makedirs('MPevaluation')
        except:
            print("Path 'MPevaluation' already exists")
        # following lines construct loop-column from given cycle numbers 'Ns' in Dataframe 'data'
        cyclenumber = 0
        loop = 1
        looplist = []
        samplenumber = 1
        for index, value in data['Ns'].items(): 
             if value == cyclenumber:
                 looplist.append([index, loop])
             if value > cyclenumber:
                 cyclenumber += 1
                 looplist.append([index, loop])
             elif value < cyclenumber:
                 loop += 1
                 cyclenumber = 0
                 looplist.append([index, loop])
        # adding determined Loop.column to dataframe
        self.loops = pd.DataFrame(looplist, columns = ['index', 'Loop'])             
        data['Loop']=self.loops['Loop']
        self.MPdf = data.copy(deep=True)
        #reducing dataframe to relevant columns
        self.MPcutdf = self.MPdf.iloc[:,[-1, 1, 2]]
        # getting Re(Z) list determined from PEIS
        correction = ReZlist.set_index('Loop').to_dict()['Re(Z)']
        self.MPcorrdf = self.MPdf.copy(deep=True)
        #adding Re(Z) values as column to every row dependent on loop
        for value in ReZlist['Loop']:            
            self.MPcorrdf.loc[self.MPcorrdf['Loop'] == value,'Re(Z)'] = correction[int(value)]
        # correcting the voltage of  Ewe by Reference elektrode, OER potential and resistivity measured with PEIS and putting the values in new column called Eta
        self.MPcorrdf['Eta'] = self.MPcorrdf['Ewe/V'] - RefElek -1.23 - (self.MPcorrdf['Re(Z)'] * self.MPcorrdf['I/mA']/1000)
        # reducing number of columns to relevant
        self.MPcorrdf = self.MPcorrdf.iloc[:,[0, 2, 5, -2, -1,  3, -3]]
        # correcting current values by area
        self.MPcorrdf['j/mA/cm^2'] = self.MPcorrdf['I/mA'].div(area)
        # reducing number of columns to relevant
        self.MPcorrcutdf = self.MPcorrdf.loc[:,'Eta':'j/mA/cm^2'].drop(['I/mA'], axis = 1)
        # calculating mean values of current and potential at end of every cycle of each loop
        MPmeans = []
        for loop in self.MPcorrdf['Loop'].unique():
            self.MPcorrdfloop = self.MPcorrdf.loc[self.MPcorrdf['Loop'] == loop]
            for cycle in self.MPcorrdf['Ns'].unique():
                MPmean = self.MPcorrdfloop.loc[self.MPcorrdfloop['Ns'] == cycle]
                MPmeanEiR = MPmean['Eta'].iloc[-21:-1].mean()
                MPmeanI = MPmean['j/mA/cm^2'].iloc[-21:-1].mean()
                MPmeans.append([loop, cycle, MPmeanEiR, MPmeanI])
        # generating new dataframe with mean values of each loop and cycle
        self.MPmeans = pd.DataFrame(MPmeans, columns = ['Loop', 'Cycle', 'Eta', 'j/mA/cm^2'])
        tafellist = []
        # calculating tafeldiagramm for each loop and determining tafelslope as smallest slope determined on 4 points within datarange  
        for loop in self.MPmeans['Loop'].unique():
            self.tafelplot = self.MPmeans.loc[self.MPmeans['Loop'] == loop].copy(deep=True)
            #calculating log10-values of each current value unless it is 0
            self.tafelplot.loc[:,'log(j/mA/cm^2)'] = [np.log10(x) if x > 0 else None for x in self.tafelplot['j/mA/cm^2']] 
            tafelfit = self.tafelplot.dropna(subset = ['log(j/mA/cm^2)']).iloc[1:-1,:].reset_index(drop=True)
            fitlist = []
            # looking for lowest slope determined by 4 points within datarange
            for value in range(len(tafelfit)-3):
                higherlimit = value + 4
                # determining dataframe rows in dataframe for linear fit
                tafelfitcut = tafelfit.iloc[value:higherlimit,:]
                # calculating linear fit for each current-potential-value-pair
                linear_model=np.polyfit(tafelfitcut['log(j/mA/cm^2)'], tafelfitcut['Eta'],1)
                fitlist.append(linear_model)
            #generating list of fitting results and sorting the values by increasing absolute slope value
            fitlist = pd.DataFrame(fitlist, columns = ['slope', 'y_at_x=0']).abs().sort_values('slope').reset_index(drop=True)
            # finding and saving only lowest determined slope within datarange
            tafelcurve = fitlist.iloc[0,:]
            LoopCurrentMax = self.MPmeans.loc[self.MPmeans['Loop'] == loop].loc[:,'j/mA/cm^2'].max()
            #determining potential values at defined current values
            if LoopCurrentMax <= 10: 
                Eat5 = getpointsnexttoj(self.tafelplot, 5)
                Eat10 = tafelcurve.iloc[0] + tafelcurve.iloc[1]
                Eat50 = None
                Eat100 = None
                Eat10valid = 'extrapolated'
            elif LoopCurrentMax <= 50:
                Eat5 = getpointsnexttoj(self.tafelplot, 5)
                Eat10 = getpointsnexttoj(self.tafelplot, 10)
                Eat50 = None
                Eat100 = None
                Eat10valid = 'valid'
            elif LoopCurrentMax <= 100:
                Eat5 = getpointsnexttoj(self.tafelplot, 5)
                Eat10 = getpointsnexttoj(self.tafelplot, 10)
                Eat50 = getpointsnexttoj(self.tafelplot, 50)
                Eat100 = None
                Eat10valid = 'valid'
            else:
                Eat5 = getpointsnexttoj(self.tafelplot, 5)
                Eat10 = getpointsnexttoj(self.tafelplot, 10)
                Eat50 = getpointsnexttoj(self.tafelplot, 50)
                Eat100 = getpointsnexttoj(self.tafelplot, 100)
                Eat10valid = 'valid'
            #putting all determined values from tafelcurve per loop into a list
            tafellist.append([loop, tafelcurve.iloc[0], tafelcurve.iloc[1], Eat5, Eat10, Eat50, Eat100, Eat10valid])
            # plot and save tafelplot + fit of every loop
            plotandfit(self.tafelplot['log(j/mA/cm^2)'], self.tafelplot['Eta'], tafelcurve.iloc[0], tafelcurve.iloc[1], 'Loop '+ str(loop)+ ' Tafelslope: ' + str(tafelcurve.iloc[0].round(3)) + ' mV/dec', 'MPevaluation','Tafelplot_loop_' + str(loop)+ '.png')            
        # making tafelplot data to DataFrame
        self.tafellist = pd.DataFrame(tafellist, columns = ['Loop', 'Tafelslope', 'y_at_x=0', r'$eta$_at_5mA/cm2', r'$eta$_at_10mA/cm2', r'$eta$_at_50mA/cm2', r'$eta$_at_100mA/cm2', r'$eta$_at_10 mA/cm2 in datarange'])   
        # ploting and saving every Eat value to a loop dependend graphic
        Eatoverloop(self.tafellist, 'MPevaluation', '_EtaoverLoop')
        #plotting and saving the determined tafelslope dependend on the loop number
        plotsave(self.tafellist['Loop'], self.tafellist['Tafelslope'], 'MPevaluation', '_Tafelslopevsloop', 'Loop', 'Tafel slope / mV/dec')
        # following lines create and save a dataframe with current-potential-value-pair for each loop as columns to be saved in subfolder
        samplenumber += 1
        safefile = pd.DataFrame()
        columnlist = []
        for loop in self.MPcorrcutdf['Loop'].unique():
            columnlist.append('Loop_' +str(int(loop))+'_Eta/V')
            columnlist.append('Loop_' +str(int(loop))+'_I/mA')
            cut = self.MPcorrcutdf[self.MPcorrcutdf['Loop'] == loop].drop(['Loop'], axis = 1).reset_index(drop=True)
            safefile = pd.concat([safefile, cut], names = columnlist, axis = 1, join = 'outer', ignore_index= True).round(5)
        safefile.columns = columnlist       
        new_filename = str(filename) + "_data.csv"
        new_file_path = os.path.join(path,'MPevaluation', new_filename)
        safefile.to_csv(new_file_path, index = None, sep= ',')
        #following lines create and save a dataframe with current-potential-mean-values for each loop and cycle as columns to be saved in subfolder
        safefile = pd.DataFrame()
        columnlist = []
        for loop in self.MPmeans['Loop'].unique():
            columnlist.append('Loop_' +str(int(loop))+'_Eta/V')
            columnlist.append('Loop_' +str(int(loop))+'_j/mA/cm^2')
            cut = self.MPmeans[self.MPmeans['Loop'] == loop].drop(['Cycle', 'Loop'], axis = 1).reset_index(drop=True).round(5)
            safefile = pd.concat([safefile, cut], names = columnlist, axis = 1, join = 'outer', ignore_index= True)
        safefile.columns = columnlist
        new_filename = str(filename) + "_MPmean.csv"
        new_file_path = os.path.join(path,'MPevaluation', new_filename)
        safefile.to_csv(new_file_path, index = None, sep= ',')
        # generating log-value-column of current values in list of mean-value-Dataframe
        self.MPmeans['log(j/mA/cm^2)'] = [np.log10(x) if x > 0 else None for x in self.MPmeans['j/mA/cm^2']]
        #plot and save mean values of tafel curve of every loop in one graphic
        plotofloops(self.MPmeans, "_MPtafelplots.png", 'MPevaluation', 'log(j/mA/cm^2)', 'Eta', 'log(j)', r'$Eta$ /mV')
        #saving dataframe with tafelslopes to text-file
        new_filename = str(filename) + "_MPtafelresults.csv"
        new_file_path = os.path.join(path,'MPevaluation', new_filename)  
        self.tafellist.round(5).to_csv(new_file_path, index = None, sep= ',')
        return

#class is used to process all imported datafiles resulting from CP-technique
class CP():
    def __init__(self):
        try:
            os.makedirs('CPevaluation')
        except: print("path 'CPevaluation' already exists")
        # these first lines find out at which lines in the dataframe correspond to each loop, by looking into the metadata (meta), because the data of every loop is appended to the datafile
        self.loops = meta['loops']['indexes']
        self.Nrpointsloop = len(dataframe)-self.loops[-2]
        # following lines generate a 'Loop' column in the DataFrame using the metadata
        counter = 1
        df = []
        for entry in self.loops[0:-1]:
            singleloop = dataframe.iloc[entry:entry+self.Nrpointsloop].copy()
            singleloop.loc[:,'Loop'] = counter
            counter += 1
            df.extend(singleloop.values.tolist())
        self.CPdf = pd.DataFrame(df, columns = ['Ns', 'time', 'control_I', 'Ewe', 'I', 'dQ', '(Q-Qo)', 'half cycle', 'Q charge/discharge', 'I Range', 'Ece', 'mode', 'ox/red', 'error', 'control changes', 'Ns changes', 'counter inc.', 'uts', 'Loop'])
        #cutting dataframe to relevant number of columns
        self.CPcutdf = self.CPdf.iloc[:,[1, 3, 4, -1]]
        # adding Re(Z)-values of each loop to column
        correction = ReZlist.set_index('Loop').to_dict()['Re(Z)']
        self.CPcorrdf = self.CPdf
        for value in ReZlist['Loop']:            
            self.CPcorrdf.loc[self.CPcorrdf['Loop'] == value, 'Re(Z)'] = correction[int(value)]
        # correcting the voltage of  Ewe by Reference elektrode, OER potential and resistivity measured with PEIS and putting the values in new column called Eta
        self.CPcorrdf['Eta'] = self.CPdf['Ewe'] - RefElek -1.23 -(self.CPdf['Re(Z)'] * Is/1000)
        # generating new time-column that starts at 0 for every new loop
        time = []
        for loop in self.CPcorrdf['Loop'].unique():
            time.extend(self.CPcorrdf.loc[self.CPcorrdf['Loop'] == loop]['time'] - self.CPcorrdf.loc[self.CPcorrdf['Loop'] == loop]['time'].min())
        self.CPcorrdf['time'] = pd.DataFrame(time).astype(int)
        #reducing number of columns to relevant
        self.CPcorrcutdf = self.CPcorrdf.iloc[:,[1, -3, 3, -2, -1, 4]]
        #plotting and saving potential vs time plots for each loop in one graphic
        plotofloops(self.CPcorrcutdf, '_CPsvstime','CPevaluation', 'time', 'Eta', 't /s', r'$Eta$ /mV')
        #determining the mean potential for each loop using the last 50 values of each loop
        EiRmean = []
        for loop in self.CPcorrdf['Loop'].unique():
            df = self.CPcorrdf.loc[self.CPcorrdf['Loop'] == loop]
            Emean = df['Eta'].iloc[-51:-1].mean()
            EiRmean.append([loop, Emean])
        self.EiRmean = pd.DataFrame(EiRmean, columns = ['Loop', 'Eta_mean of last 50 values'])
        # generating a dataframe with time and current-potential-value-paris for each loop 
        safefile = self.CPcorrcutdf.loc[self.CPcorrcutdf['Loop'] == 1].iloc[:,0]
        columnlist = ['Loop_' +str(int(loop))+'_t/s']
        for loop in self.CPcorrcutdf['Loop'].unique():
            columnlist.append('Loop_' +str(int(loop))+'_Eta/V')
            columnlist.append('Loop_' +str(int(loop))+'_I/mA')           
            cut2 = self.CPcorrcutdf[self.CPcorrcutdf['Loop'] == loop].iloc[:, 3:5].reset_index(drop=True).round(5)
            safefile = pd.concat([safefile, cut2], names = columnlist, axis = 1, join = 'outer', ignore_index= True)
        # saving this generated dataframe to file in subfolder
        safefile.columns = columnlist
        new_filename = str(filename) + "_data.csv"
        new_file_path = os.path.join(path,'CPevaluation', new_filename)
        safefile.to_csv(new_file_path, index = None, sep= ',')
        #saving mean value datafrmae to text file in subfolder
        new_filename = str(filename) + "_Eta_mean.csv"
        new_file_path = os.path.join(path,'CPevaluation', new_filename)
        self.EiRmean.round(5).to_csv(new_file_path, index = None, sep= ',')
        # plot and save mean-EiR-values vs Loops as graphic in subfolder
        plotsave(self.EiRmean['Loop'], self.EiRmean['Eta_mean of last 50 values']*1000, 'CPevaluation', '_Etameanvsloop', 'Loop', r'$Eta$ /mV')
        return



#entry of filepath
os.chdir(input('Enter full path of data folder: '))
#Data folder with .csv-data of one sample to do math
path = os.getcwd()
#creates a list of file names
all_files = glob.glob(path + '/*.mpr') 
all_files.sort()
#get samplename from current sample
samplename = str(all_files[0].rsplit('.', 1)[0].split('\\')[-1].rsplit('_', 2)[0])
#ask for potential of current reference electrode
RefElek = None # getvalue('Potential of reference electrode in V: ', float) 
#ask for area of current working electrode
area = getvalue('Geometric Electrode Area in cm^2 of "'+ str(samplename) +'": ', float)
# for each mpr file find technique and do specific class related stuff
for entry in all_files:
    filename  = str(entry.rsplit('.', 1)[0].split('\\')[-1])
    technique = str(filename.rsplit('_', 2)[1])
    method = str(filename.rsplit('_', 1)[1])
    print('Working on technique ' +technique+ ' containing ' +method+ ' method')
    # MP technique is included in galvani package an returns no meta data, therefore it has to be treated differently to other techniques
    if 'MP.mpr' in entry: 
        technique = 'MP'
        mpr_file = BioLogic.MPRfile(entry)
        data = pd.DataFrame(mpr_file.data)
        data = data.drop(['flags', 'time/s', 'I Range'], axis = 1)
        data = ModularPotentio()
    # following techniques are included in eclabfiles package and meta data (meta) and measurement data (data) can be used 
    else: 
        data, meta = ecf.process(str(entry.split('\\')[-1]))
        dataframe = pd.DataFrame(data)
        if 'Rcmp' in dataframe.columns:
            dataframe = dataframe.drop('Rcmp', axis = 1)
        technique = meta['settings']['technique']        
        if technique == 'PEIS':
            data = PEIS()
            ReZlist = data.ReZlist            
            data.ReZlist.to_csv(str(filename) + '_Re(Z)_list.txt', index = False)
        if technique == 'OCV':
            data = OCV()
        if technique == 'ZIR':
            data = ZIR()
        if technique == 'LSV':
            data = LSV()
        if technique == 'CV':
            if RefElek == None:
                RefElek = meta['params'][0]['Ei'] #technique 4
            else:
                RefElek = meta['params'][0]['Ei'] #technique 8
            data = CV()
        if technique == 'CP':
            if 'Energy charge' in dataframe.columns:
                dataframe = dataframe.drop(['Energy charge', 'Energy discharge', 'Capacitance charge', 'Capacitance discharge'], axis = 1)                   
            Is = meta['params'][0]['Is']
            print(Is)
            data = CP()
            

print('Done!')
