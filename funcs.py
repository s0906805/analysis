import numpy as np
from scipy import interpolate
import numpy.ma as ma



# function to return indices of missing values 
def getIndices(array,value):
    output = [] 
    index = 0 
    for i in array:
        if i <= value: 
            output.append(index)


        index += 1 
    return np.asarray(output) 
#-----------------------------------------------------------------------------#

# a function which 
#1.runs through the indicies of the missing data
#2.test if it's only a gap of one data point
#3.interpolates across the points either side
#4.saves points into a new array
def imput1(series):
    im1 = np.copy(series) # copy series
    inds = getIndices(series, -999)
    for elem in inds:
        pt = elem
        prev = elem -1
        sub = elem +1
        if series[prev]!= -999 and series[sub] != -999:
            x = np.linspace(0,2, num = 2)
            y = (series[prev], series[sub])
            f = interpolate.interp1d(x,y)
            xnew = np.linspace(0,2, num = 3)
            ynew = f(xnew)
            im1[pt] = ynew[1] # take middle value (the interpolated one)
           # print 'the new value of the missing point is:',im1[pt]
    return im1
#-----------------------------------------------------------------------------#

# a function which
# 1. reads in time step data, rad data and the start year
# 2. for each time slot find the matching time slots in the other years
# 3. average the available data at that time over the 5 years
# 4. fill in if average of 4 years 

def climo(time, series):
    styr = time[0,0] # get the first year
    series_orig = ma.masked_equal(series, -999)# create a masked array so can keep track of orginial entries
    idx = np.asarray(np.where(time[:,0]==styr)[0])
    av =[0]*len(np.asarray(idx)) # want an array to hold the averages
    ms =[0]*len(np.asarray(idx)) # and another to hold the info of how many missing numbers
    yr1 = time[np.asarray(idx),:]# all time data from the first year

    for elem in idx:
        midx =np.asarray(np.where(time[:,1] == yr1[elem,1])[0])#all the indexes which are in month = 1
        didx = np.asarray(np.where(time[midx,2] ==yr1[elem,2])[0])# all the indexes of these which are month = 1, day = 1
        hidx = np.asarray(np.where(time[didx,3] == yr1[elem,3])[0])# month = 1 day = 1 hour = 0 
        x = np.asarray(np.where(time[hidx,4] ==yr1[elem,4])[0]) # exactly the same time stamp
        d=  midx[didx[hidx[x]]]# index of full series where the date/time is the same
        ms[elem] = ma.count(series_orig[d]) # number of missing values
        av[elem] = np.mean(series_orig[d]) #average of non missing values

    return av, ms, yr1
#-----------------------------------------------------------------------------#
#a function which
#1.takes in the series and loops through missing value indices
#2.for each missing value, identifies the time stamp 
#3.pulls the climotology for that value and impute 
def imput2(time, series, avg, ms): #imput the time array, the series, and the outputs avg and ms from climo()
    im2 = np.copy(series) # make a copy
    missingInds = getIndices(series, -999)
    missing = [0]*len(im2)
    styr = time[0,0] # get the first year
    idx = np.asarray(np.where(time[:,0]==styr)[0])
    yr1 = time[np.asarray(idx),:]# all time data from the first year
    for elem in missingInds: # loop through the missing values of the data series
          timeStamp = time[elem, :] # time stamp of missing value
          # now look through the yr1 matrix to find the accompanying climatology
          midx =np.asarray(np.where(yr1[:,1] == timeStamp[1])[0])#month
          didx = np.asarray(np.where(yr1[midx,2] ==timeStamp[2])[0])# day
          hidx = np.asarray(np.where(yr1[didx,3] == timeStamp[3])[0])#hour
          x = np.asarray(np.where(yr1[hidx,4] ==timeStamp[4])[0]) #minute
          d=  midx[didx[hidx[x]]]# index of full series where the date/time is the same
          im2[elem] = avg[d]
          missing[elem] = ms[d]
    return im2, missing

#------------------------------------------------------------------------------#

def Count(array, number):
    return len([i for i in array if i == number])
