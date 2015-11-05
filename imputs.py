import numpy as np
from scipy import interpolate
import numpy.ma as ma
import matplotlib.pyplot as plt
import csv
from funcs import *

#read in data
ts = np.loadtxt('time_stamp.csv', delimiter = ',')
ser = np.loadtxt('olw_niamey.csv', delimiter = ',')
newfilename = 'olw_niamey_IC.csv'


print 'done reading'

#number of missing data points
n_ser = len(getIndices(ser, -999))
print 'original number of missing values', n_ser          

#impute 1    
ser_im1 = imput1(ser)
n_ser_im1 = n_ser - len(getIndices(ser_im1, -999))
print 'number of first imputations:' , n_ser_im1

#impute 2
#get climatology
output_im1 =climo(ts, ser_im1)
avg_ser = np.asarray(output_im1[0])
ms_ser = np.asarray(output_im1[1])
print 'climotology done'

#imput climotology into series
output_im2 = imput2(ts, ser_im1, avg_ser, ms_ser)
ser_im2 = np.asarray(output_im2[0])
ser_climNum = np.asarray(output_im2[1]) # gives the number of years each climatology slot was created from
#n_ser_im2 = n_ser - n_ser_im1

#print stats
#print 'number of second imputations:', n_ser_im2
print 'check', Count(ser_climNum,0)
print 'of which use a 5 year climotology:', Count(ser_climNum,5)
print 'of which use a 4 year climotology:', Count(ser_climNum,4)
print 'of which use a 3 year climotology:', Count(ser_climNum,3)
print 'of which use a 2 year climotology:', Count(ser_climNum,2)

ser_ma = ma.masked_equal(ser, -999)
ser_im1_ma = ma.masked_equal(ser, -999)


plt.figure(1)
plt.subplot(411)
plt.plot(ser)
plt.subplot(412)
plt.plot(ser_im1)
plt.subplot(413)
plt.plot(ser_im2, 'b.-', ser_im1_ma, 'gx-' , ser_ma, 'r.-')
plt.subplot(414)
plt.plot(ms_ser)
plt.show()

#save
np.savetxt(newfilename, ser_im2, delimiter=",")
print 'done writing'
