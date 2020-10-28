import re
import csv 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle 


  


def readFile(file, min_x, min_y, max_x, max_y, max, num):

	with open(file) as f:
	    data = [ele.split(',') for ele in f]
	
	limited = []
	max = max*num
	## Python will convert \n to os.linesep
	for i in range(len(data)):
		if len(data[i]) == 3: 
			for j in range(num):
				limited.append([int(data[i][0]), float(data[i][1]), float(data[i][2])+(60000)*j])

	print(len(limited))
	# max_x = (max_x)*float(num+1) 
	# max_x = (max_x)+48021.6
	max_y = (60000)*float(num)+max_y 
	# writing to csv file  
	with open('data/newData/Point_Of_Interest_modified_xx10.csv', 'w', newline='') as csvfile:  
	    # creating a csv writer object  
	    csvwriter = csv.writer(csvfile)  
	        
	    # writing the data rows  
	    csvwriter.writerow([int(len(limited))])	
	    csvwriter.writerows(limited)
	    csvwriter.writerow([])
	    csvwriter.writerow([float(min_x), float(min_y), float(max_x), float(max_y)])	

readFile('data/Point_Of_Interest_modified.csv', 0.0, 0.0, 48021.6, 55014.9, 20563, 10)		
