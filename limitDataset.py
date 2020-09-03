import re
import csv 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle 


  


def readFile(file, min_x, min_y, max_x, max_y):

	with open(file) as f:
	    data = [ele.split(',') for ele in f]
	
	limited = []
	max = 0
	## Python will convert \n to os.linesep
	for i in range(len(data)):
		if len(data[i]) == 3 and max_x > float(data[i][1]) and max_y > float(data[i][2]) and min_x <= float(data[i][1]) and min_y <= float(data[i][2]):
			limited.append([int(data[i][0]), float(data[i][1]), float(data[i][2])])
			max += 1

	# writing to csv file  
	with open('data/newData/Seattle2012_bt_2_0.csv', 'w') as csvfile:  
	    # creating a csv writer object  
	    csvwriter = csv.writer(csvfile)  
	        
	    # writing the data rows  
	    csvwriter.writerow([int(max)])	
	    csvwriter.writerows(limited)
	    csvwriter.writerow([])
	    csvwriter.writerow([float(min_x), float(min_y), float(max_x), float(max_y)])	

readFile('data/newData/Seattle2012_bt_2.csv', 10674.2, 0, 12000.0, 18243.4)		
