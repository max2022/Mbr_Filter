import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle 


  


def readFile(file):

	data = pd.read_csv(file, header=None).values

	table = np.zeros((0, 11), float)
	for i in range(len(data)):
		row = np.empty(1)
		row[0] = int(data[i, 0])
		for j in range(1, len(data[i])):
			str1 = data[i, j].replace("(", "")
			str1 = str1.replace(")", "")
			row = np.append(row, float(str1))

		table = np.vstack((table, row))
	
	
	return table;


def plotFile(table, cmbr):
	color = ['red', 'green', 'blue', 'black', 'pink']
	fig1 = plt.figure() 
	ax = fig1.add_subplot(111) 

	for i in range(len(table)):
		w = abs(table[i, 3] - table[i, 5])
		h = abs(table[i, 4] - table[i, 2])

		rect = Rectangle((table[i, 1], table[i, 2]), w, h, color = 'green', fill=False)   
		ax.add_patch(rect) 
	
	for i in range(len(cmbr)):
		w = abs(cmbr[i, 1] - cmbr[i, 3])
		h = abs(cmbr[i, 4] - cmbr[i, 6])
		print(cmbr[i, 0])
		if cmbr[i, 0] == 2:			
			rect = Rectangle((cmbr[i, 7], cmbr[i, 8]), w, h, color = 'red')   
			ax.add_patch(rect)

	plt.xlim([80, 1200]) 
	plt.ylim([80, 1200]) 
	plt.axis('tight')

	plt.show() 

plotFile(readFile('MBRAll.csv'), readFile('CMBRAll.csv'))
