#subset_sum.py
import sys
def subset(array, target):
	sol = [[False for x in range(target + 1)] for x in range(len(array) + 1)]
	for i in range(len(array)+1):
		sol[i][0] = True
	for i in range(1,(len(array)+1)):
		for j in range(1, target+1):
			if (j - array[i-1] >= 0):
				sol[i][j] = sol[i-1][j] | sol[i-1][j - array[i-1]]
			else:
				sol[i][j] = sol[i-1][j]
    
	if sol[len(array)][target]:
		return True
	else:
		return False
	'''
	if num == 0 or num < 1:
		return False
	elif len(array) == 0:
		return False
	else:
		if array[0] == num:
			return True
		else:
			return subset(array[1:],(num - array[0])) or subset(array[1:],num) 
	'''
def printSub(sol, array, target):

	if(sol[len(array)][target]):
		print("Found!")
		i = len(array)
		j = target
		while(j!=0):
			if(sol[i-1][j] == True):
				i-=1
			else:
				sys.stdout.write(str(array[i-1]))
				sys.stdout.write(' ')
				j = j - array[i-1]
	else:
		print("No combination found! ")

sol = subset([1,2,7], 6)
