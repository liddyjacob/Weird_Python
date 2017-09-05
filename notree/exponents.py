#exponents.py
#The Exponent tree - find exponents to primesets:
from tools import *
DEBUG = False

#SPECIAL TOOLS:

"""
potential_abundant
See if raising many exponents beyond indexed prime
will make a number abundant eventually:

#TODO:
	This is the problem with the algorithm.
	It fails to try rasing other exponents, because it wants
	to guarentee a number will be primative by raising its maximum exponent

	What needs to happen is some algorithm that will determine which 

"""
def potential_abundant(pset, curr_exps, index):
    return b(pset[:index], curr_exps[:index]) * b_inf(pset[index:]) > 2





"""
primative - 
determine if abundant number is primative
"""
	
def primative(pset, curr_exps, return_condition = False):
	print("Primative check")
	max_dp = get_max_del(pset, curr_exps)

	if b(pset, curr_exps) * max_dp < 2:
		add(pset, curr_exps)
		print("^PASS")
		return True
		
	else:
		print (zip(pset,curr_exps))
		print("^FAILURE")
		return return_condition




#ALGORITHM
"""
find_exp_combos
Find all exponent combinations that make this primeset primative-abundant
If none are found, return false
If one or more is found, return true
"""

def find_exp_combos(pset):
	#print("finding exp combos")
	#print(pset)

	initial_exps = []
	for i in range(0, len(pset)):
		initial_exps+=[1]

	return find_exp_recursive(pset, initial_exps)


"""
find_exp_recursive
Recursivly find all exponent combos for a primeset
Activated by find_exp_combos.

"""
def find_exp_recursive(pset, curr_exps):
	number_found = [0]
	exp_abundant(pset, curr_exps, number_found)
	if DEBUG:
		print 'number of abundants found with primeset {0} in find_exp_recursive:{1}'.format(pset, number_found)
	return number_found != [0]


"""
exp_abundant

Determine if product is abundant
Check if primative if so
Otherwise increase exponents.
"""

def exp_abundant(pset, curr_exps, number_found):
	print 'In exp_abundance for {0}'.format(zip(pset, curr_exps))
	if b(pset, curr_exps) >= 2:
		if primative(pset, curr_exps):
			number_found[0]+=1
		

	else:
		del_pos_check(pset, curr_exps, number_found)

"""
del_pos_check
see if we can make a number primative abundant based off its del pos
values
"""
def del_pos_check(pset, curr_exps, number_found, startat = 0):
	if DEBUG:
		print 'delta check for {0}'.format(zip(pset, curr_exps))


	(max_dp, index) = get_max_del_pos(pset, curr_exps)
	if b(pset, curr_exps) * max_dp > 2:
		return find_min_exp_inc(pset, curr_exps, number_found)
	else:
		#Originally, raised the exponent that increased b_n the most:
		#curr_exps[index]+=1
		#del_pos_check(pset, curr_exps)
		return raise_exponents(pset, curr_exps, number_found, startat)

"""
find_min_exp_inc
find minimum exponent increase such that n is primative
"""

def find_min_exp_inc(pset, curr_exps, number_found):
	index_list = find_all_exp_inc(pset, curr_exps)

	if DEBUG:
		print 'found {0} primatives'.format(len(index_list))	

	number_found[0]+=len(index_list)

	if len(index_list) == 0:
		return False
	else:
		return True


def find_all_exp_inc(pset, expset):
	pexpset = zip(pset, expset)
	index_list = []
	index = 0
	min_b_diff = 1
	for(prime, exp) in pexpset:
		new_expset = list(expset)
		b_ = del_pos(prime,exp) * b(pset, expset)
		if b_ > 2:
			new_expset[index]+=1
			if primative(pset, new_expset):
				index_list+= [index]
		index+=1


	return index_list


"""
raise exponents:
Raise all exponents that could eventually lead in an abundant number:

This will find the sequence of indices in order from del_pos(p_i, m_i) 
greatest to least.

That will then be passed into a function


"""
def raise_exponents(pset, curr_exps, number_found, startat):
	if DEBUG:
		print 'Raising exponents for {0}'.format(zip(pset, curr_exps))

	index_seq = i_sequence_dp(pset, curr_exps)
	return exponent_inc_algorithm(pset, curr_exps, index_seq, number_found, startat)

"""
	i_sequence_dp
	Find the sequence of indices of primes 
		such that the corresponding prime's del_pos is greatest to least

	return that sequence.
"""

def i_sequence_dp(pset, curr_exps):
	sequence = []

	for index in range(0, len(pset)):
		sequence.append((index, del_pos(pset[index], curr_exps[index])))

	#Make/find the index sequence.
	sequence.sort(key=lambda x: x[1], reverse=True)
	index_seq = [x[0] for x in sequence]

	print index_seq
	return index_seq


"""
	exponent_increase_algoritm:
	See what exponents can and cant be increased.
	Based off list that is indexes of primes with greatest 
	current influence on b(n)


	Startat determines where on the list of indexes we should start
	(To minimize repeated primatives
"""

def exponent_inc_algorithm(pset, curr_exps, index_seq, number_found, startat = 0):

	worked = False

	print "startat = {}".format(startat)

	for i in range(startat, len(index_seq)):
		indices_b = index_seq[:i]
		
		if mixed_b(pset, curr_exps, indices_b) > 2:
			new_exps = list(curr_exps)
			exp_index = index_seq[i]
			new_exps[exp_index] += 1

			if del_pos_check(pset, new_exps, number_found, i):
				worked = True

		else:
			break

	return worked


"""
	mixed_b:
	Mixed b does the operation of mixing the right b(p^m)
	with the b_inf(p)'s

	The rule is:
		if the index of p is in indices_b
		 -then we should evaluate b(p^m)
		otherwise:
			-evaluate b_inf(p)


"""
def mixed_b(pset, curr_exps, indices_b):
	product = 1.0

	for i in range(0, len(pset)):
		prime = pset[i]
		exponent = curr_exps[i]
		if i in set(indices_b):
			product *= b([prime], [exponent])
		else:
			product *= b_inf([prime])

	return product





#i_sequence_dp([3,5,7], [3,1,2])
