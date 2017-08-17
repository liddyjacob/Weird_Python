
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


import sys
import sympy
import bisect
from math import log, ceil
from sympy import sieve
from sympy import ntheory
sieve._reset()

prime = sympy.ntheory.generate.nextprime(5,1)
#TOOLS ---
def b(pset, expset):
	product = 1
	pexpset = zip(pset, expset)
	for (prime,exp) in pexpset:
		product *= (prime ** (exp + 1) - 1) / (1.0*(prime - 1) * prime ** exp)

	return product

def b_inf(pset):
	product = 1
	for prime in pset:
		product *= (prime + 0.0) / (prime - 1.0)

	return product

def b_1(pset):
	product = 1
	for prime in pset:
		product *= (prime + 1.0) / (prime + 0.0)

	return product

def get_lower_bound(pset):
	b = b_1(pset)
	lower_bound = b / (2.0 - b)

	return lower_bound

def del_(prime, exp):
	return ((prime ** exp) - 1.0 ) * prime / (prime ** (exp + 1) - 1.0)

def del_pos(prime, exp):
	return (prime ** (exp + 2) - 1.0) / (prime * (prime ** (exp + 1) - 1.0))

def get_max_del(pset, expset):
	pexpset = zip(pset, expset)
	max_d = 0
	for (prime, exp) in pexpset:
		if del_(prime, exp) > max_d:
			max_d = del_(prime, exp)

	return max_d

def get_max_del_pos(pset, expset):
	pexpset = zip(pset, expset)
	max_dp = 0
	max_dp_index = 0
	index = 0
	for (prime, exp) in pexpset:
		if del_pos(prime,exp) > max_dp:
			max_dp = del_pos(prime, exp)
			max_dp_index = index
		index +=1

	return (max_dp, max_dp_index)


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

			
			


def add(primeset, expset):
	print(zip(primeset, expset))


"""
nextprime 
---------
Retrieve the next prime. Pass it, along with pset,  
to the b_n checker


"""
def nextprime(pset, low_bound, end_divisors):
	newprime = sympy.ntheory.generate.nextprime(low_bound, 1)

	return b_inf_check(pset, newprime, end_divisors)

"""

b_inf_n_chek
Check to see where the next prime should branch 
based of abundance.

if b_inf(n) * b_inf(p) < 2,
	send to check number of divisors without abundance[deficient(def)]

if b_inf(n) * b_inf(p) >= 2,
	send to check b_1_n


"""
def b_inf_check(pset, newprime, end_divisors):

	if ((b_inf(pset + [newprime]))	> 2.0):
		return b_1_check(pset, newprime, end_divisors)
	else:
		return number_divisors_def(pset, newprime, end_divisors)



"""
number_divisors_def
Number of divisors for this deficient primeset.

If |pset| + 1 = d, return false

Otherwise, branch off in two ways:
On one branch, add the prime to pset.
On another set, don't add the prime.


"""

def number_divisors_def(pset, newprime, end_divisors):
	
	if len(pset) + 1 == end_divisors:
		return False

	else:
		return add_dontadd(pset, newprime, end_divisors)



"""
add/dontadd

split up divisors:
	Try adding. If that works, try not adding.
"""
def add_dontadd(pset, newprime, end_divisors):

	print 'in add_dontadd with pset:{0} and newprime {1}'.format(pset, newprime)	
	if nextprime(pset + [newprime], newprime, end_divisors):
		print ("PSET WORKED")
		print(pset)

		return nextprime(pset, newprime, end_divisors)

	return False
"""
b_1_n
Check if b_1_n > 2

IF it is, then this better be the last prime

Otherwise, check the next prime back at the beginning
"""
def b_1_check(pset, newprime, end_divisors):
	
	if b_1(pset + [newprime]) > 2:
		return number_divisors_b_1(pset, newprime, end_divisors)
	else:
		return number_divisors_b_x(pset, newprime, end_divisors)


"""
number_divisors_b_1
Branch off based off number of divisors
if b_1 > 2 
"""

def number_divisors_b_1(pset, newprime, end_divisors):

	if len(pset) + 1 == end_divisors:
		explist = []
		for i in range(0, end_divisors):
			explist+=[1]
		add(pset + [newprime], explist)

		return True

	else:
		return find_good_prime(pset, newprime, end_divisors)

"""
number_divisors_b_x
branch off based off number of divisors
if b_1 < 2 and b_inf > 2


"""


def number_divisors_b_x(pset, newprime, end_divisors):

	if len(pset) + 1 == end_divisors:
		if find_exp_combos(pset + [newprime]):
			return nextprime(pset, newprime, end_divisors)
		else:
			return False
	else:
		return add_dontadd(pset, newprime, end_divisors)

"""
find_good_prime
Find a good prime
if current prime makes b_1 > 2 and the divisors is less than desired

"""

def find_good_prime(pset, newprime, end_divisors):
	lower_bound = get_lower_bound(pset)
	return nextprime(pset, lower_bound, end_divisors)


"""
find_exp_combos
Find all exponent combinations that make this primeset primative-abundant
If none are found, return false
If one is found, return true
"""

def find_exp_combos(pset):
	#print("finding exp combos")
	#print(pset)

	initial_exps = []
	for i in range(0, len(pset)):
		initial_exps+=[1]

	return find_exp_recursive(pset, initial_exps)

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


"""
exp_abundant
Determine if product is abundant
Check if primative if so
Otherwise increase exponents.

"""

def exp_abundant(pset, curr_exps):
	#print("Abundance check")
	#print(zip(pset, curr_exps))
	if b(pset, curr_exps) >= 2:
		return primative(pset, curr_exps)
		
		

	else:
		return exp_increase(pset, curr_exps)

"""
exp_increase
Increase exponents to make a primative abundant


"""
def exp_increase(pset, curr_exps):
	return del_pos_check(pset, curr_exps)


"""
del_pos_check
see if any del_pos makes 

"""
def del_pos_check(pset, curr_exps):
	#print("delta check")
	(max_dp, index) = get_max_del_pos(pset, curr_exps)
	if b(pset, curr_exps) * max_dp > 2:
		return find_min_exp_inc(pset, curr_exps)
	else:
		#Originally, raised the exponent that increased b_n the most:
		#curr_exps[index]+=1
		#del_pos_check(pset, curr_exps)
		return raise_exponents(pset, curr_exps)
"""
potential_abundant
See if raising many exponents beyond indexed prime
will make a number abundant eventually:

"""
def potential_abundant(pset, curr_exps, index):
	return b(pset[:index], curr_exps[:index]) * b_inf(pset[index:]) > 2
	


"""
raise exponents:
Raise all exponents that could eventually lead in an abundant number:
"""
def raise_exponents(pset, curr_exps):

	worked = False

	for index in range(0, len(pset)):
		if potential_abundant(pset, curr_exps, index):
			(max_dp, rel_index) = get_max_del_pos(pset[index:], curr_exps[index:])
			#True index of prime with max dp
			trueindex = rel_index + index
			new_exps = list(curr_exps)
			new_exps[trueindex]+=1
			#print 'trying: {0}'.format(zip(pset, new_exps, curr_exps))
			if del_pos_check(pset, new_exps):
				worked = True

	return worked


"""
find_min_exp_inc
find minimum exponent increase such that n is primative
"""

def find_min_exp_inc(pset, curr_exps):
	index_list = find_all_exp_inc(pset, curr_exps)

	if len(index_list) == 0:
		return False
	else:
		return True
	"""
	for index in index_list:
		new_exps = list(curr_exps)
		new_exps[index]+=1
		add(pset, curr
	"""
	#return primative(pset, curr_exps)
	

""" 
"""

"""
find_exp_recursive
Recursivly find all exponent combos for a primeset
Activated by find_exp_combos.

"""
def find_exp_recursive(pset, curr_exps):
	exp_abundant(pset, curr_exps)	


nextprime([],2,3)


