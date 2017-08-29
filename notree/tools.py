#TOOLS.PY
#the tools needed to search for all primative abundant odd numbers
import sys
import sympy
import bisect
from math import log, ceil
from sympy import sieve
from sympy import ntheory
sieve._reset()


"""
b
the b(n) function, where n is a positive integer
Calculates b(n), where n is in factored form

PARAMETERS: 
pset - the set of primes that compose n
expset - the exponents of the respective primes

"""
def b(pset, expset):
    product = 1
    pexpset = zip(pset, expset)
    for (prime,exp) in pexpset:
        product *= (prime ** (exp + 1) - 1) / (1.0*(prime - 1) * prime ** exp)

    return product

"""
b_inf
the b(n) function, where n is a positive integer,
as the exponents of the prime factors of n approach infinity.

Calculates b_inf(n), where n is in factored form

PARAMETERS: 
pset - the set of primes that compose n

"""

def b_inf(pset):
    product = 1
    for prime in pset:
        product *= (prime + 0.0) / (prime - 1.0)

    return product

"""
b_1
Calculate b(p_1 * p_2 * p_3 ... )
b with all exponents as one

PARAMETERS:
pset - the set of primes that compose n
"""
def b_1(pset):
    product = 1
    for prime in pset:
        product *= (prime + 1.0) / (prime + 0.0)

    return product
                    

"""
get_lower_bound(pset)
In the case that a set of primes is nearly abundant with exponents 1
find the lowest number l such that if l is prime, b_1(pset + [l]) > 2

PARAMETERS:
pset
"""

def get_lower_bound(pset):
    b = b_1(pset)
    lower_bound = b / (2.0 - b)

    return lower_bound

"""
del_ and del_pos

Calculate the del/del_pos of a prime and exponent
"""


def del_(prime, exp):
    return ((prime ** exp) - 1.0 ) * prime / (prime ** (exp + 1) - 1.0)

def del_pos(prime, exp):
    return (prime ** (exp + 2) - 1.0) / (prime * (prime ** (exp + 1) - 1.0))

"""
get_max_del and get_max_del_pos

Out of a set of primes and exponents, find the pair with maximum del/del_pos

"""

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

def add(primeset, expset):
	with open("primatives.txt", "a") as addfile:
		addfile.write('{0}\n'.format(zip(primeset, expset)))
	print(zip(primeset, expset))
	ab = abundance(primeset, expset)
	div = divisors_below_a(primeset, expset)
	print '---ABUNDANCE:{0}'.format(ab)
#	print '---DIVISORS BELOW ^:{0}'.format(div)
#	if len(div) < 100:
#		sub = subset(div, ab)
#		print '---PSEUDOPERFECT:{0}'.format(sub)
#
#		if not sub:
#			with open("WEIRDS.txt", "a") as addfile:
#				addfile.write('{0}\n'.format(zip(primeset, expset)))

def a(primeset, expset):
	product = 1
	for (prime, exp) in zip(primeset, expset):
		product  *= (prime ** (exp + 1) - 1) / (1.0 * (prime - 1))

	
	return product

def prod(primeset, expset):
	product = 1
	for (prime, exp) in zip(primeset, expset):
		product *= (prime ** (exp))


	return product


def abundance(primeset, expset):
	return a(primeset, expset) - 2 * prod(primeset, expset)


def divisors_below_a(primeset, expset):
	#Divisors to return:
	divisors = [1]
	a = abundance(primeset, expset)
	#print abundance
	#print ("a")
	primelets = zip(primeset, expset)
	for i in range(0, len(primelets)):
		p_product = 1
		prev_d_size = len(divisors)
		for exp in (1, primelets[i][1]):
			p_product *= primelets[i][0]
			#Do not allow larger than abundance
			if p_product > abundance:
				break
	
			for div_index in range(0, prev_d_size):
				if (p_product * divisors[div_index] <= abundance):
					divisors.append(p_product * divisors[div_index])
				else:
					break
	return divisors

def subset(array, target):
    if target == 0 or target < 1:
        return False
    elif len(array) == 0:
        return False
    else:
        if array[0] == target:
            return True
        else:
            return subset(array[1:],(target - array[0])) or subset(array[1:],target) 


