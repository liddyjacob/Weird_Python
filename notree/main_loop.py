#MAIN_LOOP.PY
#The main loop that finds sets of primes
#that could lead to primative abundant numbers:
import sys
import sympy
import bisect
from math import log, ceil
from sympy import sieve
from sympy import ntheory
sieve._reset()
sys.setrecursionlimit(3000)

from tools import *
from exponents import *
DEBUG = False
DEBUG_TABS = 0
"""
nextprime:

Find the next prime, and pass on the set of current primes 
plus the new prime to a b_inf checker

PARAMETERS:
pset - the set of primes:
low_bound - the lower bound on where to start searching for primes:
end_divisors - how many divisors are wanted at the end of the cycle?

"""

def nextprime(pset, low_bound, end_divisors):
	if DEBUG:
		print("In nextprime")

	newprime = sympy.ntheory.generate.nextprime(low_bound, 1)


	return b_inf_check(pset, newprime, end_divisors)


"""
b_inf_check

Check to see if primeset could potentially be abundant, with right exponents:

If it can, 
	check b_1(pset) (RIGHT BRANCH)

If it can't, 
	find out weather to add a new prime or if no more primes can be added
	(LOWER BRANCH)

PARAMETERS:
pset
newprime - new prime gathered from nextprime
end_divisors
"""
def b_inf_check(pset, newprime, end_divisors):
	if DEBUG:
		print 'b_inf(pset + [newprime]):'.format(b_inf(pset + [newprime]))


	if ((b_inf(pset + [newprime]))  > 2.0):
		return b_1_check(pset, newprime, end_divisors)
	else:
		return number_divisors_def(pset, newprime, end_divisors)


#-----------RIGHT BRANCH--------------#
"""
number_divisors_def
For deficient primesets (3 7 11), are we at the maximum allowed divisors?
If so this primeset can't be abundant, so scrap it

Otherwise, move on to the next prime (3) -> (3 5). If that works, try the prime
after that (3) -> (3 7)

"""


def number_divisors_def(pset, newprime, end_divisors):
	if DEBUG:
		print 'In right branch, with set {0}'.format(pset + [newprime])

	if len(pset) + 1 == end_divisors:
		return False

	else:
		return add_dontadd(pset, newprime, end_divisors)


#-----------LOWER BRANCH------------#
"""
b_1_check
If b_inf > 2, then we need to check if the primeset
is abundant with all exponents 1.

If this is greater than 2, then this is a primative abundant.
	If the number of divisors is correct, then add this number
	If not, find next prime that won't do this

If less than 2, determine weather to find exponents or add primes.
"""

def b_1_check(pset, newprime, end_divisors):
	if DEBUG:
		print 'In lower branch, with {0}'.format(pset + [newprime])


	if b_1(pset + [newprime]) > 2:
		return number_divisors_b_1(pset, newprime, end_divisors)
	else:
		return number_divisors_b_x(pset, newprime, end_divisors)


"""
number_divisors_b_1
If b_1 > 1, check number of divisors
	If d, then number is primative abundamt with desired divisors.
	Add
		
	If < d, then prim. abundant with less than desired divisors.
	Find next prime that wont do this.
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
b_inf < 2 < b_1
So if the number of divisors is as desired, we should have a few primative abundants - 
	Check exponents
If not, then add/dontadd prime

"""
def number_divisors_b_x(pset, newprime, end_divisors):
	if DEBUG:
		print ' in b_x'

	if len(pset) + 1 == end_divisors:
		if find_exp_combos(pset + [newprime]):
			if DEBUG:
				print 'Exponent combos found'
			nextprime(pset, newprime, end_divisors)
			return True
		else:
			return False
	else:
		return add_dontadd(pset, newprime, end_divisors)



"""
find_good_prime
If b_1 > 2 with current prime, and divisors is less than desired,
find a new prime such that b_1 < 2
"""

def find_good_prime(pset, newprime, end_divisors):
	lower_bound = get_lower_bound(pset)
	return nextprime(pset, lower_bound, end_divisors)
#----------BOTH BRANCHES-------------#


"""
add_dontadd
Add a prime, then if that works and creats an abundant, 
	see what would happen if you did not add a prime,
	and instead add the next prime after newprime.
"""

def add_dontadd(pset, newprime, end_divisors):
	
	if DEBUG:
		print 'in add_dontadd with pset:{0} and newprime {1}'.format(pset, newprime) 

	if nextprime(pset + [newprime], newprime, end_divisors):
		print ("PSET WORKED")
		print(pset)
		nextprime(pset, newprime, end_divisors)
		return True

	return False
"""
"""
TOTAL = 0
nextprime([], 2,4)
print 'total: {0}'.format(TOTAL)
