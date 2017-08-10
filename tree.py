#from __future__ import print_function
import sys
import sympy
import bisect
from subset_sum import subset
from math import log, ceil
from sympy import sieve
from sympy import ntheory
sieve._reset()
#tree.py

#TODO: Make primesets_tried more efficient. Use a datar structar

def prod(primelets):
	product = 1
	for (prime, exp) in primelets:
		product *= prime ** exp

	return product


#A - the alpha function
def a(primelets):
	product = 1
	for (prime, exp) in primelets:
		product *= (prime ** (exp + 1) - 1) / (prime - 1)

	return product

def b(primelets):
	product = 1
	for (prime, exp) in primelets:
		product *= (prime ** (exp + 1) - 1) / (1.0 * (prime - 1) * prime ** exp)

	return product



class PAO_node:
	def __init__(self, prime, power, parent = False):
		self.prime = prime
		self.power = power
		self.children = []
		self.parent = parent

	def add(self, ch_prime, ch_power):
		newnode = PAO_node(ch_prime, ch_power, self)
		bisect.insort_left(self.children, newnode)
		#Return index of added element
		return bisect.bisect(self.children, newnode) - 1

	def __str__(self):
		return "{0}^{1}".format(self.prime, self.power)

	def __getitem__(self, index):
		return self.children[index]

	def __lt__(self, other):
		return b([(self.prime, self.power)]) > b([(other.prime, other.power)])
	
	def index(self, (prime, exp)):
		i = 0;
		for child in self.children:
			if child.prime == prime and child.power == exp:
				return i
			i+=1
		return -1


	def delta(self):
		prime = self.prime
		power = self.power
		return (prime ** power - 1) * prime / (1.0 * (prime ** (power + 1) - 1))

	def delta_up(self):
		prime = self.prime
		power = self.power
		return (prime ** (power + 2) - 1) / (1.0 * prime * (prime ** (power + 1) - 1))


	def delta_a(self):

		prime = self.prime
		power = self.power
		return (prime ** power - 1) * prime

	def delta_b(self):
		prime = self.prime
		power = self.power
		return (prime ** (power + 1) - 1)


	def primeset_above(self):
		root = self
		primelets = []
		while root.prime != 1:
			primelets.insert(0, (root.prime, root.power))
			root = root.parent

		return primelets

	def nodeset_above(self):
		root = self
		nodeset = []
		while root != False:
			nodeset.insert(0,root)
			root = root.parent

		return nodeset

	def abundance(self):
		primeset = self.primeset_above()
		return a(primeset) - 2 * prod(primeset) 

	def divisors_below_a(self):
		#Divisors to return:
		divisors = [1]
		abundance = self.abundance()
		#print abundance
		#print ("a")
		primelets = self.primeset_above()

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


def delta((prime, power)):
	return (prime ** power - 1) * prime / (1.0 * prime ** (power + 1) - 1.0)

def delta_plus((prime, power)):
	return (prime ** (power + 2) - 1) / (1.0 * prime * (prime ** (power + 1) - 1))


def abundance(primeset):
	return a(primeset) - 2 * prod(primeset) 

def divisors_below_a(primelets):
	#Divisors to return:
	divisors = [1]
	ab = abundance(primelets)
	#print ab


	for i in range(0, len(primelets)):
		p_product = 1
		prev_d_size = len(divisors)
		for exp in (1, primelets[i][1]):
			p_product *= primelets[i][0]
			#Do not allow larger than abundance
			if p_product > ab:
				break
			
			for div_index in range(0, prev_d_size):
				if (p_product * divisors[div_index] <= ab):
					#print "Appending divisors"

					divisors.append(p_product * divisors[div_index])
				else:
					break
	return divisors


#B - the beta function
#primelets - set of prime tuplets
#(PRIME, EXPONENT):
def b(primelets):
	product = 1
	for (prime, exp) in primelets:
		product *= (prime ** (exp + 1) - 1) / (1.0 * (prime - 1) * prime ** exp)

	return product


def b_l1(primelets):
	product = 1
	for (prime, exp) in primelets:
		product *= (prime ** (exp + 1) - 1)

	return product

def b_l2(primelets):
	product = 1
	for (prime, exp) in primelets:
		product *= (1.0 * (prime - 1) * prime ** exp)

	return product



#STEP 1
#Find the set of potential new prime `canadates` for abundant n:
def find(primeset):
	#Minprime is last prime on the list
	minprime = primeset[-1][0] + 1

	#Maxprime is calculated carefully:
	#FORMULA - b(n) b_inf(p) > 2
	#		   b(n) p / p - 1 > 2
	#		   (b(n) / 2) p > p - 1
	#          1 + (b(n) / 2) p > p
	#		   1 > p (1 - b(n) / 2)
	#          1 / (1 - b(n) / 2 ) > p

	maxprime = 1 / (1.0 -  (b(primeset) / 2.0))

	#primerange = sympy.primerange(minprime, maxprime + 1)
	#print set(primerange)
	return (minprime, maxprime)

#STEP 2
#Find the exponent for a potential prime - 
def find_exp(primeset, p):
	l_1 = b_l1(primeset)
	l_2 = b_l2(primeset)

	k_min = (log(1.0 * l_1) - log(1.0 * (l_1 * p - 2 * l_2 * (p - 1)) ) ) / log(p * 1.0)

	return int(ceil(k_min))

#lim_b - as n goes to infinity for a certain prime p
def lim_b(p):
	return p / (p - 1)


class PAOtree:
	def __init__(self):
		self.setup_tree()
		self.divisors = 3

	def setup_tree(self):
		self.root = PAO_node(1,1)
		
		self.root.add(3,5)
		#Here are the children of 3,5:
		child = self.root[0]
		child.add(5,2)
		child[0].add(13,1)
		
		self.root.add(3,4)
		#3,4
		child = self.root[1]
		child.add(5,3)
		child[0].add(13,1)
		child.add(5,2)
		child[1].add(13,2)
		'''
		self.root.add(3,3)
		#3,3
		child = self.root[2]
		child.add(5,3)
		child[0].add(13,2)
		child.add(5,2)
		child[1].add(11,1)
		child.add(5,1)
		child[2].add(7,1)
				
		self.root.add(3,2)
		#3,2
		child = self.root[3]
		child.add(5,2)
		child[0].add(7,1)
		child.add(5,1)
		child[1].add(7,2)
		'''

	def display(self):
		index = 0
		root = self.root
		self.rec_disp(root, index, 0)
		

	def add_power(self):
		root = self.root
		self.rec_power(root)
		self.divisors+=1


	def rec_power(self, node, deltalist = [], index = 0, tried_sets = set()):
		this_deltalist = deltalist[:]
		if node.prime != 1:
			this_deltalist.append(node.delta())
			if index == self.divisors:
			#	print 'From rec_power: {0}'.format(this_deltalist)
				self.make_sets(node.primeset_above(), this_deltalist, tried_sets)

		for child in node.children:
			self.rec_power(child, this_deltalist, index + 1, tried_sets)


	#make_sets - modify the set and attemt to branch
	def make_sets(self, pset, deltalist, tried_sets):
		#if pset in tried_sets:
		#	return
		#Break early for repeats
		set_index = 0
		new_pset = pset		

		for (prime, exp) in new_pset:
			if (exp != 1):
				#Try_down will call make sets at the end of the day
				self.try_down(pset, deltalist, set_index, tried_sets)
				
			set_index += 1
			


	def try_down(self, pset, deltalist, index, tried_sets):
		
		debug = False
		if pset[0] == (3,3):
			if pset[1] == (5,1):
				if pset[2] == (13, 2) or pset[2] == (13, 1):
					debug = True


		new_dlist = deltalist
		new_dlist.pop(index)
		primeset = []
		delta_plist = []
		#if debug:
		#	print "on"
		#	print index
		#	print new_dlist
		i = 0
		for (prime, power) in pset:
			#Append primeset ()
			if i == index:
				primeset.append((prime, power - 1))
				new_dlist.insert(index, delta((prime, power - 1)))
				delta_plist.append(delta_plus((prime, power)))
			else:
				primeset.append((prime,power))
				delta_plist.append(delta_plus((prime, power)))
			i += 1

		#print (primeset)
		#print new_dlist
		if tuple(primeset) in tried_sets:
			#Set previously tried
			return
		#Start at prime not modified - This will add exponents:
		self.attempt(tried_sets, primeset, new_dlist)
		self.up_set(primeset, delta_plist, tried_sets)
		self.make_sets(primeset, new_dlist, tried_sets)

	#Check higher exponents beyond the index
	def up_set(self, pset, delta_plist, tried_sets):	
		criteria = 2.0 / b(pset) 
		#print 'Criteria: {0}'.format(criteria)
		primeset = []
		i = 0

		for (prime, power) in pset:
			#Append primeset ()
			#print(delta_plus((prime,power)) )
			if delta_plus((prime,power)) < criteria:
				#print "Met "
				self.add_attempt(pset, i, tried_sets)
			i += 1

		#print (primeset)
		#print new_dlist
		if tuple(primeset) in tried_sets:
			#Set previously tried
			return

	def double_prime(self, pset, index, tried_sets):
		primeset = []

		for (prime, power) in pset:
			pass
			#Append all but  

	def prime_remove(self, pset, index, tried_sets):
		deltalist = []
		primeset = []

		i = 0
		for (prime, power) in pset:
			if i != index:
				primeset.append((prime, power))
				deltalist.append(delta((prime, power)))
			i += 1

	def add_attempt(self, pset, index, tried_sets):
		deltalist = []
		primeset = []

		i = 0
		for (prime, power) in pset:
			if i == index:
				primeset.append((prime, power +	1))
				deltalist.append(delta((prime, power + 1)))
			else:
				primeset.append((prime, power))
				deltalist.append(delta((prime, power)))
			i += 1
		#If it has been tried before:
		if tuple(primeset) in tried_sets:
			return

		self.attempt(tried_sets, primeset, deltalist)	

	def attempt(self, tried_sets, primeset, new_dlist):
		debug = False
		if primeset[0] == (3,3):
			if primeset[1] == (5,1):
				if primeset[2] == (13, 2) or primeset[2] == (13, 1):
					debug = True

		#TODO: THERE IS A BUG CAUSING THE DLIST TO BE INCORRECT. THIS IS 
		# A TEMPORARY AND BAD FIX
		new_dlist = []
		for pair in primeset:
			new_dlist.append(delta(pair))
			

		minprime = 0
		maxprime = 0
		
		#print "New set of primes"
		tried_sets.add(tuple(primeset))

		(minprime, maxprime) = find(primeset)

		new_primes = set(sympy.ntheory.generate.primerange(minprime,maxprime))

		for prime in new_primes:
			exp = find_exp(primeset, prime)
			newdelta = delta((prime, exp))
			maxdelta = max(new_dlist + [newdelta])
			#STEP 3
			#if debug:
			#	print "DEBUG:"
			#	print (new_dlist + [newdelta])
			#	print 'MAX: {0}'.format(maxdelta)
			

			#STEP 4 - SEE IF PRIMATIVE:			
			if (b(primeset) * b( [(prime, exp)] ) * maxdelta > 2):
				#print "Not a set: {0}".format(b(primeset) * b( [(prime, exp)] ) * maxdelta)
				pass
			else:
				self.add_branch(primeset, prime, exp)	

	def add_branch(self, primeset, nprime, nexp):
		print'primeset:{0} : prime,exp:{1}^{2}'.format(primeset, nprime, nexp)
		#if len(divisors_below_a(primeset + [(nprime, nexp)])) <= 43:
		#	if subset(divisors_below_a(primeset + [(nprime, nexp)]), abundance(primeset + [(nprime, nexp)])):
		#		print "not weird"
#		#		test = input()
		#	else:
		#		print "WEIRD: "
		#		test= input('Is this real?: ')
		root = self.root
		#Break and create children
		break_ = False
		for (prime, exp) in primeset:
			if root.index((prime,exp)) != -1 and break_ == False:
				#print "Repeat child"
				root = root.children[root.index((prime, exp))]
			else:
				break_ = True
				index_child = root.add(prime, exp)
				root = root.children[index_child]
				
				pass
				#index = root.add(root, exp )
		root.add(nprime, nexp)


	#recursive display
	def rec_disp(self, node, index, highn, f = open('numbers', 'w')):
		print "Putting ints into file"
		for i in range(0, index):
			f.write('    ')
		f.write(str(node))
		highn = max (prod(node.primeset_above()), highn)
		f.write( ": n = {0}\n".format(prod(node.primeset_above())) )

		

		if index == self.divisors:
			pass
			#print "n = {0}, abundant divisors: {1}".format(prod(node.primeset_above()), node.divisors_below_a())
		else:
			for child in node.children:
				self.rec_disp(child, index + 1, highn, f)


		return highn
