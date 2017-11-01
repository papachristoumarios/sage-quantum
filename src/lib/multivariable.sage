#!/usr/bin/sage -python

from sage.all import *
import sys; import numpy as np

# differential operators 
#TODO move to diff_op.sage

def laplacian(f):
	"""Returns the laplacian of a function f(x1,...,xn)"""
	vlist = f.variables()
	lapf = function('lapf')
	lapf = diff(f,vlist[0],2)
	for i in range(len(1,vlist)):
		lapf +=diff(f,vlist[i],2)
	return lapf
	
def divergence(F):
	"""Returns the divergence of a vector field F(x1,...,xn)"""
	vlist = get_vector_variables(F)
	assert (len(F) == len(vlist))
	n = len(F)
	divF = function('divF')
	divF = diff(F[0], vlist[0])
	for i in range(1,n):
		divF.add(diff(F[i], vlist[i]))
	return divF
	
def divergence_dot(F):
	"""Returns the divergence of a vector field F(x,y,z) using NumPy's dot(a,b) function"""
	vlist = get_vector_variables(F)
	assert (len(F) == len(vlist))
	n = len(F)
	for v in vlist:
		np.append(nabla, DiffOp(vlist[i])); np.append(farray, F[i])
	return np.dot(nabla,farray)
	
def curl_cross(F):
	"""Returns the curl of a field F(x,y,z) using NumPy's cross(a,b) function""" 
	assert(len(F)==3)
	vlist = get_vector_variables(F)
	nabla,  farray = np.array([]), np.array([])
	for v in vlist:
		np.append(nabla, DiffOp(vlist[i])); np.append(farray, F[i])
	return vector(np.cross(nabla,farray))

def curl2(F):
	"""Returns the curl of a field F(x,y,z)"""
	assert(len(F) == 3)
	if F is vector((0,0,0)):
		return vector((0,0,0))
	vlist = get_vector_variables(F)
	return vector([diff(F[2],vlist[1])-diff(F[0],vlist[2]), diff(F[0],vlist[2])-diff(F[2],vlist[0]), diff(F[1],vlist[0])-diff(F[0],vlist[1])]) 
	
def vector_desolve(F,G,n):
	"""Solve the vector equation F(n) = G(n)"""
	assert(len(F) == len(G))
	sols = []
	for i in range(len(F)):
		sols.append(desolve(F[i] == G[i], n))
	return sols
	
def vector_desolve_diff(F,G,v): 
	"""Solve differential equations for vector fields as F(x1,x2,...,xn) = G(x1,x2,...,xm)"""
	assert(len(F) == len(G)); n, R = len(F), []
	for i in range(n):
		f,g = function('f',v), function('g',v)
		f(v),g(v) = F[i], G[i]
		R.append(desolve(f(v) == g(v), g))
	return vector(R)
		
			
def get_vector_variables(F):
	n = len(F)
	vars = tuple(F[0].variables())
	for i in range(1,n):
		for v in F[i].variables():
			if not(v in vars): vars += (v)
	return vars
	
def get_characteristics_of_vector_field(F):
	return (divergence(F), curl2(F))
	
def multiple_integral(f, variables, bounds=None):
	"""Calculate the multiple integral of a function f with respect to different integrands"""
	g = function('g')
	g = f
	n = len(variables)
	if bounds is None:
		for i in range(n):
			g = integrate(g,variables[i])
	else:
		for i in range(n):
			assert(len(bounds[i]) == 2)
			g = integrate(g, variables[i], variables[i][0], variables[i][1])
	return g
	
class DiffOp:
"""Differential (derivative) operator that acts on multiplication"""

	def __init__(self, dep_var, order=1):
		self.dep_var = dep_var
		self.order = order
		
	def __mul__(self, f):
		if isinstance(f, DiffOp) is True:
			if f.dep_var is self.dep_var:
				return DiffOp(self.dep_var, self.order + f.order)
			else:
				return None
		else:
			return diff(f, self.dep_var, self.order)
			


#Field Parent Class

class SpaceField:
		def __init__(self, F, strength=vector((0,0,0)), time_dependent=True):
			assert(len(F) == 3)
			self.F, self.strength = F, strength
			self.vars, self.svars = get_vector_variables(F), get_vector_variables(strength)
			self.time_dependent = time_dependent
			if time_dependent: t = var('t')
			
		def __call__(self, x,y,z,t=0):
			if time_dependent:
				return self.F(x,y,z)
			else:
				return self.F(x,y,z,t)
			
		def __add__(self, other):
			return SpaceField(self.F + other.F, self.strength + other.strength, time_dependent=((self.time_dependent) or (other.time_dependent)))
		
				
#Maxwell's Equations			
			
class ElectricField(SpaceField):
	def __init__(self, E):
		assert (len(E) == 3)
		super(ElectricField, self).__init__(E)
		self.E = E
		self.vars = get_vector_variables(E)
		self.div, self.curl = get_characteristics_of_vector_field(E)
		
	def __call__(self, x,y,z):
		return self.E(x,y,z)
		
	def apply_Faradays_law(self, BF):
		if BF.time_dependent:
			Bt = -diff(BF.B, BF.t)
			sols = []
			for i in range(3):
				desolve(Bt[i] == self.curl, Bt[i],BF.t)
			return sols 
		else:
			eqns = []
			for i in range(3):
				eqns.append(self.curl[i] == 0)
			return solve(eqns, x,y,z)
			
class MagneticField(SpaceField):
	def __init__(self,B, H,time_dependent=True):
		assert(len(B) == 3 and len(H) == 3)
		self.B, self.H, self.vars = B, H, get_vector_variables(B)
		self.div, self.curl = get_characteristics_of_vector_field(B)
		self.time_dependent = time_dependent
		if time_dependent: self.t = var('t')
		
	def __call__(self, x,y,z,t=0):
		if not(self.time_dependent)
			return self.B(x,y,z)
		return self.B(x,y,z,t)

	def absence_of_magnetic_monopole(self, v, check=True):
		if check :
			if self.div is 0:
				return True
			return False
		else:
			solve(self.div == 0, v)
			
	def apply_maxwell_ampere_law(self, JF, DF):
		if not(DF.time_dependent):
			sols = []
			
			for i in range(3):
				sols.append(desolve(self.curl[i] == JF.J[i] + DF.Dt[i], Dt[i], DF.t))
		else:
		

class ElectricDisplacementField:
	def __init__(self, D, time_dependent=True):
		assert (len(D) == 3)
		self.D = D
		self.vars = get_vector_variables(D)
		self.div, self.curl = get_characteristics_of_vector_field(D)
		self.time_dependent = time_dependent
		if time_dependent: 
			self.t = var('t')
			self.Dt = self.D.diff(self.t)
			for i in range(3):
				self.D[i] = function('D[{0}]'.format(i), self.t)
				self.Dt[i] = function('Dt[{0}]'.format(i), self.t)
			
	def __call__(self, x,y,z,t=0):
		if not(self.time_dependent)
			return self.B(x,y,z)
		return self.B(x,y,z,t)
		
	def apply_gauss_law(self, v, rho=None, check=True):
		if rho is None:
			rho = var('rho')
			return solve(self.div == rho, v)
		elif rho is not None:
			return solve(self.div == rho, v)	
		elif check and rho !=None:
			if self.div is rho:
				return True
			return False
		else:
			raise Exception("Cannot apply Gauss's Law")
			
def apply_Faradays_law2(EF,BF):
	return EF.apply_Faradays_law(BF)

			
if __name__ == '__main__':
	print 'Tests'
	x = var('x')
	y = var('y')
	#f = function('f')
	f(x,y) = sin(x) + 2*x*y 
	gradf = diff(f)
	lapf = laplacian(f)
	print lapf
	show(lapf)

