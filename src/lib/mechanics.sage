%load("./lib/multivariable.sage")

#planck constant / reduced planck constant
global hbar; hbar = var('hbar')
assume(hbar > 0)
global h; h = var('h')
assume(h>0)

#quantum physics
class Hamiltonian:
"""The hamiltonian operator of given potential function U(r)"""
	def __init__(self, U):
		self.potential = U
		
	def __mul__(self, psi):
		hbar = var('hbar')
		m = var('m')
		return -self.hbar/(2*self.m)*laplacian(psi) + self.potential

class PositionOperator:
	
	def __init__(self,v):
		self.v = v
	
	def __mul__(self, f):
		return f*v

class MomentumOperator:
	
	def __init__(self, ndim):
		self.ndim = ndim
		
	def __mul__(self, f):
		return -I*hbar*diff(f)
		
class Wavefunction:
	def __init__(self, ndim, mass = None, time_dependent=False, cartesian_notation=False):
		self.Psi = function('Psi')
		self.ndim = ndim
		self.time_dependent = time_dependent
		if time_dependent: self.t = var('t')
		if mass is None:
			self.mass = var('m')
			assume(m > 0)
		elif isinstance(mass, Integer) or isinstance(mass, Float):
			assert(m > 0)
			self.mass = mass
		self.vars = None
		if cartesian_notation:
			assert(ndim is 3); self.vars = (var('x'),var('y'),var('z'))
		else:
			self.vars = []
			for i in range(ndim):
				self.vars.append(var('x{0}'.format(i)))
			self.vars = tuple(self.vars)
		
	
	def __call__(self, values):
		assert(len(values) is self.ndim)
		return self._psi(*values)
			
	def solve_tise(self, V, E):
		assert (len(V.variables()) == len(E.variables()) and len(E.variables()) == self.ndim)
		_H = Hamiltonian(V)
		self.Psi(*self.vars) = 0
		for i,x in enumerate(self.vars):
			psi = function('psi',x)
			se = _H*psi == E*psi
			psi = desolve(se,psi,ivar=x)
			self.Psi = self.Psi + var('a_{0}'.format(i))*psi
		
	def solve_tdse(self,V,E):
		assert (len(V.variables()) == len(E.variables()) and len(E.variables()) == self.ndim + 1)
		self.Psi(*self.vars) = 0
		_H = Hamiltonian(V)
		for i,x in enumerate(self.vars):
			f = function('f',x)
			g = function('g',self.t)
			g = desolve (diff(g,t) == E/(I*hbar)*g(t), g, ivar=var('t'))
			f = desolve (_H*f = E*f, f, ivar=x)
			self.Psi = self.Psi + var('a_{0}'.format(i))*f(x)*g(t)
			
	def get_normalization_condition(self, bounds):
		assert(bounds is not None)	
		try:
			p = multiple_integral(self.Psi.norm(), self.vars, bounds=bounds)
			self.normalization_condition = p == 1
			return p == 1
		except ValueError:
			self.normalization_condition = None
			
