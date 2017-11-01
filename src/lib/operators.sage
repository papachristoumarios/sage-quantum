import numpy as np
from matplotlib import pyplot as plt

check_normalization = lambda _M : isinstance(inner_product(_M, _M), sage.symbolic.expression.Expression) or ((inner_product(_M,_M)) == 1)
inner_product = lambda A,B: (A.T.conjugate()*B)[0,0]

class Operator(object):
	
	def __init__(self, vars = None, expr=None):
		self.expr = expr
		if not(vars == None):
			var(*vars)
		self.vars = list(vars)
		self.is_matrix_operator = hasattr(self, 'M') # 'M' represents the corresponding matrix of the operator
		
	def __call__(self, f):
		pass
		
	def __mul__(self, f):
		if not isinstance(f, Operator):
			return self.__call__(f)
		#recheck
		elif isinstance(f, Operator):
			g = define_arbitrary_function([self, f])
			return self.__call__(f*g)/g
			
	def __repr__(self):
		f = define_arbitrary_function([self])
		return self.__call__(f)
	
	def __power__(self, N, f):
		assert (N>=1 and N == int(N)) 
		_s = self.__call__(f)
		for i in range(N-1):
			_s = self.__call__(_s)
		return _s
	
	def eigenvalues(self):
		#solve eigenvalue equation
		#TODO add exceptions
		if self.is_matrix_operator:
			try:
				return self.M.eigenvalues()
			except (ValueError, TypeError, NotImplementedError):
				pass
				
		a = var('a'); function('f')(*self.vars)
		if len(self.vars) is 1:
			eigenvalue_equation = self.__call_(f(*self.vars)) == a*f(*self.vars)
			try:
				g(*self.vars) = desolve_laplace(eigenvalue_equation, f(*self.vars), ivar=self.vars[0])			
				eigenvalue_equation = eigenvalue_equation.substitute_function(f,g)
				result =  (eigenvalue_equation.factor()/g(x)).simplify()
				if not(self.vars[0] in result.variables()):
					return solve(result, a)
				else:
					raise NotImplementedError
			except NotImplementedError:
				pass
		else:
			raise NotImplementedError
		return None

class VectorOperator(Operator):
	
	def __init__(self, directions):
		super(VectorOperator, self).__init__(directions)
		self.components = np.array([])
	
	def __len__(self):
		return len(self.components)
		
	def __call__(self, f):
		if isinstance(f, np.ndarray):
			assert (len(f) is self.__len__())
			result = np.array([])
			for i in range(len(f)):
				np.append(result, self.components[i](f[i]))
			return result
		elif len(f) is 1:
			result = np.zeros(len(self.components))
			for i in range(len(result)):
				result[i] = self.components[i](f)
				return result
		else:
			pass
			
	def __mul__(self, f):
		return self.__call__(f)

	def cross_product(self, other):
		return cross_product_operators(self,other)
	
def define_arbitrary_function(list_of_operators):
	if len(list_of_operators) is 1:
		return function('f')(var(list_of_operators[0]))
	else:
		l = []	
		for O in list_of_operators:
			l += O.vars
		l = list(set(l))
		return function('f')(*var(l)) 
	
def commute(A,B, g = None):	
	if g is None:
		f = define_arbitrary_function([A,B])
	else:
		f = g
	eq = A(B(f)) - B(A(f))
	#print eq
	try:
		eq = eq.factor()/f
	except ValueError:
		pass 
	finally:	
		return eq

def anticommute(A,B, g=None):
	if g is None:
		f = define_arbitrary_function([A,B])
	else:
		f = g
	eq = A(B(f)) + B(A(f))
	try:
		eq = eq.factor()/f
	except ValueError:
		pass
	finally:
		return eq

def get_uncertainty(A,B, g=None):
	rel = commute(A,B, g)/2
	return sqrt(rel.norm()).simplify()
	
def cross_product_operators(A,B, g=None):
	assert(len(A) == len(B) and len(A) == 3)
	f = g
	C = np.cross(A*f, B*f)
	try:
		for c in C:
			c = c.factor()
	except ValueError:
		pass
	return C

class PartialDifferentialOperator(Operator):
	
	def __init__(self, direction):
		super(PartialDifferentialOperator, self).__init__([direction])
		self.direction = direction
		
	def __call__(self,f):
		v = var(self.direction)
		return diff(f, v) 	
	
class PositionOperator(Operator):
# the \hat x operator
	def __init__(self, direction):
		super(PositionOperator, self).__init__([direction])
		self.direction = direction
		
	def __call__(self, f):
		v = var(self.direction)
		return v*f
		
class MomentumOperator1D(Operator):
# the \hat p_x, p_y, p_z operator (linear momentum)
	def __init__(self, direction):
		super(MomentumOperator1D, self).__init__([direction])
		self.direction = direction
		
	def __call__(self, f):
		var('h')
		assume(h>0)
		return I*h/(2*pi)*f.diff(var(self.direction))
		
class SpatialMomentumOperator(VectorOperator):
#linear momentum operator as a vector	
	def __init__(self, directions):
		super(SpatialMomentumOperator, self).__init__(directions)
		self.components = np.array([])
		for direction in directions:
			np.append (self.components, MomentumOperator1D(direction))
			
class KineticEnergyOperator(Operator):
	
	def __init__(self, directions, mass=None):
		self.mass = mass
		super(KineticEnergyOperator, self).__init__(directions)
		self.p = SpatialMomentumOperator(directions)
		if self.mass is None:
			self.mass = var('m')
		
	def __call__(self, f):
		return self.p(self.p(f))/(2*mass)

class HamiltonianOperator(Operator):

	def __init__(self, directions, V, mass=None):
		self._T = KineticEnergyOperator(directions, mass)
		self._V = PotentialEnergyOperator(V)

	def __call__(self, f):
		return self._T(f) + self._V(f)

class PotentialEnergyOperator(Operator):

	def __init__(self, U):
		self.U = U
		#self.separable = self.is_separable()

	# TODO add separability conditions

	def is_separable(self):
		raise NotImplementedError
		separability_conditions = [




		]
		result = False

		for condition in separability_conditions:
			result = result or condition
		return result

#quantum computing

class BasicQuantumStates:
	ZERO = QuantumState(1, [1,0])
	ONE = QuantumState(1, [0,1])
	PLUS = QuantumState(1, [1/sqrt(2), 1/sqrt(2)])
	MINUS = QuantumState(1, [1/sqrt(2), -1/sqrt(2)])
	BELL = staticmethod(lambda x,y: BasicQuantumGates.CNOT(BasicQuantumGates.HADAMARD1(QuantumState(1, [kronecker_delta(x,0), kronecker_delta(x,1) ] )), QuantumState(1, [kronecker_delta(y,0), kronecker_delta(y,1)])))
	class BellStates:
		def __getitem__(self, (x,y)):
			return BasicQuantumStates.BELL(x,y)
	
class BasisSpace:
	
	def __init__(self, k, basis_dict = {}):
		self.k = k
		self.matrix_space = MatrixSpace(SR, k,1)
		self.basis = {}
		if basis_dict == {}:
			for i,basis in enumerate(self.matrix_space.basis()):
				elem = '{0}{1}'.format( ceil(N(log(k,2)))*'0' - len(bin(i)[2:])  , bin(i)[2:])
				self.basis[elem] = basis
		else:
			self.basis = basis_dict
			
	@property		
	def is_orthogonal(self): return self._is_orthogonal
	
	@is_orthogonal.getter
	def is_orthogonal(self):
		self._is_orthogonal =  orthogonality_test(self.basis)
		return self._is_orthogonal
	
	@staticmethod
	def orthogonality_test(d):	
		result = True
		K = len(d.values())
		for i in range(K):
			for j in range(i+1, K-1):
				p = inner_product(d[i], d[j])
				result = self._is_orthogonal and (isinstance(p, sage.symbolic.expression.Expression) or p == 0)			
		return result
		  
	def get_arbitrary_element(self):
		vars = [var('v{0}'.format(i)) for i in range(self.k)]
		return prod([vars[i]*self.basis.values()[i] for i in range(self.k)])
		
	def __pow__(self, n):
		return BasisSpace(self.k ^ n)
	
	def __getitem__(self, key):
		if isinstance(key, Integer):
			return self.basis[self.basis.keys()[key]]
		return self.basis[key]
	
	def __setitem__(self, key, value):
		if isinstance(key, Integer):
			self.basis[self.basis.keys()[key]] = value
		self.basis[key] = value
		
	def switch(self, D, *args):
		assert(len(D) == self.k)
		new_basis = BasisSpace(self.k)
		new_basis.basis = {}
		
		for key in D.keys():
			new_basis[key] = D[key](*args) 
		assert (new_basis.is_orthogonal)
		
		#self = new_basis; 
		return new_basis
				 
	def switch_to_signum_basis(self):
		assert(self.k == 2)
		self.switch({'+': lambda: 1/sqrt(2)*(self[0] + self[1]), '-': lambda: 1/sqrt(2)*(self[0]) - self[1]})
			
	def rotate(self, f):
		assert(self.k == 2)
		rotation_matrix = matrix(2,2, [cos(f), -sin(f), sin(f), cos(f)])
		for key in self.basis.keys():
			self.basis[key] = rotation_matrix*self.basis[key] 
	
	def apply_fcn(self, f, *args):
		for key in self.basis.keys():
			self.basis[key] = f(self.basis[key], *args)
			
	def __eq__(self, other):
		return (self.basis == other.basis)
		
	def __len__(self, other):
		return len(self.basis) == len(other.basis)
		
	def __mul__(self, other):
		if isinstance(other, QuantumState):
			return QuantumState(other.number_of_qubits, [other[i,0]*self.__getitem__[i] for i in range(2^other.number_of_qubits)])
		elif isinstance(other, BasisSpace):
			new_basis_dict = {}
			for i in range(self.__len__()):
				for j in range(len(other)):
					new_basis_dict['{0}{1}'.format(self.basis.keys()[i], other.basis.keys()[j])] = self.__getitem__[i].tensor_product(other.__getitem__[j])
			return BasisSpace(len(new_basis_dict), basis_dict=new_basis_dict)
	
class BasicComputationalBases:	
	SIGNUM_BASIS = BasisSpace(2, basis_dict = {'+': matrix(2,1, [1/sqrt(2), 1/sqrt(2)]), '-': matrix(2,1, [1/sqrt(2), -1/sqrt(2)]) })
	STANDARD_BASIS = BasisSpace(2)	
							
class QuantumState:
	
	def __init__(self, number_of_qubits=1, coefficients_array = None):
		assert(number_of_qubits >= 1)
		self.number_of_qubits = number_of_qubits
		self.computational_bases = BasisSpace(2^number_of_qubits)
		def _generate_computational_bases_and_coefficients(k, generate_variables=True):
			coefficients = []
			for j in range(2^k):
				_s = 'c{0}{1}'.format((k - len(bin(j)) + 2) * '0',str(bin(j))[2:])
				if generate_variables:
					__s = var(_s); coefficients.append(__s)
			coefficients = matrix(coefficients)
			if not(generate_variables):
				coefficients = zero_matrix(1,2^k)[0]
			return matrix(coefficients).T
		self.computational_coefficients = _generate_computational_bases_and_coefficients(number_of_qubits, coefficients_array is None)
		if not(coefficients_array is None): self.computational_coefficients = matrix(coefficients_array).T
			
	def __eq__(self, other):
		assert(isinstance(other, QuantumState))
		return (self.number_of_qubits == other.number_of_qubits) and (self.computational_bases == other.computational_bases) and (self.computational_coefficients == other.computational_coefficients)

	def __repr__(self):
		return self.computational_coefficients.__repr__()

	def tensor_product(self, other):
		return self.__mul__(other)

	def __mul__(self, other):
		#add tensor product
		result_computational_coefficients = self.computational_coefficients.tensor_product(other.computational_coefficients)
		result = QuantumState(self.number_of_qubits+other.number_of_qubits)
		result.computational_coefficients = result_computational_coefficients
		return result
	
	def inner_product(self, other):
			return (self.computational_coefficients.conjugate().T * other.computational_coefficients)[0,0]
			
	def __pow__(self, power):
		assert(power >= 1 and isinstance(power, int))
		return prod([self for _i in range(power)])
	
	@property
	def is_normalized(self): return self._is_normalized
	
	@is_normalized.getter
	def is_normalized(self):	
		return check_normalization(self.computational_coefficients)
		
	def __getitem__(self, (I,J)): return self.computational_coefficients[I,J]
	
	def __setitem__(self, (I,J), value):
		backup = self.computational_coefficients
		self.computational_coefficients[I,J] = value
		if not(self.is_normalized): self.computational_coefficients = backup
		
	@property
	def nrows(self):
		return self.computational_coefficients.nrows()
		
	@nrows.getter
	def nrows(self):
		return self.computational_coefficients.nrows()
	
	@property
	def ncols(self):
		return self.computational_coefficients.ncols()
		
	@ncols.getter
	def ncols(self):
		return self.computational_coefficients.ncols()
	
	def decompose(self):
		if self.number_of_qubits == 1:
			return list(self)
		
		n = self.nrows / 2
		vars, states, result = [], [], {}
		for i in range(n):
			vars.append(var('a{0}'.format(i)))
			vars.append(var('b{0}'.format(i)))
			states.append(QuantumState(1, [vars[2*i], vars[2*i + 1] ]))
		L = prod(states)
		eqns = [L[i,0] == self[i,0] for i in range(2*n)]
		for i in range(n):
			eqns.append(vars[2*i]*vars[2*i].conjugate() + vars[2*i + 1]*vars[2*i + 1].conjugate() == 1)
		
		solns = solve(eqns, *vars, solution_dict=True)	
		if solns == {}: raise DecompositionError('State is entangled')
		for i,soln in enumerate(solns):
			soln_values = soln.values()
			result[i-1] = [QuantumState(1, [soln_values[2*i], soln_values[2*i+1]]) for i in range(n)]
		return result	
			
	@property
	def probability_matrix(self):
		return self._probability_matrix		
			
	@probability_matrix.getter
	def probability_matrix(self):
		self._probability_matrix = self.computational_coefficients
		for _i in range(self.ncols):
			for _j in range(self.nrows):
				self._probability_matrix[_i,_j] = self._probability_matrix[_i,_j] * self._probability_matrix[_i,_j].conjugate()		
		return self._probability_matrix
		
	def switch_computational_basis(self, basis):
		assert (len(basis) == len(self.computational_bases))
		self.computational_bases = basis
		for i in range(2^self.number_of_qubits):
			self.computational_coefficients[i,0] = self.computational_coefficients[i,0]*self.computational_bases.basis.values()[i]
			
	def probability(self, J):
		return inner_product(self.computational_bases[J], self.computational_coefficients).norm()
		
	@property
	def is_entangled(self): return self._is_entangled
	
	@is_entangled.getter
	def is_entangled(self):
		try:
			self.decompose()
		except DecompositionError:
			return True
		return False	
			
class DecompositionError(Exception):
	"""Is raised if state is entangled"""
	
	def __init__(self): super(DecompositionError, self).__init__()		
			
class BasicQuantumGates:
	class PauliMatrices:
		_PAULIX = matrix(2,2, [0,1,1,0])
		_PAULIY = matrix(2,2, [0,-I,I,0])
		_PAULIZ = matrix(2,2, [1,0,0,-1])
	PAULIX = QuantumGate(PauliMatrices._PAULIX)
	PAULIY = QuantumGate(PauliMatrices._PAULIY)
	PAULIZ = QuantumGate(PauliMatrices._PAULIZ)
	NOT = PAULIX
	SWAP = QuantumGate(matrix(4,4, [1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1]))		
	CNOT = ControlledGate(NOT)
	CONTROLLEDX = CNOT
	CONTROLLEDY = ControlledGate(PAULIY)
	CONTROLLEDZ = ControlledGate(PAULIZ)
	phi_ = var('phi') 
	PHASESHIFT = staticmethod(lambda phi_=phi: QuantumGate(matrix(2,2, [1,0,0,exp(phi_*i)])))
	CONTROLLEDPHASESHIFT = ControlledGate(PHASESHIFT)
	SQUAREROOTSWAP = QuantumGate(matrix(4,4, [1,0,0,0,0,1/2*(1 + I), 1/2*(1-I), 0, 0, 1/2*(1-I), 1/2*(1+I), 0, 0,0,0,1]))
	TOFFOLI = ControlledGate(NOT, iters=2)
	ISWAP = QuantumGate(matrix(4,4, [1,0,0,0,0,0,I,0,0,I,0,0,0,0,0,1]))
	SIGN = QuantumGate(matrix(2,2, [1,0,0,-1]))
	CSIGN = ControlledGate(SIGN)
	DEUTSCH = staticmethod(lambda phi=phi_: ControlledGate(QuantumGate(matrix(2,2, [I*cos(phi), sin(phi), sin(phi), I*cos(phi)]))))
	HADAMARD1 = HadamardGate(1)
	
	@staticmethod
	def QFT_Matrix(n):
		omega = exp(2*pi*I/n)
		M = matrix(SR, n,n)
		for k in range(M.nrows()):
			for l in range(M.ncols()):
				M[i,j] = omega^(k*l)
		return 1/sqrt(n) * M

	QFT = staticmethod(lambda n=1: QuantumGate(QFT_Matrix(n)))

class QuantumGate(Operator):
	
	@property
	def applicable_qubits(self):
		return self._applicable_qubits
		
	@applicable_qubits.setter
	def applicable_qubits(self, qubits):
		assert(self._applicable_qubits is None or 2^len(qubits) is self.M.nrows())
		self._applicable_qubits = qubits	

	def __init__(self, M):
		self._applicable_qubits = None
		self.M = M
		def _check_if_unitary(X):
			assert (X.nrows() is X.ncols())
			return (X.T.conjugate()*X) is (X.inverse()*X)
		self.is_unitary = _check_if_unitary(self.M)
		try:
			vars = self.M.variables()
		except AttributeError:
			vars = None
		finally:
			super(QuantumGate, self).__init__(vars=vars, expr=self.M)
			
		
	def __call__(self, *states):
		assert(isinstance(prod(states), QuantumState))
		if len(states) == 1:
			states[0].computational_coefficients = self.__mul__(states[0])
			return states[0]
		else:
			state = prod(states)
			return self.__call__(state)

	def __add__(self, other):
		return QuantumGate(self.M + other.M)
	
	def __sub__(self, other):
		return QuantumGate(self.M - other.M)
					
	def __mul__(self, state):
		#if self._applicable_qubits is None:
		assert(isinstance(state, QuantumState))
		return self.M * state.computational_coefficients 

	def __eq__(self,other):
		assert(isinstance(other, QuantumGate))
		return self.M is other.M
		
	def __repr__(self):
		return self.M.__repr__()
		
	def __len__(self):
		return self.M.nrows()
		
	def show(self):
		self.__repr__()
		
	@property
	def inverse(self):
		return self._inverse
	
	@inverse.getter
	def inverse(self):
		try:
			self._inverse = self.M.inverse()
			return QuantumGate(self._inverse)
		except ZeroDivisionError:
			print 'Gate is singular'
			self._inverse = None
			
def ControlledMatrix(M, iters=1):
	if iters == 1:
		C = matrix(SR, 2*M.nrows(), 2*M.ncols(), 1)
		f = lambda q: q + M.nrows()
		for i in range(M.nrows()):
			for j in range(M.ncols()):
				C [f(i),f(j)] = M[i,j]
		return C
	elif iters > 1:
		C = ControlledMatrix(M) 
		for i in range(iters-1):
			C = ControlledMatrix(C)
		return C

ControlledGate = lambda U,iters=1: QuantumGate(ControlledMatrix(U.M, iters=iters))
HadamardGate = lambda  k: QuantumGate((1/2^(k/2)) * sage.combinat.matrices.hadamard_matrix.hadamard_matrix(2^k))
SimpleQuantumState = lambda coefficients_array=None: QuantumState(1, coefficients_array=coefficients_array)
		
class ClassicalBit:
	ZERO = 0
	ONE  = 1
	
	@property
	def value(self):
		return self._value
		
	@value.setter
	def value(self, v):
		assert (v is ClassicalBit.ZERO or v is ClassicalBit.ONE)
		self._value = v
		
class QuantumCircuitExceptions:

	class EmptyCircuitException(Exception):
		def __init__(self):
			pass

#class QuantumCircuit:

	#def __init__(self, states_array):
		#assert(all(isinstance(states_array_elem, SimpleQuantumState) for states_array_elem in states_array))
		#self.quantum_register_array = states_array
		#self.initial_quantum_register_array = states_array
		#self.quantum_gates = []

	#def add_gate(self, U, *argv):
		##assert(isinstance(U, QuantumGate) and (2^state.number_of_qubits is U.cols())) 
		#assert((isinstance(U, QuantumGate) or isinstance(U, ControlledGate)) and (2^len(argv) is U.cols()) and all(counter < len(states_array) for counter in range(len(states_array))) and len(set(list(argv))) is len(list(argv)))
		#U.applicable_qubits = list(argv)
		#self.quantum_gates.append(U)
	
	#@staticmethod
	#def apply_gate(G):
		
			
		
	#def simulate(self, steps=len(self.quantum_gates)):
		#assert(steps <= len(self.quantum_gates))
		#if self.quantum_gates is not []:
			#for j in range(steps):
				
				
		
			#for j in range(steps):
				#self.quantum_register = self.quantum_gates[i](self.quantum_register)
		
		
		#else:
			#raise QuantumCircuitExceptions.EmptyCircuitException()

	#def revert(self):
		#self.quantum_register = self.initial_quantum_register

	#def is_symbolic(self):
		#return isinstance(sum(self.quantum_register.computational_coefficients), sage.symbolic.expression.Expression)

	#def is_numeric(self):
		#return not(self.is_symbolic())

	#def get_probability_matrix(self):
		#probability_matrix = self.quantum_register
		#for _x in self.quantum_register:
			#for _y in _x:
				#_y = _y.norm()
		#return probability_matrix

	#def plot(self, probability=True):
		#fig = plt.figure(); ax = fig.add_subplot(111)
		#if probability:
			#Y_data = self.get_probability_matrix()
			#X_data = []
			#for _i in len(Y_data):
				#for _j in len(Y_data[0])
					#X_data.append('{0}{1}'.format(_i,_j))
			#Y_data = Y_data.list()
			#N = len(probability_matrix)^2
			#ind = np.arange(N)
			#rects1 = ax.bar(ind, menMeans, width,
				#color='black',
				#yerr=menStd,
				#error_kw=dict(elinewidth=2,ecolor='red'))

class QuantumCircuit:
	
	def __init__(self):
		self.registers = {}
		
	def add_register(self, name, noq):
		self.registers[name] = QuantumRegister(noq)
		
	def __len__(self):
		return len(self.registers.keys())
	
	def __getitem__(self, key):
		return self.registers[name]
		
	def __setitem__(self, key, value):
		assert(isinstance(value, QuantumRegister))
		self.registers[key] = value
	
class QuantumRegister:
	
	def __init__(self, noq, initial_state=BasicQuantumStates.ZERO):
		assert(isinstance(initial_state, QuantumState))
		self.noq = noq
		self.qubits = []
		for j in range(noq):
			self.qubits.append(initial_state)
		self.state = self.get_state()
		self.gates = GatesDictionary(noq)
						
	def get_state(self):
		self.state = prod(self.qubits)
		return self.state
		
	def add_gate(self, U, *args):
		self.gates.add(U, *args)
		
	def __getitem__(self, key):
		return  self.qubits[key]
		
	def __setitem__(self, key, value):
		assert (check_normalization(value))
		self.qubits[key].computational_coefficients = value
		
	def __repr__(self):
		self.get_state(); return self.state.__repr__()
			
	def __add__(self, other):
		assert(isinstance(other, QuantumRegister))
		result = QuantumRegister(self.noq + other.noq)
		for i in range(len(result.qubits)):
			if i < self.noq:
				result.qubits[i] = self.qubits[i]
			else:
				result.qubits[i] = other.qubits[i - self.noq + 1]
		return result	
		
class GatesDictionary:
	#dict-type class used for single steps in a quantum circuit
	def __init__(self, n):
		self.d = {}
		self.values = []
		
	def add(self, U, *args):
		for arg in args:
			assert ((arg <= n-1) and not(arg in self.values))
			self.values.append(arg)
		try:
			self.d[U]
		except KeyError:
			self.d[U] = []
		finally:	
			self.d[U].append(args)
		
	def __getitem__(self, u): return self.d[u]
	
	def __setitem__(self, U, args): self.add(U, *args)
	
	
	
	
	
	
	
	
	
	
	
	
	

		

				
if __name__ == '__main__':
	

	#an example case
	px = MomentumOperator1D('x')
	py = MomentumOperator1D('y')
	pz = MomentumOperator1D('z')
	X, Y, Z = PositionOperator('x'), PositionOperator('y'), PositionOperator('z')
	var('x,y,z')
	f  = function('f')(x,y,z) #arbitrary function
	r = np.array([X,Y,Z]) #position operator as a vector
	R = np.array([x,y,z]) #spatial variables
	p = np.array([px,py,pz]) #vector linear momentum operator
	L = cross_product_operators(r,p, f(x,y,z)) #angular momentum operator
