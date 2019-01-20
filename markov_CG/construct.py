import cplex
import pyomo
#from pyomo.opt import SolverFactory
#from pyomo.core import Var
from pyomo.environ import *

def MP(m):
	m.K = Set()
	m.H = Set()
	m.KH = Set(within=m.K*m.H)
	m.H_bar = Set()
	m.D = Set(within=m.K*m.H*m.H_bar)
	m.N = Set()

	m.c_hat = Param(m.K, m.H, default=0)
	m.finv_LO2 = Param(m.N, m.H_bar)
	m.finv_LN2 = Param(m.N, m.H_bar)
	m.V_LO2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.V_LN2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.c_LO2 = Param(m.N, default=[55, 237, 427, 621, 951])
	m.c_LN2 = Param(m.N, default=[50, 215, 388, 565, 864])
	m.pn_LO2 = Param(default=1000)
	m.pn_LN2 = Param(default=1000)
#NonNegativeReals
	m.z = Var(m.K, m.H, domain=NonNegativeReals, bounds=(0,1))
	m.z_bar = Var(m.H_bar, domain=NonNegativeReals, bounds=(0,1))
	#m.z = Var(m.K, m.H, domain=Boolean)
	#m.z_bar = Var(m.H_bar, domain=Boolean)
	m.x_LO2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.x_LN2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.rfinv_LO2 = Var(m.N, m.H_bar, domain=NonNegativeReals)
	m.rfinv_LN2 = Var(m.N, m.H_bar, domain=NonNegativeReals)

	def logic_1(m, k):#multiplier: u
		return sum(m.z[k,h] for h in m.H if (k,h) in m.KH) >= 1
	m.logic1 = Constraint(m.K, rule=logic_1)
	def logic_2(m):
		return sum(m.z_bar[h_bar] for h_bar in m.H_bar) >= 1
	m.logic2 = Constraint(rule=logic_2)
	def logic_3(m,k,h,h_bar):#multiplier: v
		return m.z_bar[h_bar] <= m.z[k,h]
	m.logic3 = Constraint(m.D, rule=logic_3)
	def logic_4(m,h_bar):#multipler: w
		return m.z_bar[h_bar] >= sum(m.z[k,h] for (k,h,h_bar_prime) in m.D if h_bar_prime == h_bar) - len(m.K) + 1
	m.logic4 = Constraint(m.H_bar, rule=logic_4)

	def logic_LO2(m):
		return sum(m.x_LO2[n] for n in m.N) >= 1
	m.logicLO2 = Constraint(rule=logic_LO2)
	def logic_LN2(m):
		return sum(m.x_LN2[n] for n in m.N) >= 1
	m.logicLN2 = Constraint(rule=logic_LN2)

	def inv_LO2_1(m, n, h_bar):
		return m.rfinv_LO2[n,h_bar] <= m.z_bar[h_bar]*m.finv_LO2[n,h_bar]
	m.invLO21 = Constraint(m.N, m.H_bar, rule=inv_LO2_1)
	def inv_LO2_2(m, n, h_bar):
		return m.rfinv_LO2[n,h_bar] <= m.x_LO2[n]*m.finv_LO2[n,h_bar]
	m.invLO22 = Constraint(m.N, m.H_bar, rule=inv_LO2_2)
	def inv_LO2_3(m, n, h_bar):
		return m.rfinv_LO2[n,h_bar] >= (m.z_bar[h_bar]+m.x_LO2[n]-1)*m.finv_LO2[n,h_bar]
	m.invLO23 = Constraint(m.N, m.H_bar, rule=inv_LO2_3)

	def inv_LN2_1(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] <= m.z_bar[h_bar]*m.finv_LN2[n,h_bar]
	m.invLN21 = Constraint(m.N, m.H_bar, rule=inv_LN2_1)
	def inv_LN2_2(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] <= m.x_LN2[n]*m.finv_LN2[n,h_bar]
	m.invLN22 = Constraint(m.N, m.H_bar, rule=inv_LN2_2)
	def inv_LN2_3(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] >= (m.z_bar[h_bar]+m.x_LN2[n]-1)*m.finv_LN2[n,h_bar]
	m.invLN23 = Constraint(m.N, m.H_bar, rule=inv_LN2_3)

	def netcost(m):
		return sum(m.c_hat[k,h]*m.z[k,h] for (k,h) in m.KH) + sum(m.x_LO2[n]*m.c_LO2[n] for n in m.N)+ sum(m.x_LN2[n]*m.c_LN2[n] for n in m.N)+ sum(m.rfinv_LO2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)+ sum(m.rfinv_LN2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)
	m.obj = Objective(rule=netcost, sense=minimize)
	return m

#---------------------------------------------------------------------------------------------------------------------

def MP1(m):
	m.K = Set()
	m.H = Set()
	m.KH = Set(within=m.K*m.H)
	m.H_bar = Set()
	m.D = Set(within=m.K*m.H*m.H_bar)
	m.N = Set()

	m.c_hat = Param(m.K, m.H, default=0)
	m.finv_LO2 = Param(m.N, m.H_bar)
	m.finv_LN2 = Param(m.N, m.H_bar)
	m.V_LO2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.V_LN2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.c_LO2 = Param(m.N, default=[55, 237, 427, 621, 951])
	m.c_LN2 = Param(m.N, default=[50, 215, 388, 565, 864])
	m.pn_LO2 = Param(default=1000)
	m.pn_LN2 = Param(default=1000)
#NonNegativeReals
	m.z = Var(m.K, m.H, domain=NonNegativeReals, bounds=(0,1))
	m.z_bar = Var(m.H_bar, domain=NonNegativeReals, bounds=(0,1))
	#m.z = Var(m.K, m.H, domain=Boolean)
	#m.z_bar = Var(m.H_bar, domain=Boolean)
	m.x_LO2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.x_LN2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.rfinv_LO2 = Var(m.N, m.H_bar, domain=NonNegativeReals)
	m.rfinv_LN2 = Var(m.N, m.H_bar, domain=NonNegativeReals)

	def logic_1(m, k):#multiplier: u
		return sum(m.z[k,h] for h in m.H if (k,h) in m.KH) >= 1
	m.logic1 = Constraint(m.K, rule=logic_1)
	#def logic_2(m):
	#	return sum(m.z_bar[h_bar] for h_bar in m.H_bar) >= 1
	#m.logic2 = Constraint(rule=logic_2)
	#def logic_3(m,k,h,h_bar):#multiplier: v
	#	return m.z_bar[h_bar] <= m.z[k,h]
	#m.logic3 = Constraint(m.D, rule=logic_3)
	#def logic_4(m,h_bar):#multipler: w
	#	return m.z_bar[h_bar] >= sum(m.z[k,h] for (k,h,h_bar_prime) in m.D if h_bar_prime == h_bar) - len(m.K) + 1
	#m.logic4 = Constraint(m.H_bar, rule=logic_4)

	def logic_LO2(m):
		return sum(m.x_LO2[n] for n in m.N) >= 1
	m.logicLO2 = Constraint(rule=logic_LO2)
	def logic_LN2(m):
		return sum(m.x_LN2[n] for n in m.N) >= 1
	m.logicLN2 = Constraint(rule=logic_LN2)

	def inv_LO2_1(m, n, k, h, h_bar):
		return m.rfinv_LO2[n,h_bar] <= m.z[k,h]*m.finv_LO2[n,h_bar]
	m.invLO21 = Constraint(m.N, m.D, rule=inv_LO2_1)
	def inv_LO2_2(m, n, h_bar):
		return m.rfinv_LO2[n,h_bar] <= m.x_LO2[n]*m.finv_LO2[n,h_bar]
	m.invLO23 = Constraint(m.N, m.H_bar, rule=inv_LO2_2)
	def inv_LO2_3(m, n, h_bar):
		return m.rfinv_LO2[n,h_bar] >= (sum(m.z[k,h] for (k,h,h_bar_prime) in m.D if h_bar_prime == h_bar) - len(m.K) +m.x_LO2[n])*m.finv_LO2[n,h_bar]
	m.invLO24 = Constraint(m.N, m.H_bar, rule=inv_LO2_3)

	def inv_LN2_1(m, n, k, h, h_bar):
		return m.rfinv_LN2[n,h_bar] <= m.z[k,h]*m.finv_LN2[n,h_bar]
	m.invLN21 = Constraint(m.N, m.D, rule=inv_LN2_1)
	def inv_LN2_2(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] <= m.x_LN2[n]*m.finv_LN2[n,h_bar]
	m.invLN22 = Constraint(m.N, m.H_bar, rule=inv_LN2_2)
	def inv_LN2_3(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] >= (sum(m.z[k,h] for (k,h,h_bar_prime) in m.D if h_bar_prime == h_bar) - len(m.K) +m.x_LN2[n])*m.finv_LN2[n,h_bar]
	m.invLN23 = Constraint(m.N, m.H_bar, rule=inv_LN2_3)

	def netcost(m):
		return sum(m.c_hat[k,h]*m.z[k,h] for (k,h) in m.KH) + sum(m.x_LO2[n]*m.c_LO2[n] for n in m.N)+ sum(m.x_LN2[n]*m.c_LN2[n] for n in m.N)+ sum(m.rfinv_LO2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)+ sum(m.rfinv_LN2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)
	m.obj = Objective(rule=netcost, sense=minimize)
	return m
