import cplex
import pyomo
#from pyomo.opt import SolverFactory
#from pyomo.core import Var
from pyomo.environ import *

def MP_Indep(m):#initial independent optimization/equilibrium checking
	m.K = Set()
	m.H = Set()
	m.KH = Set(within=m.K*m.H)
	m.N = Set()

	m.c_hat = Param(m.KH, default=0)
	m.finv_LO2 = Param(m.N, m.KH)
	m.finv_LN2 = Param(m.N, m.KH)
	m.V_LO2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.V_LN2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.c_LO2 = Param(m.N, default=[55, 237, 427, 621, 951])
	m.c_LN2 = Param(m.N, default=[50, 215, 388, 565, 864])
#NonNegativeReals
	#m.z = Var(m.K, m.H, domain=NonNegativeReals, bounds=(0,1))
	m.z = Var(m.KH, domain=Boolean)
	m.x_LO2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.x_LN2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.rfinv_LO2 = Var(m.N, m.KH, domain=NonNegativeReals)
	m.rfinv_LN2 = Var(m.N, m.KH, domain=NonNegativeReals)

	def logic_1(m, k):#multiplier: u
		return sum(m.z[k,h] for h in m.H if (k,h) in m.KH) >= 1
	m.logic1 = Constraint(m.K, rule=logic_1)

	def logic_LO2(m):
		return sum(m.x_LO2[n] for n in m.N) >= 1
	m.logicLO2 = Constraint(rule=logic_LO2)
	def logic_LN2(m):
		return sum(m.x_LN2[n] for n in m.N) >= 1
	m.logicLN2 = Constraint(rule=logic_LN2)

	def inv_LO2_1(m, n, k, h):
		return m.rfinv_LO2[n,k,h] <= m.z[k,h]*m.finv_LO2[n,k,h]
	m.invLO21 = Constraint(m.N, m.KH, rule=inv_LO2_1)
	def inv_LO2_2(m, n, k, h):
		return m.rfinv_LO2[n,k,h] <= m.x_LO2[n]*m.finv_LO2[n,k,h]
	m.invLO22 = Constraint(m.N, m.KH, rule=inv_LO2_2)
	def inv_LO2_3(m, n, k, h):
		return m.rfinv_LO2[n,k,h] >= (m.z[k,h]+m.x_LO2[n]-1)*m.finv_LO2[n,k,h]
	m.invLO23 = Constraint(m.N, m.KH, rule=inv_LO2_3)

	def inv_LN2_1(m, n, k, h):
		return m.rfinv_LN2[n,k,h] <= m.z[k,h]*m.finv_LN2[n,k,h]
	m.invLN21 = Constraint(m.N, m.KH, rule=inv_LN2_1)
	def inv_LN2_2(m, n, k, h):
		return m.rfinv_LN2[n,k,h] <= m.x_LN2[n]*m.finv_LN2[n,k,h]
	m.invLN22 = Constraint(m.N, m.KH, rule=inv_LN2_2)
	def inv_LN2_3(m, n, k, h):
		return m.rfinv_LN2[n,k,h] >= (m.z[k,h]+m.x_LN2[n]-1)*m.finv_LN2[n,k,h]
	m.invLN23 = Constraint(m.N, m.KH, rule=inv_LN2_3)

	def netcost(m):
		return (sum(m.c_hat[k,h]*m.z[k,h] for (k,h) in m.KH) 
			+ sum(m.x_LO2[n]*m.c_LO2[n] for n in m.N)+ sum(m.x_LN2[n]*m.c_LN2[n] for n in m.N)
			+ sum(m.rfinv_LO2[n,k,h] for (n,k,h) in m.N*m.KH)+ sum(m.rfinv_LN2[n,k,h] for (n,k,h) in m.N*m.KH))
	m.obj = Objective(rule=netcost, sense=minimize)
	return m

def MP_extend(m):#extend based on current optimum
	m.H_bar = Set()
	m.N = Set()
	m.c_hat = Param(m.H_bar, default=0)
	m.finv_LO2 = Param(m.N, m.H_bar)
	m.finv_LN2 = Param(m.N, m.H_bar)
	m.c_LO2 = Param(m.N, default=[55, 237, 427, 621, 951])
	m.c_LN2 = Param(m.N, default=[50, 215, 388, 565, 864])
#NonNegativeReals
	#m.z = Var(m.K, m.H, domain=NonNegativeReals, bounds=(0,1))
	m.z_bar = Var(m.H_bar, domain=Boolean)
	m.x_LO2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.x_LN2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.rfinv_LO2 = Var(m.N, m.H_bar, domain=NonNegativeReals)
	m.rfinv_LN2 = Var(m.N, m.H_bar, domain=NonNegativeReals)
	def logic_1(m):#multiplier: u
		return sum(m.z_bar[h_bar] for h_bar in m.H_bar if h_bar in m.H_bar) >= 1
	m.logic1 = Constraint(m.H_bar, rule=logic_1)
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
		return m.rfinv_LO2[n,h_bar] >= (m.z_bar[h_bar] + m.x_LO2[n] - 1)*m.finv_LO2[n,h_bar]
	m.invLO23 = Constraint(m.N, m.H_bar, rule=inv_LO2_3)

	def inv_LN2_1(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] <= m.z_bar[h_bar]*m.finv_LN2[n,h_bar]
	m.invLN21 = Constraint(m.N, m.H_bar, rule=inv_LN2_1)
	def inv_LN2_2(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] <= m.x_LN2[n]*m.finv_LN2[n,h_bar]
	m.invLN22 = Constraint(m.N, m.H_bar, rule=inv_LN2_2)
	def inv_LN2_3(m, n, h_bar):
		return m.rfinv_LN2[n,h_bar] >= (m.z_bar[h_bar] + m.x_LN2[n] - 1)*m.finv_LN2[n,h_bar]
	m.invLN23 = Constraint(m.N, m.H_bar, rule=inv_LN2_3)

	def netcost(m):
		return sum(m.c_hat[h_bar]*m.z_bar[h_bar] for h_bar in m.H_bar) + sum(m.x_LO2[n]*m.c_LO2[n] for n in m.N)+ sum(m.x_LN2[n]*m.c_LN2[n] for n in m.N)+ sum(m.rfinv_LO2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)+ sum(m.rfinv_LN2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)
	m.obj = Objective(rule=netcost, sense=minimize)
	return m

'''
#-------------------------------------------------------------
def MP_extend(m):#extend based on current optimum
	m.K = Set()
	m.H = Set()
	m.KH = Set(within=m.K*m.H)
	m.HR = Set()
	m.KHR = Set(within=m.K*m.HR)
	m.KHHR = Set(within=m.K*m.H*m.HR)
	m.N = Set()

	m.c_hat = Param(m.K, m.H, default=0)
	m.finv_LO2 = Param(m.N, m.K, m.H, m.HR)
	m.finv_LN2 = Param(m.N, m.K, m.H, m.HR)
	m.V_LO2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.V_LN2 = Param(m.N, default=[100, 400, 700, 1000, 1500])
	m.c_LO2 = Param(m.N, default=[55, 237, 427, 621, 951])
	m.c_LN2 = Param(m.N, default=[50, 215, 388, 565, 864])
#NonNegativeReals
	#m.z = Var(m.K, m.H, domain=NonNegativeReals, bounds=(0,1))
	m.z = Var(m.K, m.H, m.HR, domain=Boolean)
	m.x_LO2 = Var(m.N, domain=Boolean)
	m.x_LN2 = Var(m.N, domain=Boolean)
	m.rfinv_LO2 = Var(m.N, m.K, m.H, m.HR, domain=NonNegativeReals)
	m.rfinv_LN2 = Var(m.N, m.K, m.H, m.HR, domain=NonNegativeReals)

	def logic_LO2(m):
		return sum(m.x_LO2[n] for n in m.N) >= 1
	m.logicLO2 = Constraint(rule=logic_LO2)
	def logic_LN2(m):
		return sum(m.x_LN2[n] for n in m.N) >= 1
	m.logicLN2 = Constraint(rule=logic_LN2)

	def inv_LO2_1(m, n, k, h, hr):
		return m.rfinv_LO2[n,k,h,hr] <= m.z[k,h,hr]*m.finv_LO2[n,k,h,hr]
	m.invLO21 = Constraint(m.N, m.KHHR, rule=inv_LO2_1)
	def inv_LO2_2(m, n, k, h, hr):
		return m.rfinv_LO2[n,k,h,hr] <= m.x_LO2[n]*m.finv_LO2[n,k,h,hr]
	m.invLO22 = Constraint(m.N, m.KHHR, rule=inv_LO2_2)
	def inv_LO2_3(m, n, k, h, hr):
		return m.rfinv_LO2[n,k,h,hr] >= (m.z[k,h,hr]+m.x_LO2[n]-1)*m.finv_LO2[n,k,h,hr]
	m.invLO23 = Constraint(m.N, m.KHHR, rule=inv_LO2_3)

	def inv_LN2_1(m, n, k, h, hr):
		return m.rfinv_LN2[n,k,h,hr] <= m.z[k,h,hr]*m.finv_LN2[n,k,h,hr]
	m.invLN21 = Constraint(m.N, m.KHHR, rule=inv_LN2_1)
	def inv_LN2_2(m, n, k, h, hr):
		return m.rfinv_LN2[n,k,h,hr] <= m.x_LN2[n]*m.finv_LN2[n,k,h,hr]
	m.invLN22 = Constraint(m.N, m.KHHR, rule=inv_LN2_2)
	def inv_LN2_3(m, n, k, h, hr):
		return m.rfinv_LN2[n,k,h,hr] >= (m.z[k,h,hr]+m.x_LN2[n]-1)*m.finv_LN2[n,k,h,hr]
	m.invLN23 = Constraint(m.N, m.KHHR, rule=inv_LN2_3)
	#def logic_1(m,k,h):
	#	return sum(m.z[k,h,hr] for hr in m.HR if (k,hr) in m.KHR) >= 1
	#m.logic1 = Constraint(m.KH, rule=logic_1)
	#def logic_2(m,k,hr):
	#	return sum(m.z[k,h,hr] for h in m.H if (k,h) in m.KH) >= 1
	#m.logic2 = Constraint(m.KHR,rule=logic_2)
	def logic_1(m, k):
		return sum(m.z[k,h,hr] for (h,hr) in m.H*m.HR if (k,h,hr) in m.KHHR) >= 1
	m.logic1 = Constraint(m.K, rule=logic_1)

	def netcost(m):
		return (sum(m.c_hat[k,h]*m.z[k,h,hr] for (k,h,hr) in m.KHHR)
		+ sum(m.x_LO2[n]*m.c_LO2[n] for n in m.N)+ sum(m.x_LN2[n]*m.c_LN2[n] for n in m.N)
		+ sum(m.rfinv_LO2[n,k,h,hr] for (n,k,h,hr) in m.N*m.KHHR )
		+ sum(m.rfinv_LN2[n,k,h,hr] for (n,k,h,hr) in m.N*m.KHHR ))
	m.obj = Objective(rule=netcost, sense=minimize)
	return m
'''