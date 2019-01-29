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
	m.x_LO2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.x_LN2 = Var(m.N, domain=NonNegativeReals, bounds=(0,1))
	m.rfinv_LO2 = Var(m.N, m.H_bar, domain=NonNegativeReals)
	m.rfinv_LN2 = Var(m.N, m.H_bar, domain=NonNegativeReals)

	def logic_1(m, k):#multiplier: u
		return sum(m.z[k,h] for h in m.H if (k,h) in m.KH) == 1
	m.logic1 = Constraint(m.K, rule=logic_1)
	def logic_2(m):
		return sum(m.z_bar[h_bar] for h_bar in m.H_bar) == 1
	m.logic2 = Constraint(rule=logic_2)
	def logic_3(m,k,h,h_bar):#multiplier: v
		return m.z_bar[h_bar] <= m.z[k,h]
	m.logic3 = Constraint(m.D, rule=logic_3)
	def logic_4(m,h_bar):#multipler: w
		return m.z_bar[h_bar] >= sum(m.z[k,h] for (k,h,h_bar_prime) in m.D if h_bar_prime == h_bar) - len(m.K) + 1
	m.logic4 = Constraint(m.H_bar, rule=logic_4)

	def logic_LO2(m):
		return sum(m.x_LO2[n] for n in m.N) == 1
	m.logicLO2 = Constraint(rule=logic_LO2)
	def logic_LN2(m):
		return sum(m.x_LN2[n] for n in m.N) == 1
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
		return sum(m.c_hat[k,h]*m.z[k,h] for (k,h) in m.KH) 
		+ sum(m.x_LO2[n]*m.c_LO2[n] for n in m.N)
		+ sum(m.x_LN2[n]*m.c_LN2[n] for n in m.N)
		+ 3650*m.pn_LO2*sum(m.rfinv_LO2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)
		+ 3650*m.pn_LN2*sum(m.rfinv_LN2[n,h_bar] for (n,h_bar) in m.N*m.H_bar)
	m.obj = Objective(rule=netcost, sense=minimize)
	return m

def PP(m):
	m.H = Set()
	m.S = Set()
	m.SF = Set(within=m.S)
	m.HS = Set(within=m.H*m.S)

	m.inst = Param(m.H)
	m.Q = Param(m.S*m.S)
	m.M = Param(default=1000)#the big M
	m.u = Param(default=1000)
	m.v = Param(default=1000)

	m.z = Var(m.H, domain=Binary)
	m.c = Var(domain=NonNegativeReals)
	m.pi = Var(m.S, domain=NonNegativeReals, bounds=(0,1))

	def Cost(m):
		return c >= sum(m.z[h]*m.inst[h] for h in m.H)
	m.cost = Constraint(rule=Cost)
	def logic_1(m):
		return sum(m.z[h] for h in m.H)
	m.logic1 = Constraint(rule=logic_1)
	def logic_2(m,s):
		return m.pi[s] <= sum(m.z[h] for (h,s) in m.HS)
	m.logic2 = Constraint(m.S, rule=logic_2)
	def prob_1(m):
		return sum(m.pi[s] for s in m.S) == 1
	m.prob1 = Constraint(rule=prob_1)
	def prob_2(m,s):
		return sum(m.pi[s]*m.Q[s,t] for t in m.S) <= m.M*(1-sum(m.z[h] for (h,s) in m.HS))
	m.prob2 = Constraint(rule=prob_2)
	def prob_3(m,s):
		return sum(m.pi[s]*m.Q[s,t] for t in m.S) >= m.M*(sum(m.z[h] for (h,s) in m.HS)-1)
	m.prob3 = Constraint(rule=prob_3)
	def subOBJ(m):
		return 
	m.obj = Objective(rule=subOBJ)


