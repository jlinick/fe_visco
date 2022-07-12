#!/usr/bin/env python3

import os
import math
import sys
from fenics import *
from ufl import nabla_div
import numpy as np
import matplotlib.pyplot as plt

print('------------initializing simulation-------------')
# set parameters
TOTAL_DAYS = 10.            # total length of simulation (days)
SEC_PER_DAY = 86400
T = TOTAL_DAYS*SEC_PER_DAY  # final time (sec)
num_steps = 100             # number of time steps
dt = T / num_steps          # time step size
#print('Running over {} days, with {} increments'.format(TOTAL_DAYS,num_steps))

#LAME_MU = Constant(0.8354E9)     # shear modulus (G or mu)
E = YOUNGS_MODULUS = 9.1E9 # E = 9.1 GPa
K = 8.9E9 # 8.9 GPa empirical (Schulson and Duval, 2009)
LAMBDA = 3*K*(3*K-E)/(9*K-E)
G = SHEAR_MODULUS = MU = 3*K*E/(9*K-E)
v = (3*K-E)/(6*K)
print('Moduli for Glacial Ice')
print('Bulk Modulus:   {:2.2f} GPa'.format(K/1E9))
print('Lambda:         {:2.2f} GPa'.format(LAMBDA/1E9))
print('Shear Modulus:  {:2.2f} GPa'.format(G/1E9))
print('Poisson Ratio:  {:2.2f}'.format(v))
print()

ICE_DENSITY = Constant(917.)      # mass density of ice kg/m3
WATER_DENSITY = Constant(997.)    # mass density of water kg/m3
GRAVITY = Constant(9.8)
# viscoelastic parameters
ETA_M = 8.5E5 # empirical
ETA_K = Constant(0.0)

# Create mesh and define function space
LENGTH = 1000.
WIDTH = 1000.
HEIGHT = 100.
LAKE_RADIUS = 100.
MAX_LAKE_DEPTH = 5.

# times for filling/draining lake
T_START_FILL = T*.1#T * .2
T_START_DRAIN = T*.6#T_START_FILL + 5. * SEC_PER_DAY # +5 days
T_END_DRAIN = T*.7#T_START_DRAIN + 1.* SEC_PER_DAY   # +1 day

# number of points on mesh
N_POINTS_LENGTH = 40
N_POINTS_WIDTH = 40
N_POINTS_HEIGHT = 5

# Generalized-alpha method parameters
alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
gamma   = Constant(0.5+alpha_f-alpha_m)
beta    = Constant((gamma+0.5)**2/4.)

ROLLER_BOUNDARIES = False # set to True/False

output_products_dir = 'products'
if not os.path.exists(output_products_dir):
    os.makedirs(output_products_dir)
output_prefix = os.path.join(output_products_dir, 'result')

# END PARAMETERS


# Mesh and Vector Function Space
mesh = BoxMesh(
    Point(0.0, 0.0, 0.0),
    Point(LENGTH, WIDTH, HEIGHT),
    N_POINTS_LENGTH,
    N_POINTS_WIDTH,
    N_POINTS_HEIGHT,
)
# function space for displacement, velocity, and acceleration
V = VectorFunctionSpace(
    mesh,
    "Lagrange",
    1,
)
# function space for stresses
Vsig = TensorFunctionSpace(mesh, "DG", 0)

# Test and trial functions
du = TrialFunction(V)
u_ = TestFunction(V)
# Current (unknown) displacement
u = Function(V, name="Displacement")
# Fields from previous time step (displacement, velocity, acceleration)
u_old = Function(V)
v_old = Function(V)
a_old = Function(V)

# set boundary conditions
tol = 1E-14
def left(x, on_boundary):
    return on_boundary and abs(x[0]) < tol

def right(x, on_boundary):
    return on_boundary and abs(x[0] - LENGTH) < tol

def back(x, on_boundary):
    return on_boundary and abs(x[1]) < tol

def front(x, on_boundary):
    return on_boundary and abs(x[1] - WIDTH) < tol

def bottom(x, on_boundary):
	return on_boundary and abs(x[2] < tol)

def lake(x, on_boundary):
	return on_boundary and abs(x[2] - HEIGHT) < tol and (x[0] - LENGTH/2) ** 2 + (x[1] - WIDTH/2) ** 2 < LAKE_RADIUS ** 2

if ROLLER_BOUNDARIES:
    bcl = DirichletBC(V.sub(0), Constant((0)), left)
    bcr = DirichletBC(V.sub(0), Constant((0)), right) 
    bcb = DirichletBC(V.sub(0), Constant((0)), back) 
    bcf = DirichletBC(V.sub(0), Constant((0)), front) 
else:
    zero = Constant((0.,0.,0.))
    bcl = DirichletBC(V, zero, left)
    bcr = DirichletBC(V, zero, right) 
    bcb = DirichletBC(V, zero, back) 
    bcf = DirichletBC(V, zero, front) 
bcs = [bcl, bcr, bcb, bcf]

# Create mesh function over the cell facets
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)

buoyancy_boundary = AutoSubDomain(bottom)
buoyancy_boundary.mark(boundary_subdomains, 3)

lake_boundary = AutoSubDomain(lake)
lake_boundary.mark(boundary_subdomains, 4)

# Define measure for boundary condition integral
dss = ds(subdomain_data=boundary_subdomains)

# recursive building of lake drainage equation
eq1 = ('t < T_START_FILL' , '0.0') # before lake fill
eq2 = ('t <= T_START_DRAIN', 'GRAVITY * WATER_DENSITY * MAX_LAKE_DEPTH / (T_START_DRAIN - T_START_FILL) * (t - T_START_FILL)') # lake filling
eq3 = ('t <= T_END_DRAIN', 'GRAVITY * WATER_DENSITY * (-MAX_LAKE_DEPTH / (T_END_DRAIN - T_START_DRAIN) * (t - T_START_DRAIN) + MAX_LAKE_DEPTH)')  # lake drainage
eq4 = ('t > T_END_DRAIN', '0.0') # no lake

def build_expression(explist):
    '''builds the lake drainage expression'''
    if len(explist) <= 1 :
        return '{}'.format(explist[0][1])
    return '{} ? {} : ({})'.format(explist[0][0], explist[0][1], build_expression(explist[1:]))

lake_eq = build_expression([eq1,eq2,eq3,eq4])
fl = Expression( ('0.','0.', lake_eq), GRAVITY=GRAVITY, MAX_LAKE_DEPTH=MAX_LAKE_DEPTH, WATER_DENSITY=WATER_DENSITY, T_START_FILL=T_START_FILL, T_START_DRAIN=T_START_DRAIN, T_END_DRAIN=T_END_DRAIN, t=0., degree=1)


# expression for buoyancy
#fb_expression = Expression( ('0.','0.', 'displacement >= 0. ? 0. : GRAVITY * WATER_DENSITY * displacement'), displacement=0., GRAVITY=GRAVITY, WATER_DENSITY=WATER_DENSITY, degree=1)
class Buoyancy(UserExpression):
    def eval(self, value, x):
        if value[2] < 0.:
            return np.array((0., 0., -GRAVITY * WATER_DENSITY * value[2]))
        return np.array([0.,0.,0.])
        #return fb_expression(displacement=value[2])
    def value_shape(self):
        return (3,)
fb = Buoyancy()
#fb = Constant((0.0, 0.0, -WATER_DENSITY * GRAVITY))


class Gravity(UserExpression):
    def eval(self, value, x):
        return np.array((0., 0., GRAVITY * WATER_DENSITY))
    def value_shape(self):
        return (3,)
fg = Gravity()
#eqb = 'u < 0. ? GRAVITY * DENSITY * -u : 0.0'
#fb = Expression(('0.0', '0.0', eqb), DENSITY=DENSITY, GRAVITY=GRAVITY, u=u_.vector(), degree=0)
#fb = Constant((0.0, 0.0, WATER_DENSITY * GRAVITY))
grav_forcing = Constant((0.0, 0.0, ICE_DENSITY * GRAVITY))

# Define strain and stress
def epsilon(u):
    engineering_strain = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
    return engineering_strain

# stress tensor
def sigma(u):
    cauchy_stress = (LAMBDA * tr(epsilon(u)) * Identity(3) + 2 * MU * epsilon(u))
    return cauchy_stress

# stress tensor using nabla
def st(u):
    return LAMBDA * nabla_div(u) * Identity(3) + 2 * MU * epsilon(u)

def deviatoric(u):
    return sigma(u) - (1./3)*tr(sigma(u))*Identity(3)  # deviatoric stress

# Mass form
def m(u, u_):
    return ICE_DENSITY*inner(u, u_)*dx

# Elastic stiffness form
def k(u, u_):
    return inner(sigma(u), sym(grad(u_)))*dx

# viscous damping form
def c(u, u_):
    return ETA_M*m(u, u_) + ETA_K*k(u, u_)

# Work of gravity
def Wgrav(u_):
    return dot(u_, grav_forcing)*dx()

# Work of buoyancy
def Wbuoy(u_):
    return dot(u_, fb)*dss(3) # buoyancy only applies to the bottom (dss)

# Work of lake loading
def Wlake(u_):
	return dot(u_, fl)*dss(4)


# the system variational form is generated by expressing the new acceleration and velocity as a function of du using update_a and update_v.
# Intermediate averages are generated from avg. The weak form evolution equation is then written using all these quantities.

# Update formula for acceleration
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
def update_a(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)
    return (u-u_old-dt_*v_old)/beta_/dt_**2 - (1-2*beta_)/2/beta_*a_old

# Update formula for velocity
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
def update_v(a, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)
    return v_old + dt_*((1-gamma_)*a_old + gamma_*a)

def update_fields(u, u_old, v_old, a_old):
    """Update fields at the end of each time step.""" 

    # Get vectors (references)
    u_vec, u0_vec  = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector() 

    # use update functions using vector arguments
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # Update (u_old <- u)
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()

def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new

# Residual
a_new = update_a(du, u_old, v_old, a_old, ufl=True)
v_new = update_v(a_new, u_old, v_old, a_old, ufl=True)

res = m(avg(a_old, a_new, alpha_m), u_) + c(avg(v_old, v_new, alpha_f), u_) + k(avg(u_old, du, alpha_f), u_) + Wlake(u_) + Wgrav(u_) + Wbuoy(u_)

a_form = lhs(res)
L_form = rhs(res)

# Define solver for reusing factorization
K, res = assemble_system(a_form, L_form, bcs)
solver = LUSolver(K, "mumps")
solver.parameters["symmetric"] = True

# starting the time stepping loop
time = np.linspace(0, T, num_steps+1)

middle_displacement = np.zeros((num_steps+1,)) # at midpoint
lake_depth = np.zeros((num_steps+1,)) # at midpoint
midpoint_von_mises = np.zeros((num_steps+1,))
energies = np.zeros((num_steps+1, 7))

E_damp = 0
E_ext = 0
E_grav = 0
E_lake = 0
E_buoy = 0
sig = Function(Vsig, name="sigma")

data_file = f"{output_prefix}_data.xdmf"
xdmf_file = XDMFFile(data_file)
xdmf_file.parameters["flush_output"] = True
xdmf_file.parameters["functions_share_mesh"] = True
xdmf_file.parameters["rewrite_function_mesh"] = False

# loading is first evaluated at :math:`t=t_{n+1-\alpha_f}`. The right-hand side is assembled and the system is solved, updating quantities and post-processing stresses are computed
def local_project(v, V, u=None):
    """Element-wise projection using LocalSolver"""
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

for (i, dt) in enumerate(np.diff(time)):

    t = time[i+1]


    # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
    fl.t = t-float(alpha_f*dt)

    # Solve for new displacement
    res = assemble(L_form)
    for bc in bcs:
        bc.apply(res)
    solver.solve(K, u.vector(), res)


    # Update old fields with new quantities
    update_fields(u, u_old, v_old, a_old)

    # Save solution to XDMF format
    xdmf_file.write(u, t)

    # Compute stresses and save to file
    local_project(sigma(u), Vsig, sig)
    xdmf_file.write(sig, t)

    # update t for lake forcing
    fl.t = t
    # Record displacement and compute energies
    # Note: Only works in serial
    if MPI.comm_world.size == 1:
        # get displacement
        md = u(LENGTH/2., WIDTH/2., 0.)[2]
        middle_displacement[i+1] = md
        # get von-mises stress
        s  = sig(LENGTH/2., WIDTH/2., 0.) # sigma at midpoint (bottom of ice)
        vm = sqrt(3./2*np.inner(s, s)) # von mises stress
        midpoint_von_mises[i+1] = vm/1E6
        temp_ld = float(fl(None)[2]) / float(WATER_DENSITY * GRAVITY) # lake depth placeholder
        lake_depth[i+1] = temp_ld

        print('--------')
        print('time: {:5.1f} / {:5.1f} days'.format(t/SEC_PER_DAY,T/SEC_PER_DAY))
        print('lake depth:             {:10.2f} m'.format(temp_ld))
        print('von mises stress        {:10.2f} MPa'.format(vm/1E6))
        print('midpoint displacement:  {:10.2f} m'.format(md))
        
    E_elas = assemble(0.5*k(u_old, u_old))
    E_kin = assemble(0.5*m(v_old, v_old))
    E_damp += dt*assemble(c(v_old, v_old))
    E_grav += assemble(Wgrav(u-u_old))
    E_buoy += assemble(Wbuoy(u-u_old))
    E_lake += assemble(Wlake(u-u_old))
    E_tot = E_elas+E_kin+E_damp-E_grav-E_buoy-E_lake
    energies[i+1, :] = np.array([E_elas, E_kin, E_damp, E_grav, E_buoy, E_lake, E_tot])/1E9

#print(f'middle displacement: {middle_displacement}')
#print(f'lake depth: {lake_depth}')

print(f'results saved to : {data_file}')

time = time/SEC_PER_DAY # correct time to be in days

output_file = f'{output_prefix}_midpoint_displacement.png'
print(f'saving {output_file}')
plt.plot(time[40:], middle_displacement[40:])
plt.title('Midpoint Displacement')
plt.xlabel('time (days)')
plt.ylabel('displacement (m)')
plt.savefig(output_file)
plt.close()

output_file = f'{output_prefix}_lake_depth.png'
print(f'saving {output_file}')
plt.plot(time[40:], lake_depth[40:])
plt.title('Lake Depth')
plt.xlabel('time (days)')
plt.ylabel('lake depth (m)')
plt.savefig(output_file)
plt.close()

output_file = f'{output_prefix}_midpoint_stress.png'
print(f'saving {output_file}')
plt.plot(time[40:], midpoint_von_mises[40:])
plt.title('Effective Stress at Midpoint of Lake')
plt.xlabel('time (days)')
plt.ylabel('effective stress (MPa)')
plt.savefig(output_file)
plt.close()


output_file = f'{output_prefix}_energies.png'
print(f'saving {output_file}')
plt.plot(time, energies)
plt.title('Energies during Simulation')
plt.xlabel('time (days)')
plt.ylabel('Energy (GJ)')
plt.legend(['elastic', 'kinetic', 'damping', 'gravitational', 'buoyancy', 'lake gravitational', 'total'])
plt.savefig(output_file)
plt.close()
