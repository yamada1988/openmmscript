import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

os.chdir('../')

# SYSTEM INPUT
sysdir = 'SYS/'
mddir = 'MD/'
sysgro = mddir + 'nvt.gro'
systop = sysdir +'topol.top'

gro = GromacsGroFile(sysgro)
top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=AllBonds)

# SIMULATION SETTING
temperature = 300.0*kelvin
fric_const = 1/picosecond
dt = 0.002*picosecond
pressure = 1*bar
nvtstep =   5000
nptstep =  50000
mdstep  =  50000
totstep = nptstep + mdstep 
rectime = 1000

# SYSTEM OUTPUT
nptlog = mddir +'npt.log'
nptpdb = mddir + 'npt.pdb'
nptgro = mddir +'npt.gro'

print('\nNPT Equilibrating...')
system.addForce(MonteCarloBarostat(pressure, temperature))
integrator = LangevinIntegrator(temperature, fric_const, dt)
simulation_npt = Simulation(top.topology, system, integrator)
simulation_npt.context.setPositions(gro.positions)

simulation_npt.reporters.append(PDBReporter(nptpdb, rectime))
simulation_npt.reporters.append(StateDataReporter(nptlog, 1000, time=True,
        totalEnergy=True, temperature=True, density=True,
        progress=True, remainingTime=True, speed=True, totalSteps=nptstep, separator='\t'))

simulation_npt.step(nptstep)

t_npt = md.load(nptpdb)
t_npt[-1].save_gro(nptgro)

print('Done!')
