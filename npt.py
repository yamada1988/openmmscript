import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

# SYSTEM INPUT
sysdir = 'SYS/'
mddir = 'MD/'
sysgro = mddir + 'nvt.gro'
systop = sysdir +'topol.top'

gro = GromacsGroFile(sysgro)
top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(hydrogenMass=4*amu, nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
         constraints=AllBonds)

# SIMULATION SETTING
temperature = 300.0*kelvin
fric_const = 1/picosecond
dt = 0.0050*picosecond
pressure = 1*bar
nvtstep =     5000
nptstep =    10000
mdstep  =  4000000
totstep = nptstep + mdstep 
rectime =    10000

# SYSTEM OUTPUT
nptlog = mddir +'npt.log'
nptpdb = mddir + 'npt.pdb'
nptgro = mddir +'npt.gro'
nptchk = mddir +'npt.chk'

print('\nNPT Equilibrating...')
system.addForce(MonteCarloBarostat(pressure, temperature))
integrator = LangevinIntegrator(temperature, fric_const, dt)
simulation_npt = Simulation(top.topology, system, integrator)
simulation_npt.context.setPositions(gro.positions)

simulation_npt.reporters.append(PDBReporter(nptpdb, rectime))
simulation_npt.reporters.append(StateDataReporter(nptlog, rectime, time=True,
        totalEnergy=True, temperature=True, density=True,
        progress=True, remainingTime=True, speed=True, totalSteps=totstep, separator='\t'))
simulation_npt.reporters.append(CheckpointReporter(nptchk, rectime))

simulation_npt.step(nptstep)

print('\nRunning Production...')
mdpdb = mddir + 'md.pdb'
mdlog = mddir + 'md.log'
mdchk = mddir + 'md.chk'

simulation_npt.reporters.append(PDBReporter(mdpdb, rectime))
simulation_npt.reporters.append(StateDataReporter(mdlog, rectime, time=True,
        totalEnergy=True, temperature=True, density=True, 
        progress=True, remainingTime=True, speed=True, totalSteps=totstep, separator='\t'))
simulation_npt.reporters.append(CheckpointReporter(nptchk, rectime))

simulation_npt.step(mdstep)
print('Done!\n')


