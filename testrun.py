import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

os.chdir('../')

# SYSTEM INPUT
sysdir = 'SYS/'
sysgro = sysdir + 'system.gro'
systop = sysdir +'topol.top'

gro = GromacsGroFile(sysgro)
top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=AllBonds)

# SYSTEM OUTPUT
mddir = 'MD'
emdir = mddir + 'em'
nvtdir = mddir + 'nvt'
#os.mkdir(mddir)


# SIMULATION SETTING
temperature = 300.0*kelvin
fric_const = 1/picosecond
dt = 0.002*picosecond
pressure = 1*bar
totstep = 50000
nvtstep =  5000
nptstep = 50000

integrator = LangevinIntegrator(temperature, fric_const, dt)

# EM NVT SIMULATION START
simulation_nvt = Simulation(top.topology, system, integrator)
simulation_nvt.context.setPositions(gro.positions)

print('Minimizing...')
simulation_nvt.minimizeEnergy()

print('NVT Equilibrating...')
nvtpdb = 'nvt.pdb'
simulation_nvt.context.setVelocitiesToTemperature(temperature)
simulation.reporters.append(PDBReporter(nvtpdb, nvtstep))
simulation_nvt.step(nvtstep)

# NPT SIMULATION SETTING
nvt = md.load_pdb(nvtpdb)

pdb = GromacsGroFile(nvtpdb)
top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=AllBonds)

print('NPT Equilibrating...')
system.addForce(MonteCarloBarostat(pressure, temperature))
integrator = LangevinIntegrator(temperature, fric_const, dt)
simulation_npt = Simulation(top.topology, system, integrator)
simulation_npt.
simulation.step(5000)

print('Running Production...')
integrator = LangevinIntegrator(temperature, fric_const, dt)
simulation = Simulation(top.topology, system, integrator)

simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        totalEnergy=True, temperature=True, density=True, 
        progress=True, remainingTime=True, speed=True, totalSteps=totstep, separator='\t'))

simulation.step(totstep)
print('Done!')
