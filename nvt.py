import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

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

# SIMULATION SETTING
temperature = 300.0*kelvin
fric_const = 1/picosecond
dt = 0.002*picosecond
pressure = 1*bar
nvtstep =   5000
nptstep =  50000
mdstep  = 500000
totstep = nptstep + mdstep 
recstep = 1000

integrator = LangevinIntegrator(temperature, fric_const, dt)

# EM NVT SIMULATION START
simulation_nvt = Simulation(top.topology, system, integrator)
simulation_nvt.context.setPositions(gro.positions)

print('Minimizing...')
simulation_nvt.minimizeEnergy()

print('NVT Equilibrating...')
nvtpdb = 'MD/nvt.pdb'
nvtlog = 'MD/nvt.log'
simulation_nvt.context.setVelocitiesToTemperature(temperature)
simulation_nvt.reporters.append(PDBReporter(nvtpdb, nvtstep))
simulation_nvt.reporters.append(StateDataReporter(nvtlog, 1000, step=True,
        totalEnergy=True, temperature=True, density=True, 
        progress=True, remainingTime=True, speed=True, totalSteps=nvtstep, separator='\t'))

simulation_nvt.step(nvtstep)

print('Done!\n')
