import mdtraj as md
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

# CUDA INPUT
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'double', 'CudaDeviceIndex': '0,1', 'UseCpuPme': 'True'}

# SYSTEM INPUT
sysdir = 'SYS/'
sysgro = sysdir + 'system.gro'
systop = sysdir +'topol.top'

gro = GromacsGroFile(sysgro)
top = GromacsTopFile(systop, periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(hydrogenMass=4*amu,nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=AllBonds)

# SYSTEM OUTPUT
mddir = 'MD/'

# SIMULATION SETTING
temperature = 300.0*kelvin
fric_const = 1/picosecond
dt = 0.0050*picosecond
pressure = 1*bar
nvtstep =   1000
nptstep =  20000
mdstep  =  20000
totstep = nptstep + mdstep 
recstep =   1000

integrator = LangevinIntegrator(temperature, fric_const, dt)

# EM NVT SIMULATION START
simulation_nvt = Simulation(top.topology, system, integrator, platform, properties)
simulation_nvt.context.setPositions(gro.positions)

print('Minimizing...')
empdb = mddir + 'em.pdb'
simulation_nvt.minimizeEnergy()

print('Saving...')
positions = simulation_nvt.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation_nvt.topology, positions, open(empdb, 'w'))

print('NVT Equilibrating...')
nvtpdb = mddir + 'nvt.pdb'
nvtlog = mddir + 'nvt.log'
simulation_nvt.context.setVelocitiesToTemperature(temperature)
simulation_nvt.reporters.append(StateDataReporter(nvtlog, recstep, step=True,
        totalEnergy=True, temperature=True, density=True, 
        progress=True, remainingTime=True, speed=True, totalSteps=nvtstep, separator='\t'))

simulation_nvt.step(nvtstep)

print('Saving...')
positions = simulation_nvt.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation_nvt.topology, positions, open(nvtpdb, 'w'))

print('\nNPT Equilibrating...')
nptlog = mddir +'npt.log'
nptpdb = mddir + 'npt.pdb'
nptchk = mddir +'npt.chk'
system.addForce(MonteCarloBarostat(pressure, temperature))
integrator = LangevinIntegrator(temperature, fric_const, dt)
simulation_npt = Simulation(top.topology, system, integrator, platform, properties)
simulation_npt.context.setPositions(positions)

simulation_npt.reporters.append(PDBReporter(nptpdb, recstep))
simulation_npt.reporters.append(StateDataReporter(nptlog, recstep, time=True,
        totalEnergy=True, temperature=True, density=True,
        progress=True, remainingTime=True, speed=True, totalSteps=totstep, separator='\t'))

simulation_npt.step(nptstep)

print('\nRunning Production...')
mdpdb = mddir + 'md.pdb'
mdlog = mddir + 'md.log'
mdchk = mddir + 'md.chk'

simulation_npt.reporters.append(PDBReporter(mdpdb, recstep))
simulation_npt.reporters.append(StateDataReporter(mdlog, recstep, time=True,
        totalEnergy=True, temperature=True, density=True, 
        progress=True, remainingTime=True, speed=True, totalSteps=totstep, separator='\t'))
simulation_npt.reporters.append(CheckpointReporter(nptchk, recstep))

simulation_npt.step(mdstep)
print('Done!\n')

