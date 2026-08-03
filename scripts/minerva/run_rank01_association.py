#!/usr/bin/env python3
"""Run an unbiased rank01 association trajectory from a separated state.

The prepared ff19SB/OPC solvated system is reused.  The binder coordinates
are translated away from the DUSP COM before minimization; no binder-target
distance or positional restraint is present during minimization, equilibration,
or production.  This is a kinetic association attempt, not a guaranteed
binding event or a free-energy calculation.
"""
from __future__ import annotations
import argparse, gc, json, time, math
from pathlib import Path
import openmm
from openmm import unit, XmlSerializer, Vec3
from openmm.app import PDBFile, Simulation, StateDataReporter, CheckpointReporter, XTCReporter

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--prepared-dir',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--seed',type=int,required=True)
    ap.add_argument('--production-ns',type=float,default=100.0)
    ap.add_argument('--separation-a',type=float,default=35.0)
    ap.add_argument('--device-index',default='0')
    ap.add_argument('--platform',default='CUDA',choices=('CUDA','OpenCL'))
    args=ap.parse_args(); out=args.output_dir; out.mkdir(parents=True,exist_ok=True)
    status_path=out/'status.json'
    if status_path.exists() and json.loads(status_path.read_text()).get('status')=='completed': return
    pdb=PDBFile(str(args.prepared_dir/'solvated.pdb'))
    protein_pdb=PDBFile(str(args.prepared_dir/'protein_protonated.pdb'))
    system=XmlSerializer.deserialize((args.prepared_dir/'system_production.xml').read_text())
    # Production must remain free of positional/distance restraints.
    bad=[f.getName() for f in system.getForces() if isinstance(f,openmm.CustomExternalForce)]
    if bad: raise RuntimeError(f'production restraints detected: {bad}')
    inpcrd=openmm.app.AmberInpcrdFile(str(args.prepared_dir/'solvated.inpcrd'))
    positions=inpcrd.positions
    # Solvation can collapse chain IDs in the full PDB.  Use the audited
    # protein-only topology to map the first protein atom indices reliably.
    binder=[a.index for a in protein_pdb.topology.atoms() if a.residue.chain.id=='A']
    target=[a.index for a in protein_pdb.topology.atoms() if a.residue.chain.id=='B']
    if not binder or not target: raise RuntimeError('could not identify chain A binder and chain B target')
    def mean(ids): return sum((positions[i] for i in ids), positions[ids[0]]*0)/len(ids)
    bc,tc=mean(binder),mean(target); vec=bc-tc
    # Convert all geometry to numeric nanometers before constructing the
    # displacement, avoiding mixed Å/nm Quantity arithmetic.
    vec_nm=vec.value_in_unit(unit.nanometer)
    vnm=[vec_nm.x,vec_nm.y,vec_nm.z]
    norm=math.sqrt(sum(x*x for x in vnm))
    if norm < 1e-8: vnm=[1.0,0.0,0.0]; norm=1.0
    direction=[x/norm for x in vnm]; desired=args.separation_a/10.0
    shift=Vec3(*(direction[i]*desired-vnm[i] for i in range(3)))*unit.nanometer
    pos=list(positions)
    for i in binder: pos[i]=pos[i]+shift
    platform=openmm.Platform.getPlatformByName(args.platform); props={'DeviceIndex':args.device_index,'Precision':'mixed'}
    if args.platform == 'CUDA': props['CudaCompiler']='nvcc'
    integrator=openmm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,2*unit.femtosecond); integrator.setRandomNumberSeed(args.seed)
    for f in system.getForces():
        if isinstance(f,openmm.MonteCarloBarostat): f.setRandomNumberSeed(args.seed+10000)
    sim=Simulation(pdb.topology,system,integrator,platform,props); sim.context.setPositions(pos)
    sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=20000)
    sim.context.setVelocitiesToTemperature(300*unit.kelvin,args.seed)
    (out/'separated_minimized.pdb').write_text('')
    with (out/'separated_minimized.pdb').open('w') as h: PDBFile.writeFile(pdb.topology,sim.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(),h,keepIds=True)
    # Unrestrained NPT pre-equilibration from the separated state.
    sim.step(500000)  # 1 ns at 2 fs; no distance/position restraints.
    state=sim.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
    (out/'equilibration.state.xml').write_text(XmlSerializer.serialize(state))
    with (out/'separated_equilibrated.pdb').open('w') as h:
        PDBFile.writeFile(pdb.topology,state.getPositions(),h,keepIds=True)
    sim.reporters.append(StateDataReporter(str(out/'association.csv'),50000,step=True,time=True,potentialEnergy=True,kineticEnergy=True,totalEnergy=True,temperature=True,volume=True,density=True,speed=True,separator=','))
    sim.reporters.append(XTCReporter(str(out/'association_protein.xtc'),5000,atomSubset=binder+target,enforcePeriodicBox=False))
    sim.reporters.append(CheckpointReporter(str(out/'association.chk'),500000))
    steps=int(round(args.production_ns*500000)); started=time.time(); sim.step(steps)
    final=sim.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
    (out/'association_final.state.xml').write_text(XmlSerializer.serialize(final))
    with (out/'association_final.pdb').open('w') as h: PDBFile.writeFile(pdb.topology,final.getPositions(),h,keepIds=True)
    status={'status':'completed','seed':args.seed,'production_ns':args.production_ns,'production_steps':steps,'initial_separation_a':args.separation_a,'production_restraints':False,'openmm_version':openmm.__version__,'platform':platform.getName(),'elapsed_s':time.time()-started,'binder_atoms':len(binder),'target_atoms':len(target)}
    status_path.write_text(json.dumps(status,indent=2)+'\n'); del sim; gc.collect(); print(json.dumps(status,indent=2))
if __name__=='__main__': main()
