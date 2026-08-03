#!/usr/bin/env python3
"""Analyze rank01 separated-state association trajectories and make movie frames."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

HOTSPOTS=(45,47,48,50,52,56)  # chain-relative USP15 DUSP residues

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--topology',type=Path,required=True); ap.add_argument('--trajectory',type=Path,required=True); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--movie-frames',type=int,default=121); args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    u=mda.Universe(str(args.topology),str(args.trajectory)); binder=u.select_atoms('segid A'); target=u.select_atoms('segid B')
    if len(binder)==0 or len(target)==0: binder=u.select_atoms('chainID A'); target=u.select_atoms('chainID B')
    if len(binder)==0 or len(target)==0: raise RuntimeError('Could not identify binder chain A and target chain B')
    target_res={r.resid:r for r in target.residues}
    # The prepared protein-only PDB renumbers chain B consecutively after
    # chain A. Map the DUSP-local hotspot numbers onto that chain explicitly.
    target_ids=sorted(target_res)
    offset=(target_ids[0]-1) if target_ids else 0
    hs=[(x, target_res.get(offset+x)) for x in HOTSPOTS]
    rows=[]; n=len(u.trajectory); stride=max(1,n//max(1,args.movie_frames-1)); frame_dir=args.out/'movie_frames'; frame_dir.mkdir(exist_ok=True)
    for j,ts in enumerate(u.trajectory):
        b=binder.positions; t=target.positions; d=distance_array(b,t); com=float(np.linalg.norm(b.mean(0)-t.mean(0))); min_d=float(d.min()); ca=len(distance_array(binder.select_atoms('name CA'),target.select_atoms('name CA')) <= 8.0)
        hs_counts=[]
        for _, r in hs:
            hs_counts.append(int(r is not None and (distance_array(b,r.atoms) <= 4.5).any()))
        rows.append({'frame':j,'time_ns':float(ts.time)/1000.0,'binder_target_com_a':com,'min_heavy_atom_distance_a':min_d,'ca_contacts_le_8a':ca,'hotspots_contacting':sum(hs_counts),**{f'hotspot_{x}_contact':v for x,v in zip(HOTSPOTS,hs_counts)}})
        if j % stride==0 or j==n-1:
            u.atoms.write(str(frame_dir/f'frame_{len(list(frame_dir.glob("frame_*.pdb"))):04d}.pdb'))
    with (args.out/'association_metrics.csv').open('w',newline='') as h: w=csv.DictWriter(h,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    arr={k:np.array([float(r[k]) for r in rows]) for k in rows[0] if k not in ('frame','time_ns')}
    summary={'frames':n,'duration_ns':rows[-1]['time_ns'],'com_distance_a':{'min':float(arr['binder_target_com_a'].min()),'median':float(np.median(arr['binder_target_com_a'])),'max':float(arr['binder_target_com_a'].max())},'min_heavy_atom_distance_a':{'min':float(arr['min_heavy_atom_distance_a'].min()),'median':float(np.median(arr['min_heavy_atom_distance_a']))},'hotspot_occupancy':{str(x):float(arr[f'hotspot_{x}_contact'].mean()) for x in HOTSPOTS},'note':'Association attempt from a separated state; no binding event is assumed if the trajectory remains separated.'}
    (args.out/'association_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
