load protein_protonated.pdb, complex
load association_protein.xtc, complex
hide everything, complex
show cartoon, complex
color orange, complex and chain A
color slate, complex and chain B
select hotspots, complex and chain B and resi 45+47+48+50+52+56
show sticks, hotspots
color yellow, hotspots
set cartoon_sampling, 14
set ray_trace_frames, 0
set movie_fps, 12
orient complex
mview store, 1
mview store, 121
set movie_panel, 1
mpng association_frames, width=1280, height=720, ray=0
quit
