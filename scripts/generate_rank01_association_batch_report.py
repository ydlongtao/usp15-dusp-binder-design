#!/usr/bin/env python3
"""Generate an English, self-contained report for the three rank01 association runs."""
from __future__ import annotations
import base64, json, html, csv
from pathlib import Path

ROOT = Path("local_results/rank01_association_dell")
OUT = ROOT / "USP15_rank01_association_results_en.html"
SEQ = "MKIKLVFSDGTEVEVEVDPSDTVLELKKKIEELTGYKPEQLLLFHKGKKLEDGKSLTYHGVKEGDTIHVNIVKEEE"

def data_uri(path: Path, mime: str) -> str:
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode()

summaries = {}
metrics = {}
for i in range(3):
    d = ROOT / f"seed{i}" / "analysis"
    summaries[f"seed{i}"] = json.loads((d / "association_summary.json").read_text())
    metrics[f"seed{i}"] = list(csv.DictReader((d / "association_metrics.csv").open()))
video = data_uri(ROOT / "seed0/analysis/animation/rank01_association.mp4", "video/mp4")
gif = data_uri(ROOT / "seed0/analysis/animation/rank01_association.gif", "image/gif")
payload = json.dumps({"summary": summaries, "metrics": metrics}, separators=(",", ":"))
rows = []
for seed, s in summaries.items():
    c = s["com_distance_a"]
    m = s["min_heavy_atom_distance_a"]
    occ = s["hotspot_occupancy"]
    rows.append(f"<tr><td>{seed}</td><td>{s['frames']}</td><td>{s['duration_ns']:.1f}</td><td>{c['min']:.2f}–{c['max']:.2f}</td><td>{c['median']:.2f}</td><td>{m['min']:.2f}</td><td>{', '.join(f'{k}: {100*v:.1f}%' for k,v in occ.items())}</td></tr>")

doc = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>USP15 rank01 Association MD Results</title>
<style>body{{font:15px system-ui,-apple-system,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#172033}}h1{{color:#123b68}}h2{{margin-top:2rem}}.note{{background:#fff7df;border-left:4px solid #d97706;padding:1rem;margin:1rem 0}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem}}.card{{border:1px solid #ccd6e2;border-radius:8px;padding:1rem;background:#fff}}video,img{{max-width:100%;border-radius:6px}}svg{{width:100%;height:290px;background:#f8fafc;border:1px solid #e5e7eb}}table{{border-collapse:collapse;width:100%;font-size:.9rem}}th,td{{border:1px solid #d7dee8;padding:.45rem;text-align:left;vertical-align:top}}th{{background:#f3f4f6}}pre{{white-space:pre-wrap;word-break:break-all;background:#f3f4f6;padding:1rem;border-radius:6px}}.legend span{{margin-right:1rem}}.orange{{color:#ea580c}}.slate{{color:#64748b}}.yellow{{color:#b58900}}</style></head><body>
<h1>USP15 DUSP–rank01 Association MD Results</h1>
<p><b>System:</b> USP15 DUSP and rank01 binder; three independent 100 ns NPT association attempts from an approximately 35 Å separated state. OpenMM 8.5.2, OpenCL V100, 300 K, 1 bar, 2 fs, no production distance or position restraints.</p>
<div class="note"><b>Conclusion boundary:</b> These trajectories did not show a global binder–target association event. The binder–target center-of-mass distance remained approximately 32–37 Å in all three runs. The animation is therefore a real unbiased association attempt and should not be described as proof of binding. Transient local atom contacts and hotspot occupancies are reported as diagnostics only.</div>
<h2>Binder sequence</h2><pre>{SEQ}</pre><p>76 aa; no cysteine. Chain A: binder (orange). Chain B: USP15 DUSP (slate blue). Hotspot residues are highlighted in yellow.</p>
<h2>Real trajectory animation (seed0 representative)</h2><video controls loop preload="metadata" src="{video}"></video><p><img alt="rank01 association animation" src="{gif}"></p>
<h2>Run summary</h2><table><tr><th>Run</th><th>Frames</th><th>Duration (ns)</th><th>COM distance range (Å)</th><th>COM median (Å)</th><th>Minimum heavy-atom distance (Å)</th><th>Hotspot occupancy</th></tr>{''.join(rows)}</table>
<div class="grid"><div class="card"><h2>COM distance over time</h2><svg id="dist" viewBox="0 0 900 290"></svg></div><div class="card"><h2>Hotspot contacts over time (seed0)</h2><svg id="hot" viewBox="0 0 900 290"></svg></div></div>
<h2>Interpretation</h2><ul><li>COM distance is the distance between binder and target centers of mass; it is the primary global association diagnostic.</li><li>Minimum heavy-atom distance detects local close approaches and does not by itself establish a bound complex.</li><li>Hotspot occupancy is the fraction of frames with a hotspot heavy-atom contact within 4.5 Å.</li><li>No KD, kon/koff, absolute binding free energy, or experimental affinity is inferred from these simulations.</li></ul>
<h2>Reproducibility</h2><p>Raw XTC trajectories, checkpoints, final PDB/state files, CSV time series, corrected analysis tables, and the original failed hotspot-mapping analysis directories are included alongside this report.</p>
<script id="data" type="application/json">{payload}</script><script>
const P=JSON.parse(document.getElementById('data').textContent);const C=['#2563eb','#ea580c','#16a34a'];
function line(id,series,key,labels){{const s=document.getElementById(id),W=900,H=290,L=55,R=15,T=20,B=35,all=series.flatMap(x=>x.data.map(p=>+p[key])),ymin=Math.min(...all),ymax=Math.max(...all),xmax=Math.max(...series.flatMap(x=>x.data.map(p=>+p.time_ns)));let z=`<line x1="${{L}}" y1="${{H-B}}" x2="${{W-R}}" y2="${{H-B}}" stroke="#64748b"/><line x1="${{L}}" y1="${{T}}" x2="${{L}}" y2="${{H-B}}" stroke="#64748b"/>`;series.forEach((q,j)=>{{const pts=q.data.map(p=>`${{L+(+p.time_ns)/(xmax)*(W-L-R)}},${{H-B-(+p[key]-ymin)/(ymax-ymin||1)*(H-T-B)}}`).join(' ');z+=`<polyline points="${{pts}}" fill="none" stroke="${{C[j]}}" stroke-width="1.5"/><text x="${{W-100}}" y="${{T+16*j}}" fill="${{C[j]}}">${{q.name}}</text>`}});z+=`<text x="8" y="${{T+10}}">${{ymin.toFixed(1)}}–${{ymax.toFixed(1)}}</text><text x="${{W-70}}" y="${{H-8}}">time (ns)</text>`;s.innerHTML=z}}
line('dist',Object.entries(P.metrics).map(([name,data])=>({{name,data}})),'binder_target_com_a');const h=P.metrics.seed0.filter((_,i)=>i%30===0);line('hot',[{{name:'seed0',data:h}}],'hotspots_contacting');
</script></body></html>'''
OUT.write_text(doc)
print(OUT)
