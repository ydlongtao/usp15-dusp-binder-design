#!/usr/bin/env python3
"""Build a self-contained HTML report for a rank01 association trajectory.

The report deliberately labels the result as an unbiased association attempt;
it does not infer a binding event unless the supplied trajectory metrics show
close approach and persistent contacts.
"""
from __future__ import annotations
import argparse, csv, json, html
from pathlib import Path

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--video", type=str, default="")
    ap.add_argument("--gif", type=str, default="")
    ap.add_argument("--sequence", default="MKIKLVFSDGTEVEVEVDPSDTVLELKKKIEELTGYKPEQLLLFHKGKKLEDGKSLTYHGVKEGDTIHVNIVKEEE")
    args = ap.parse_args()
    rows = list(csv.DictReader(args.metrics.open()))
    summary = json.loads(args.summary.read_text())
    data = json.dumps(rows, separators=(",", ":"))
    seq = html.escape(args.sequence)
    video = html.escape(args.video)
    gif = html.escape(args.gif)
    duration = summary.get("duration_ns", "n/a")
    html_text = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>USP15 rank01 真实结合动力学</title>
<style>body{{font:15px system-ui,sans-serif;max-width:1180px;margin:2rem auto;padding:0 1rem;color:#172033}}h1{{color:#123b68}}.note{{background:#fff7df;border-left:4px solid #d97706;padding:1rem}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem}}.card{{border:1px solid #ccd6e2;border-radius:8px;padding:1rem;background:#fff}}svg{{width:100%;height:270px;background:#f8fafc}}video,img{{max-width:100%;border-radius:6px}}code,pre{{background:#f3f4f6;padding:.3rem;word-break:break-all}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #d7dee8;padding:.45rem;text-align:left}}@media(max-width:800px){{.grid{{grid-template-columns:1fr}}}}</style></head><body>
<h1>USP15 DUSP–rank01 真实关联动力学报告</h1>
<div class="note"><b>解释边界：</b>本报告来自从约 35 Å 分离态开始的无位置/距离约束 NPT 轨迹。它展示真实热运动和接近过程；只有在轨迹指标显示持续接触时，才能称为观察到结合。该模拟不能替代 SPR/MST 实验。</div>
<h2>Binder 序列</h2><pre>{seq}</pre>
<h2>三维关联动画</h2><p>橙色：binder chain A；蓝灰色：USP15 DUSP chain B；黄色：热点残基。</p>
{('<video controls loop preload="metadata" src="'+video+'"></video>' if video else '<p>动画尚未挂载。</p>')}
{('<p><img alt="association animation" src="'+gif+'"></p>' if gif else '')}
<h2>轨迹指标</h2><div class="grid"><div class="card"><h3>Binder–target COM 距离 (Å)</h3><svg id="dist" viewBox="0 0 720 270"></svg></div><div class="card"><h3>热点接触数</h3><svg id="hot" viewBox="0 0 720 270"></svg></div><div class="card"><h3>最小重原子距离 (Å)</h3><svg id="mind" viewBox="0 0 720 270"></svg></div><div class="card"><h3>Cα 接触数 (≤8 Å)</h3><svg id="ca" viewBox="0 0 720 270"></svg></div></div>
<h2>摘要</h2><table><tr><th>项目</th><th>值</th></tr><tr><td>轨迹时长</td><td>{html.escape(str(duration))} ns</td></tr><tr><td>COM 距离范围</td><td>{summary.get('com_distance_a',{})}</td></tr><tr><td>最小重原子距离</td><td>{summary.get('min_heavy_atom_distance_a',{})}</td></tr><tr><td>热点占有率</td><td><pre>{html.escape(json.dumps(summary.get('hotspot_occupancy',{}),ensure_ascii=False,indent=2))}</pre></td></tr></table>
<script>const D={data};function plot(id,key,color){{const s=document.getElementById(id),w=720,h=270,p=35,x=D.map(r=>+r.time_ns),y=D.map(r=>+r[key]),xmin=x[0],xmax=x[x.length-1]||1,ymin=Math.min(...y),ymax=Math.max(...y);if(ymax===ymin)ymax=ymin+1;const X=v=>p+(v-xmin)/(xmax-xmin||1)*(w-2*p),Y=v=>h-p-(v-ymin)/(ymax-ymin)*(h-2*p);let path=y.map((v,i)=>(i?'L':'M')+X(x[i]).toFixed(1)+','+Y(v).toFixed(1)).join('');s.innerHTML='<path d="'+path+'" fill="none" stroke="'+color+'" stroke-width="2"/><text x="'+p+'" y="18">'+ymin.toFixed(2)+'–'+ymax.toFixed(2)+'</text><text x="'+(w-p-40)+'" y="'+(h-8)+'">ns</text>'}}plot('dist','binder_target_com_a','#ea580c');plot('hot','hotspots_contacting','#eab308');plot('mind','min_heavy_atom_distance_a','#2563eb');plot('ca','ca_contacts_le_8a','#64748b');</script></body></html>'''
    args.out.write_text(html_text)

if __name__ == "__main__":
    main()
