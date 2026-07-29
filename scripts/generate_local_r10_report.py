#!/usr/bin/env python3
"""Build a self-contained preliminary HTML report from downloaded R10 analyses."""
import csv, json, html
from pathlib import Path

root = Path(__file__).resolve().parents[1] / "local_results" / "r10_current"
out = root / "USP15_R10_current_results.html"
runs = []
for p in sorted((root / "analysis").glob("rank*/seed*")):
    summary = p / "summary.json"
    frame = p / "per_frame.csv"
    sasa_frame = p / "buried_sasa_timeseries.csv"
    if not summary.exists() or not frame.exists():
        continue
    s = json.loads(summary.read_text())
    rows = list(csv.DictReader(frame.open()))
    # Keep the browser payload small while retaining the time-course shape.
    step = max(1, len(rows) // 240)
    series = [{k: float(r[k]) for k in ("time_ns", "binder_rmsd_a", "target_rmsd_a", "buried_sasa_a2", "native_contact_fraction") if k in r}
              for r in rows[::step]]
    sasa_series = []
    if sasa_frame.exists():
        sasa_rows = list(csv.DictReader(sasa_frame.open()))
        sasa_step = max(1, len(sasa_rows) // 240)
        sasa_series = [{"time_ns": float(r["time_ns"]), "buried_sasa_a2": float(r["buried_sasa_a2"])}
                       for r in sasa_rows[::sasa_step]]
    mmp = p / "mmpbsa" / "FINAL_RESULTS_MMPBSA.dat"
    runs.append({"name": p.relative_to(root / "analysis").as_posix(), "summary": s, "series": series, "sasa_series": sasa_series,
                 "mmpbsa": mmp.read_text(errors="replace")[:12000] if mmp.exists() else "not available"})

payload = json.dumps(runs, ensure_ascii=False).replace("</", "<\\/")
page = r'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>USP15 R10 当前 OpenMM MD 结果</title>
<style>body{font-family:system-ui,-apple-system,sans-serif;max-width:1280px;margin:2rem auto;padding:0 1rem;color:#1f2937}h1{margin-bottom:.2rem}.note{background:#fff7ed;border-left:4px solid #f97316;padding:1rem;margin:1rem 0}table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:.9rem}th,td{border:1px solid #d1d5db;padding:.45rem;text-align:right}th:first-child,td:first-child{text-align:left}th{background:#f3f4f6}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(480px,1fr));gap:1rem}.plot{border:1px solid #d1d5db;border-radius:6px;padding:.5rem}.plot h3{margin:.2rem}.plot svg{width:100%;height:260px;background:#fafafa}.legend{font-size:.8rem}.small{color:#6b7280;font-size:.9rem}pre{white-space:pre-wrap;background:#111827;color:#e5e7eb;padding:1rem;overflow:auto;font-size:.75rem}</style></head>
<body><h1>USP15 DUSP R10：当前 OpenMM MD 结果</h1><p class="small">本地快照：已完成重复的结构分析与相对 MM/GBSA。数据来自 Minerva 正式数组 258390007。</p>
<div class="note"><b>解释限制：</b>当前只有部分重复完成。本页用于查看轨迹诊断，不代表最终 MD 稳定性、绝对结合自由能、KD 或实验结合结论。MM/GBSA 仅作为相对比较指标。</div>
<h2>已完成重复</h2><table id="summary"><thead><tr><th>重复</th><th>轨迹 ns</th><th>帧数</th><th>Target RMSD 中位数 (Å)</th><th>Binder RMSD 中位数 (Å)</th><th>Buried SASA 中位数 (Å²)</th><th>热点占有率 ≥0.5</th><th>原生接触中位数</th></tr></thead><tbody></tbody></table>
<h2>时间序列图</h2><div class="grid"><div class="plot"><h3>Target / binder RMSD</h3><svg id="rmsd" viewBox="0 0 720 260"></svg><div class="legend">蓝：target RMSD；橙：binder RMSD</div></div><div class="plot"><h3>Buried SASA</h3><svg id="sasa" viewBox="0 0 720 260"></svg><div class="legend">每条曲线对应一个重复</div></div><div class="plot"><h3>Native contact fraction</h3><svg id="contacts" viewBox="0 0 720 260"></svg><div class="legend">每条曲线对应一个重复</div></div></div>
<h2>相对 MM/GBSA 输出</h2><div id="mmp"></div>
<script id="payload" type="application/json">''' + payload + r'''</script>
<script>
const runs=JSON.parse(document.getElementById('payload').textContent), colors=['#2563eb','#ea580c','#16a34a','#9333ea','#0891b2'];
const tbody=document.querySelector('#summary tbody');
for(const x of runs){const s=x.summary; const tr=document.createElement('tr'); tr.innerHTML=`<td>${x.name}</td><td>${s.duration_ns_observed}</td><td>${s.frames}</td><td>${s.target_rmsd_a.median.toFixed(2)}</td><td>${s.binder_rmsd_a.median.toFixed(2)}</td><td>${s.buried_sasa_a2.median.toFixed(1)}</td><td>${s.hotspots_with_occupancy_ge_0p5}</td><td>${s.native_contact_fraction.median.toFixed(3)}</td>`; tbody.appendChild(tr);}
function chart(id,key,ylabel){const svg=document.getElementById(id),W=720,H=260,L=52,R=12,T=14,B=32;const getSeries=x=>(key==='buried_sasa_a2'?x.sasa_series:x.series);let all=runs.flatMap(x=>getSeries(x).map(p=>p[key])).filter(Number.isFinite), ymax=Math.max(...all),ymin=Math.min(...all);if(!all.length){svg.innerHTML='<text x="20" y="40" fill="#b45309">逐帧 SASA 数据尚未生成</text>';return;}if(ymax===ymin)ymax=ymin+1; const x0=runs.flatMap(x=>getSeries(x).map(p=>p.time_ns));const xmax=Math.max(...x0),xmin=Math.min(...x0);let z=`<line x1="${L}" y1="${H-B}" x2="${W-R}" y2="${H-B}" stroke="#6b7280"/><line x1="${L}" y1="${T}" x2="${L}" y2="${H-B}" stroke="#6b7280"/><text x="5" y="${T+10}" font-size="11">${ylabel}</text><text x="${W-70}" y="${H-5}" font-size="11">time (ns)</text>`;runs.forEach((x,i)=>{let pts=getSeries(x).map(p=>`${L+(p.time_ns-xmin)/(xmax-xmin)*(W-L-R)},${H-B-(p[key]-ymin)/(ymax-ymin)*(H-T-B)}`).join(' ');z+=`<polyline points="${pts}" fill="none" stroke="${colors[i]}" stroke-width="1.5"/><text x="${W-120}" y="${T+14+i*14}" fill="${colors[i]}" font-size="11">${x.name}</text>`});svg.innerHTML=z;}
chart('rmsd','target_rmsd_a','Å'); chart('sasa','buried_sasa_a2','Å²'); chart('contacts','native_contact_fraction','fraction');
const m=document.getElementById('mmp'); for(const x of runs){m.innerHTML+=`<h3>${x.name}</h3><pre>${x.mmpbsa.replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]))}</pre>`;}
</script></body></html>'''
out.write_text(page, encoding="utf-8")
print(out)
