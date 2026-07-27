#!/usr/bin/env python3
"""Add R10 complex figures, prospective MD parameters, and assay guidance."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PROPERTIES = {
    1: {"mw_da": 8653.81, "pI": 5.193, "net_charge_pH7_4": -5.623},
    2: {"mw_da": 8563.77, "pI": 4.904, "net_charge_pH7_4": -4.733},
    3: {"mw_da": 8592.80, "pI": 4.730, "net_charge_pH7_4": -6.688},
    4: {"mw_da": 8441.53, "pI": 4.910, "net_charge_pH7_4": -4.693},
    5: {"mw_da": 8473.76, "pI": 4.799, "net_charge_pH7_4": -4.731},
    6: {"mw_da": 8638.97, "pI": 5.888, "net_charge_pH7_4": -1.705},
    7: {"mw_da": 8321.37, "pI": 4.591, "net_charge_pH7_4": -7.691},
    8: {"mw_da": 8571.62, "pI": 4.630, "net_charge_pH7_4": -7.729},
    9: {"mw_da": 8510.52, "pI": 4.728, "net_charge_pH7_4": -5.707},
    10: {"mw_da": 8434.64, "pI": 5.205, "net_charge_pH7_4": -3.319},
}

MANAGED_BLOCK_IDS = {
    "structure-figures-finding",
    "structure-figures-grid",
    "structure-figures-grid-1",
    "structure-figures-grid-2",
    "structure-figures-grid-3",
    "structure-figures-grid-4",
    "structure-figures-grid-5",
    "structure-input-table-intro",
    "structure-input-table-block",
    "md-protocol-section",
    "md-protocol-table-block",
    "md-analysis-section",
    "spr-mst-section",
    "assay-table-block",
}

MANAGED_TABLE_IDS = {
    "structure-input-table",
    "md-protocol-table",
    "assay-starting-conditions-table",
}

MANAGED_SOURCE_IDS = {
    "src-best-complexes",
    "src-md-protocol",
    "src-assay-plan",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def data_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/jpeg;base64,{encoded}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--source-output-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def source_records(executed_at: str) -> list[dict]:
    branch_root = (
        "https://github.com/ydlongtao/usp15-dusp-binder-design"
        "/tree/agent/record-parameter-optimization-plan"
    )
    return [
        {
            "id": "src-best-complexes",
            "label": "R10 best USP15 predicted complexes",
            "path": "sources/md/structure_inputs.json",
            "href": f"{branch_root}/docs/structures/USP15_R10",
            "query": {
                "description": (
                    "每个 R10 候选的最佳 USP15 AF2 正向 seed PDB，以及由这些 "
                    "PDB 统一渲染的整体复合物和热点界面图。"
                ),
                "language": "PDB/JPEG",
                "engine": "DuckDB",
                "sql": (
                    "SELECT * FROM read_json_auto('sources/md/structure_inputs.json') "
                    "ORDER BY rank"
                ),
                "executed_at": executed_at,
                "tables_used": [
                    "r10/final_candidates/manifest.json",
                    "docs/structures/USP15_R10",
                    "docs/figures/USP15_R10",
                ],
                "filters": [
                    "one best positive USP15 seed per exported candidate",
                    "binder chain A; USP15 DUSP chain B",
                    "geometry-conditioned AF2 predictions",
                ],
                "metric_definitions": [
                    "MW, pI and pH 7.4 charge are Biopython ProteinAnalysis estimates for the untagged 76-aa sequence.",
                    "Source hotspots A50/A52/A53/A55/A57/A61 map to prediction-chain B45/B47/B48/B50/B52/B56.",
                ],
            },
        },
        {
            "id": "src-md-protocol",
            "label": "Prospective USP15 R10 OpenMM MD protocol",
            "path": "sources/md/md_protocol_table.json",
            "href": (
                "https://github.com/ydlongtao/usp15-dusp-binder-design/blob/"
                "agent/record-parameter-optimization-plan/config/"
                "usp15_r10_openmm_md.json"
            ),
            "query": {
                "description": (
                    "尚未执行的显式溶剂 OpenMM 参数、预平衡、生产采样、分析指标"
                    "和预先定义的计算分流条件。"
                ),
                "language": "JSON",
                "engine": "DuckDB",
                "sql": (
                    "SELECT * FROM read_json_auto('sources/md/md_protocol_table.json') "
                    "ORDER BY \"order\""
                ),
                "executed_at": executed_at,
                "tables_used": ["config/usp15_r10_openmm_md.json"],
                "filters": [
                    "OpenMM >=8.5",
                    "AMBER ff19SB with OPC water",
                    "3 independent 100-ns replicates per complex for initial triage",
                    "no binder-target restraints in production",
                ],
            },
        },
        {
            "id": "src-assay-plan",
            "label": "R10 SPR/MST validation starting plan",
            "path": "sources/md/assay_starting_conditions.json",
            "href": (
                "https://github.com/ydlongtao/usp15-dusp-binder-design/blob/"
                "agent/record-parameter-optimization-plan/docs/"
                "USP15_R10_STRUCTURE_MD_AND_SPR_MST_PLAN.md"
            ),
            "query": {
                "description": (
                    "与候选大小、电荷和无 Cys 特征相匹配的 SPR/MST 起始条件；"
                    "所有范围均需用真实样品优化。"
                ),
                "language": "Markdown",
                "engine": "DuckDB",
                "sql": (
                    "SELECT * FROM "
                    "read_json_auto('sources/md/assay_starting_conditions.json') "
                    "ORDER BY \"order\""
                ),
                "executed_at": executed_at,
                "tables_used": [
                    "docs/USP15_R10_STRUCTURE_MD_AND_SPR_MST_PLAN.md"
                ],
                "filters": [
                    "tag-free binder sequence unless explicitly stated",
                    "parallel USP4/USP11 counterscreen",
                    "starting conditions, not validated assay parameters",
                ],
            },
        },
    ]


def image_pair_blocks(figures: list[tuple[int, int, Path]]) -> list[dict]:
    blocks = []
    for index in range(0, len(figures), 2):
        pair = figures[index : index + 2]
        images = []
        for rank, _seed, path in pair:
            images.append(
                f'<img src="{data_uri(path)}" loading="lazy" decoding="async" '
                f'alt="USP15 DUSP 与 Rank {rank:02d} binder 的预测复合物和热点界面图" '
                'style="display:block;height:220px;width:auto;max-width:49%;'
                'object-fit:contain;border-radius:8px;">'
            )
        blocks.append(
            {
                "id": f"structure-figures-grid-{index // 2 + 1}",
                "type": "html",
                "body": (
                    '<div style="height:224px;display:flex;align-items:center;'
                    'justify-content:center;gap:10px;overflow:hidden;background:#fff;">'
                    + "".join(images)
                    + "</div>"
                ),
            }
        )
    return blocks


def main() -> None:
    args = parse_args()
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    export_manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    records = sorted(export_manifest["records"], key=lambda item: item["rank"])
    if len(records) != 10:
        raise ValueError(f"Expected 10 candidates, found {len(records)}")

    structure_rows = []
    figures = []
    for record in records:
        rank = int(record["rank"])
        seed = int(record["best_positive_seed"])
        pdb_name = f"USP15_R10_rank{rank:02d}_best_seed{seed}.pdb"
        figure_name = f"USP15_R10_rank{rank:02d}_complex.jpg"
        pdb_path = args.structures / pdb_name
        figure_path = args.figures / figure_name
        if not pdb_path.is_file() or not figure_path.is_file():
            raise FileNotFoundError(f"Missing rank {rank} PDB or figure")
        properties = PROPERTIES[rank]
        full_pdb_sha = sha256(pdb_path)
        full_figure_sha = sha256(figure_path)
        structure_rows.append(
            {
                "rank": rank,
                "best_seed": seed,
                **properties,
                "pdb_file": pdb_name,
                "pdb_sha256": full_pdb_sha,
                "pdb_sha12": full_pdb_sha[:12],
                "figure_file": figure_name,
                "figure_sha256": full_figure_sha,
            }
        )
        figures.append((rank, seed, figure_path))

    md_protocol_rows = [
        {"order": 1, "stage": "Preparation", "parameter": "Input", "value": "best USP15 AF2 seed; chain A binder + chain B DUSP", "status": "fixed"},
        {"order": 2, "stage": "Preparation", "parameter": "Protonation", "value": "pH 7.4; record His states per system", "status": "prospective"},
        {"order": 3, "stage": "Model", "parameter": "Force field / water", "value": "AMBER ff19SB / OPC", "status": "fixed"},
        {"order": 4, "stage": "Solvent", "parameter": "Box / salt", "value": "dodecahedron; 1.2 nm padding; 0.15 M NaCl", "status": "fixed"},
        {"order": 5, "stage": "Nonbonded", "parameter": "Electrostatics / cutoff", "value": "PME / 1.0 nm", "status": "fixed"},
        {"order": 6, "stage": "Integrator", "parameter": "Temperature / step", "value": "LangevinMiddle 300 K; 1 ps⁻¹; 2 fs", "status": "fixed"},
        {"order": 7, "stage": "Pressure", "parameter": "Barostat", "value": "MonteCarloBarostat; 1 bar; every 25 steps", "status": "fixed"},
        {"order": 8, "stage": "Minimization", "parameter": "Tolerance / iterations", "value": "10 kJ mol⁻¹ nm⁻¹ / 20,000 max", "status": "fixed"},
        {"order": 9, "stage": "Equilibration", "parameter": "NVT", "value": "0.5 ns; heavy-atom restraint 1000 kJ mol⁻¹ nm⁻²", "status": "fixed"},
        {"order": 10, "stage": "Equilibration", "parameter": "NPT", "value": "0.5 ns at 100 + 1.0 ns at 10 kJ mol⁻¹ nm⁻²", "status": "fixed"},
        {"order": 11, "stage": "Production", "parameter": "Initial sampling", "value": "3 independent replicas × 100 ns per complex", "status": "not run"},
        {"order": 12, "stage": "Production", "parameter": "Selected extension", "value": "extend each selected replica to 500 ns total", "status": "not run"},
        {"order": 13, "stage": "Output", "parameter": "Trajectory / state / checkpoint", "value": "10 ps / 10 ps / 1 ns", "status": "fixed"},
        {"order": 14, "stage": "Analysis", "parameter": "Burn-in", "value": "first 20 ns of each production replica", "status": "prospective"},
        {"order": 15, "stage": "Analysis", "parameter": "Primary readouts", "value": "binder RMSD; hotspot occupancy; buried SASA; contacts; H-bonds", "status": "prospective"},
    ]

    assay_rows = [
        {"order": 1, "method": "SPR", "parameter": "Orientation", "starting_condition": "capture/immobilize USP15 DUSP; binder as analyte", "purpose": "reduce small-analyte response and preserve binder sequence"},
        {"order": 2, "method": "SPR", "parameter": "Running buffer", "starting_condition": "10 mM HEPES pH 7.4, 150 mM NaCl, 0.05% Tween-20", "purpose": "starting condition; optimize against target stability"},
        {"order": 3, "method": "SPR", "parameter": "Binder series", "starting_condition": "0.5 nM–10 μM, 12–16 two-fold points", "purpose": "cover unknown affinity before narrowing range"},
        {"order": 4, "method": "SPR", "parameter": "Contact time", "starting_condition": "120 s association; 300–600 s dissociation", "purpose": "support steady-state and kinetic inspection"},
        {"order": 5, "method": "SPR", "parameter": "Surface / controls", "starting_condition": "about 100–500 RU; blank surface; double reference", "purpose": "limit mass transport and nonspecific signal"},
        {"order": 6, "method": "MST", "parameter": "Label strategy", "starting_condition": "label USP15 DUSP; avoid altering the no-Cys binder", "purpose": "retain the designed binder sequence"},
        {"order": 7, "method": "MST", "parameter": "Labeled target", "starting_condition": "10–50 nM after fluorescence/capillary scan", "purpose": "stay within instrument signal window"},
        {"order": 8, "method": "MST", "parameter": "Binder series", "starting_condition": "0.3 nM–10 μM, 16 two-fold points", "purpose": "broad first-pass KD coverage"},
        {"order": 9, "method": "MST", "parameter": "Buffer", "starting_condition": "10 mM HEPES pH 7.4, 150 mM NaCl, 0.01–0.05% Tween-20", "purpose": "reduce capillary adsorption; add 0.05% BSA only if needed"},
        {"order": 10, "method": "Both", "parameter": "Counterscreen", "starting_condition": "USP4 and USP11 DUSP under matched conditions", "purpose": "test whether computational selectivity transfers experimentally"},
    ]

    args.source_output_dir.mkdir(parents=True, exist_ok=True)
    for filename, rows in [
        ("structure_inputs.json", structure_rows),
        ("md_protocol_table.json", md_protocol_rows),
        ("assay_starting_conditions.json", assay_rows),
    ]:
        (args.source_output_dir / filename).write_text(
            json.dumps(rows, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

    snapshot = artifact["snapshot"]
    snapshot["datasets"]["structure_inputs"] = structure_rows
    snapshot["datasets"]["md_protocol"] = md_protocol_rows
    snapshot["datasets"]["assay_starting_conditions"] = assay_rows

    manifest = artifact["manifest"]
    manifest["tables"] = [
        table for table in manifest["tables"] if table["id"] not in MANAGED_TABLE_IDS
    ]
    manifest["tables"].extend(
        [
            {
                "id": "structure-input-table",
                "title": "十个结构图与 MD 初始 PDB",
                "subtitle": "每个候选使用 R10 manifest 指定的最佳 USP15 正向 seed；理化值为未加标签序列估算",
                "dataset": "structure_inputs",
                "sourceId": "src-best-complexes",
                "density": "dense",
                "defaultSort": {"field": "rank", "direction": "asc"},
                "columns": [
                    {"field": "rank", "label": "Rank", "format": "number"},
                    {"field": "best_seed", "label": "Best seed", "format": "number"},
                    {"field": "mw_da", "label": "MW (Da)", "format": "number"},
                    {"field": "pI", "label": "pI", "format": "number"},
                    {"field": "net_charge_pH7_4", "label": "Charge at pH 7.4", "format": "number"},
                    {"field": "pdb_file", "label": "PDB"},
                    {"field": "pdb_sha12", "label": "PDB SHA-256 (12)"},
                ],
            },
            {
                "id": "md-protocol-table",
                "title": "预先定义的 OpenMM MD 参数",
                "subtitle": "这些参数尚未执行；表中 not run 不得被解释为 MD 结果",
                "dataset": "md_protocol",
                "sourceId": "src-md-protocol",
                "density": "dense",
                "defaultSort": {"field": "order", "direction": "asc"},
                "columns": [
                    {"field": "order", "label": "#", "format": "number"},
                    {"field": "stage", "label": "Stage"},
                    {"field": "parameter", "label": "Parameter"},
                    {"field": "value", "label": "Fixed/proposed value"},
                    {"field": "status", "label": "Status"},
                ],
            },
            {
                "id": "assay-starting-conditions-table",
                "title": "SPR 与 MST 起始条件",
                "subtitle": "用于首轮方法开发；真实浓度窗口、固定密度和添加剂必须按样品行为优化",
                "dataset": "assay_starting_conditions",
                "sourceId": "src-assay-plan",
                "density": "dense",
                "defaultSort": {"field": "order", "direction": "asc"},
                "columns": [
                    {"field": "order", "label": "#", "format": "number"},
                    {"field": "method", "label": "Method"},
                    {"field": "parameter", "label": "Parameter"},
                    {"field": "starting_condition", "label": "Starting condition"},
                    {"field": "purpose", "label": "Purpose / caveat"},
                ],
            },
        ]
    )

    executed_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )
    new_sources = source_records(executed_at)
    manifest["sources"] = [
        source
        for source in manifest["sources"]
        if source["id"] not in MANAGED_SOURCE_IDS
    ] + new_sources
    artifact["sources"] = [
        source
        for source in artifact.get("sources", [])
        if source["id"] not in MANAGED_SOURCE_IDS
    ]
    artifact["sources"].extend(
        {
            "id": source["id"],
            "label": source["label"],
            "path": source["path"],
            **({"href": source["href"]} if source.get("href") else {}),
            "query": {
                key: value
                for key, value in source["query"].items()
                if key
                in {
                    "description",
                    "language",
                    "engine",
                    "sql",
                    "executed_at",
                    "tables_used",
                    "filters",
                    "metric_definitions",
                }
            },
        }
        for source in new_sources
    )

    blocks = [
        block for block in manifest["blocks"] if block["id"] not in MANAGED_BLOCK_IDS
    ]
    insertion_index = next(
        index
        for index, block in enumerate(blocks)
        if block["id"] == "candidate-table-block"
    ) + 1
    new_blocks = [
        {
            "id": "structure-figures-finding",
            "type": "markdown",
            "sourceId": "src-best-complexes",
            "body": (
                "## 十个候选均定位于同一 USP15 DUSP 热点表面\n\n"
                "**下列十张图使用每个候选的最佳 USP15 正向 AF2 seed。** 左图显示"
                "整体 DUSP–binder 构象，右图放大六个预定义热点。预测文件中 chain A "
                "是 76-aa binder，chain B 是 129-aa USP15 DUSP；3T9L 源热点 "
                "A50/A52/A53/A55/A57/A61 对应预测链 B45/B47/B48/B50/B52/B56。\n\n"
                "所有图片使用同一渲染协议，因此可比较姿势与界面位置；但它们是 "
                "interface-template 条件化结构预测，不是 MD 轨迹快照或实验结构，"
                "图片外观不能用于推断 KD。"
            ),
        },
        *image_pair_blocks(figures),
        {
            "id": "structure-input-table-intro",
            "type": "markdown",
            "sourceId": "src-best-complexes",
            "body": (
                "## 可下载 PDB 与序列理化值支持后续建模\n\n"
                "每个 PDB 均来自 R10 manifest 的最佳 USP15 seed，并保留在 GitHub "
                "`docs/structures/USP15_R10`。分子量、理论 pI 和 pH 7.4 净电荷由 "
                "Biopython `ProteinAnalysis` 对未加标签序列计算；实验标签、末端修饰"
                "和缓冲条件会改变实际值。"
            ),
        },
        {
            "id": "structure-input-table-block",
            "type": "table",
            "tableId": "structure-input-table",
        },
        {
            "id": "md-protocol-section",
            "type": "markdown",
            "sourceId": "src-md-protocol",
            "body": (
                "## MD 参数已预先定义，但尚未产生轨迹结果\n\n"
                "**建议协议使用 OpenMM ≥8.5、AMBER ff19SB 和 OPC 显式水。** 每个"
                "复合物先运行 3 个独立的 100-ns NPT 重复，总初筛采样为 3.0 μs；"
                "实验优先候选的每个重复再延长到总计 500 ns。生产阶段不得保留 "
                "binder–target 距离或位置约束。\n\n"
                "下面的 `not run` 是关键状态：报告目前只有初始结构和拟定参数，没有"
                "轨迹、RMSD、接触占有率、自由能或 MD 稳定性结论。"
            ),
        },
        {
            "id": "md-protocol-table-block",
            "type": "table",
            "tableId": "md-protocol-table",
        },
        {
            "id": "md-analysis-section",
            "type": "markdown",
            "sourceId": "src-md-protocol",
            "body": (
                "## 后续 MD 应按重复报告界面保持性，而不是只给一条 RMSD\n\n"
                "生产轨迹前 20 ns 作为 burn-in。以 USP15 DUSP chain B Cα 对齐后，"
                "分别计算 binder Cα RMSD、DUSP/binder RMSF、质心距离、4.5 Å 接触、"
                "六热点占有率、氢键、盐桥、buried SASA 和二级结构保持率。每个重复"
                "单独展示，并报告中位数和 95% bootstrap CI。\n\n"
                "预先定义的计算分流条件是：无持续分离、binder RMSD 中位数 ≤2.5 Å、"
                "至少 70% 帧 RMSD ≤3.0 Å、至少 4/6 热点接触占有率 ≥50%、buried "
                "SASA 中位数 ≥600 Å²，且三个重复结论一致。这些阈值仅用于安排实验"
                "顺序，未经 SPR/MST 校准，不能作为结合阳性的证据。MM/GBSA 只能作"
                "相对诊断，不能直接换算 KD。"
            ),
        },
        {
            "id": "spr-mst-section",
            "type": "markdown",
            "sourceId": "src-assay-plan",
            "body": (
                "## SPR/MST 首轮应标记或固定 USP15 DUSP，并平行反筛 USP4/USP11\n\n"
                "binder 仅约 8.3–8.7 kDa 且全部无 Cys。SPR 建议低密度、定向捕获 "
                "USP15 DUSP，把 binder 作为 analyte；MST 建议标记 USP15 DUSP，"
                "避免为了检测临时改造原始 binder。先用宽浓度范围寻找浓度依赖，再"
                "根据真实 KD 和非特异吸附缩窄范围。\n\n"
                "下面条件是方法开发起点而非已验证 SOP。USP4 和 USP11 DUSP 必须在"
                "相同缓冲液、标记/固定策略和浓度范围下进行 counterscreen；不要只"
                "报告单点信号。"
            ),
        },
        {
            "id": "assay-table-block",
            "type": "table",
            "tableId": "assay-starting-conditions-table",
        },
    ]
    blocks[insertion_index:insertion_index] = new_blocks
    manifest["blocks"] = blocks

    for block in manifest["blocks"]:
        if block["id"] == "technical-summary":
            addition = (
                "\n\n**本版新增：** 十个最佳正向 seed 的整体/热点结构图、可下载 "
                "PDB、候选理化值、尚未执行的显式溶剂 MD 参数，以及 SPR/MST 起始"
                "条件。新增内容不改变 R10 的 geometry-conditioned 解释边界。"
            )
            if "**本版新增：**" not in block["body"]:
                block["body"] += addition

    args.output.write_text(
        json.dumps(artifact, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
