#!/usr/bin/env node

/**
 * Render consistent two-panel USP15 DUSP–binder figures from AF2 PDB files.
 *
 * Requirements:
 *   - Node.js
 *   - Playwright with an installed Chromium
 *   - a local copy of ngl.js (NGL 2.x)
 *
 * The AF2 exports use chain A for the 76-aa binder and chain B for the
 * 129-aa USP15 DUSP target. Source 3T9L residues A50/A52/A53/A55/A57/A61
 * therefore correspond to AF2 output residues B45/B47/B48/B50/B52/B56.
 */

import { createServer } from "node:http";
import { readFile, readdir, mkdir } from "node:fs/promises";
import { extname, join, resolve } from "node:path";
import process from "node:process";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const { chromium } = require("playwright");

function parseArgs(argv) {
  const args = {};
  for (let index = 2; index < argv.length; index += 2) {
    const key = argv[index];
    const value = argv[index + 1];
    if (!key?.startsWith("--") || value === undefined) {
      throw new Error(`Invalid argument near ${key ?? "<end>"}`);
    }
    args[key.slice(2)] = value;
  }
  for (const required of ["input-dir", "output-dir", "ngl"]) {
    if (!args[required]) throw new Error(`Missing --${required}`);
  }
  return args;
}

function contentType(path) {
  return {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".pdb": "chemical/x-pdb",
  }[extname(path)] ?? "application/octet-stream";
}

function rankAndSeed(filename) {
  const match = filename.match(/rank(\d+)_best_seed(\d+)\.pdb$/);
  if (!match) throw new Error(`Unexpected PDB filename: ${filename}`);
  return { rank: Number(match[1]), seed: Number(match[2]) };
}

function pageHtml({ pdbFile, rank, seed }) {
  const hotspotSelection = "45,47,48,50,52,56:B";
  return `<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>USP15 R10 Rank ${String(rank).padStart(2, "0")}</title>
  <script src="/ngl.js"></script>
  <style>
    * { box-sizing: border-box; }
    html, body { margin: 0; width: 1400px; height: 900px; overflow: hidden; }
    body {
      background: #f5f7fa;
      color: #162033;
      font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }
    header {
      height: 96px;
      padding: 18px 28px 12px;
      background: #ffffff;
      border-bottom: 1px solid #dbe2ea;
      display: flex;
      align-items: center;
      justify-content: space-between;
    }
    h1 { margin: 0 0 5px; font-size: 26px; font-weight: 720; letter-spacing: -0.02em; }
    .subtitle { color: #607089; font-size: 15px; }
    .badge {
      padding: 8px 12px;
      border-radius: 999px;
      background: #fff1e8;
      color: #9a4318;
      border: 1px solid #f0c5a9;
      font: 700 14px ui-monospace, SFMono-Regular, Menlo, monospace;
    }
    main {
      height: 724px;
      padding: 18px 22px;
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 18px;
    }
    .panel {
      min-width: 0;
      background: #ffffff;
      border: 1px solid #dbe2ea;
      border-radius: 14px;
      overflow: hidden;
      box-shadow: 0 6px 20px rgba(38, 55, 77, 0.07);
    }
    .panel-title {
      height: 48px;
      padding: 14px 18px;
      border-bottom: 1px solid #e7ecf2;
      font-size: 15px;
      font-weight: 700;
    }
    .viewport { width: 100%; height: 638px; }
    footer {
      height: 80px;
      padding: 14px 28px;
      background: #ffffff;
      border-top: 1px solid #dbe2ea;
      display: flex;
      justify-content: space-between;
      align-items: flex-start;
      gap: 20px;
      font-size: 14px;
      color: #52627a;
    }
    .legend { display: flex; flex-wrap: wrap; gap: 16px; }
    .legend-item { display: flex; align-items: center; gap: 7px; white-space: nowrap; }
    .swatch { width: 14px; height: 14px; border-radius: 4px; border: 1px solid rgba(0,0,0,.15); }
    .note { max-width: 630px; line-height: 1.45; text-align: right; }
  </style>
</head>
<body>
  <header>
    <div>
      <h1>USP15 DUSP–Binder 预测复合物 · Rank ${String(rank).padStart(2, "0")}</h1>
      <div class="subtitle">AF2 model_2_ptm interface-template · best positive seed ${seed} · geometry-conditioned</div>
    </div>
    <div class="badge">76 aa binder</div>
  </header>
  <main>
    <section class="panel">
      <div class="panel-title">整体结合构象</div>
      <div id="overview" class="viewport"></div>
    </section>
    <section class="panel">
      <div class="panel-title">DUSP 热点界面近景</div>
      <div id="interface" class="viewport"></div>
    </section>
  </main>
  <footer>
    <div class="legend">
      <span class="legend-item"><span class="swatch" style="background:#6d8dad"></span>USP15 DUSP</span>
      <span class="legend-item"><span class="swatch" style="background:#e36b35"></span>Binder</span>
      <span class="legend-item"><span class="swatch" style="background:#d7266d"></span>6 个设计热点</span>
    </div>
    <div class="note">热点按 3T9L 源编号为 A50/A52/A53/A55/A57/A61；在本 AF2 输出的目标链 B 中对应 B45/B47/B48/B50/B52/B56。该图是模板条件化结构预测，不是 MD 轨迹快照或实验结构。</div>
  </footer>
  <script>
    const pdbUrl = "/pdb/${encodeURIComponent(pdbFile)}";
    const hotspotSelection = ${JSON.stringify(hotspotSelection)};
    const stageOptions = {
      backgroundColor: "white",
      cameraType: "orthographic",
      quality: "high",
      sampleLevel: 2
    };

    function addCoreRepresentations(component, closeup) {
      component.addRepresentation("cartoon", {
        sele: ":B",
        colorScheme: "uniform",
        colorValue: "#6d8dad",
        opacity: closeup ? 0.55 : 1.0,
        smoothSheet: true
      });
      component.addRepresentation("cartoon", {
        sele: ":A",
        colorScheme: "uniform",
        colorValue: "#e36b35",
        smoothSheet: true
      });
      component.addRepresentation("ball+stick", {
        sele: hotspotSelection,
        colorScheme: "uniform",
        colorValue: "#d7266d",
        multipleBond: "symmetric",
        radiusScale: 1.15,
        aspectRatio: 1.8
      });
      if (closeup) {
        component.addRepresentation("surface", {
          sele: ":B",
          colorScheme: "uniform",
          colorValue: "#8fa9c1",
          opacity: 0.28,
          surfaceType: "av",
          probeRadius: 1.4,
          scaleFactor: 1.5
        });
        component.addRepresentation("licorice", {
          sele: ":A",
          colorScheme: "element",
          radiusScale: 0.55,
          multipleBond: "symmetric"
        });
        component.addRepresentation("label", {
          sele: hotspotSelection + " and .CA",
          colorScheme: "uniform",
          colorValue: "#8f164a",
          labelType: "resno",
          labelGrouping: "residue",
          showBackground: true,
          backgroundColor: "white",
          backgroundOpacity: 0.78,
          borderColor: "#8f164a",
          borderWidth: 0.4,
          fontWeight: "bold",
          radiusScale: 0.65
        });
      }
    }

    async function buildStage(elementId, closeup) {
      const stage = new NGL.Stage(elementId, stageOptions);
      const component = await stage.loadFile(pdbUrl, { ext: "pdb" });
      addCoreRepresentations(component, closeup);
      if (closeup) {
        component.autoView(":A or " + hotspotSelection, 0);
      } else {
        component.autoView("protein", 0);
      }
      stage.viewer.requestRender();
      return stage;
    }

    Promise.all([
      buildStage("overview", false),
      buildStage("interface", true)
    ]).then(() => {
      setTimeout(() => {
        document.documentElement.dataset.ready = "true";
      }, 1800);
    }).catch((error) => {
      document.body.textContent = String(error?.stack || error);
      document.documentElement.dataset.ready = "error";
    });
  </script>
</body>
</html>`;
}

async function main() {
  const args = parseArgs(process.argv);
  const inputDir = resolve(args["input-dir"]);
  const outputDir = resolve(args["output-dir"]);
  const nglPath = resolve(args.ngl);
  await mkdir(outputDir, { recursive: true });

  const pdbFiles = (await readdir(inputDir))
    .filter((name) => name.endsWith(".pdb"))
    .sort();
  if (!pdbFiles.length) throw new Error("No PDB files found");

  const server = createServer(async (request, response) => {
    try {
      const url = new URL(request.url, "http://127.0.0.1");
      if (url.pathname === "/ngl.js") {
        const body = await readFile(nglPath);
        response.writeHead(200, { "Content-Type": contentType(nglPath) });
        response.end(body);
        return;
      }
      if (url.pathname.startsWith("/pdb/")) {
        const name = decodeURIComponent(url.pathname.slice("/pdb/".length));
        if (!pdbFiles.includes(name)) throw new Error("Unknown PDB");
        const body = await readFile(join(inputDir, name));
        response.writeHead(200, { "Content-Type": "chemical/x-pdb" });
        response.end(body);
        return;
      }
      if (url.pathname === "/") {
        const name = url.searchParams.get("pdb");
        if (!name || !pdbFiles.includes(name)) throw new Error("Missing PDB");
        const { rank, seed } = rankAndSeed(name);
        response.writeHead(200, { "Content-Type": "text/html; charset=utf-8" });
        response.end(pageHtml({ pdbFile: name, rank, seed }));
        return;
      }
      response.writeHead(404);
      response.end("Not found");
    } catch (error) {
      response.writeHead(500, { "Content-Type": "text/plain; charset=utf-8" });
      response.end(String(error?.stack || error));
    }
  });

  await new Promise((resolveListen) => server.listen(0, "127.0.0.1", resolveListen));
  const address = server.address();
  const origin = `http://127.0.0.1:${address.port}`;

  const browser = await chromium.launch({
    headless: true,
    executablePath: process.env.CHROME_EXECUTABLE || undefined,
    args: [
      "--disable-dev-shm-usage",
      "--enable-webgl",
      "--ignore-gpu-blocklist",
      "--use-angle=swiftshader",
    ],
  });

  try {
    const page = await browser.newPage({
      viewport: { width: 1400, height: 900 },
      deviceScaleFactor: 1,
    });
    for (const pdbFile of pdbFiles) {
      const { rank, seed } = rankAndSeed(pdbFile);
      await page.goto(`${origin}/?pdb=${encodeURIComponent(pdbFile)}`, {
        waitUntil: "networkidle",
      });
      await page.waitForFunction(
        () => document.documentElement.dataset.ready === "true",
        undefined,
        { timeout: 45_000 },
      );
      const output = join(
        outputDir,
        `USP15_R10_rank${String(rank).padStart(2, "0")}_complex.jpg`,
      );
      await page.screenshot({
        path: output,
        type: "jpeg",
        quality: 84,
        fullPage: false,
      });
      process.stdout.write(`rendered rank ${rank} seed ${seed}: ${output}\n`);
    }
  } finally {
    await browser.close();
    await new Promise((resolveClose) => server.close(resolveClose));
  }
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
