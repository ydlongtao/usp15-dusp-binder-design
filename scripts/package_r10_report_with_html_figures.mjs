#!/usr/bin/env node

/**
 * Package the R10 portable report when it contains sandboxed HTML image blocks.
 *
 * The canonical delivery helper currently waits for window.load before booting
 * the enhanced reader. Its own lazy srcdoc iframe fallback for an HTML block can
 * keep that event from reaching the chart-extraction startup window. This
 * wrapper still uses the canonical package builder, chart extractor, and
 * structural verifier:
 *
 * 1. extract static chart SVGs from an otherwise identical artifact with only
 *    the HTML image blocks omitted;
 * 2. rebuild the complete artifact, including all image blocks, with those
 *    canonical static charts;
 * 3. verify exact embedded-payload equality and all required portable roots.
 *
 * Browser QA of the HTML-block report remains a documented renderer limitation;
 * the semantic reader stays complete and the high-resolution JPEGs/PDBs remain
 * separate tracked deliverables.
 */

import { mkdtemp, readFile, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { pathToFileURL } from "node:url";
import process from "node:process";

function parseArgs(argv) {
  const values = {};
  for (let index = 2; index < argv.length; index += 2) {
    const key = argv[index];
    const value = argv[index + 1];
    if (!key?.startsWith("--") || value === undefined) {
      throw new Error(`Invalid argument near ${key ?? "<end>"}`);
    }
    values[key.slice(2)] = value;
  }
  for (const required of ["input", "output", "plugin-root"]) {
    if (!values[required]) throw new Error(`Missing --${required}`);
  }
  return values;
}

async function importScript(pluginRoot, name) {
  const path = join(
    pluginRoot,
    "skills",
    "build-report",
    "scripts",
    name,
  );
  return import(pathToFileURL(path).href);
}

async function main() {
  const args = parseArgs(process.argv);
  const inputPath = resolve(args.input);
  const outputPath = resolve(args.output);
  const pluginRoot = resolve(args["plugin-root"]);
  const artifact = JSON.parse(await readFile(inputPath, "utf8"));

  const { buildPortableArtifact } = await importScript(
    pluginRoot,
    "build_portable_artifact.mjs",
  );
  const { extractPortableChartSvgs } = await importScript(
    pluginRoot,
    "extract_portable_chart_svgs.mjs",
  );
  const { verifyPortableArtifactStructure } = await importScript(
    pluginRoot,
    "verify_portable_artifact.mjs",
  );

  const extractionArtifact = structuredClone(artifact);
  extractionArtifact.manifest.blocks =
    extractionArtifact.manifest.blocks.filter((block) => block.type !== "html");

  const temporaryDirectory = await mkdtemp(join(tmpdir(), "usp15-r10-report-"));
  const extractionHtml = join(temporaryDirectory, "chart-extraction.html");
  try {
    await writeFile(
      extractionHtml,
      buildPortableArtifact(extractionArtifact),
      "utf8",
    );
    const staticCharts = await extractPortableChartSvgs({
      htmlPath: extractionHtml,
      readyTimeoutMs: 10_000,
      actionTimeoutMs: 5_000,
    });
    await writeFile(
      outputPath,
      buildPortableArtifact(artifact, { staticCharts }),
      "utf8",
    );
    const verification = verifyPortableArtifactStructure({
      artifactPath: inputPath,
      htmlPath: outputPath,
    });
    process.stdout.write(
      `${JSON.stringify({
        ok: true,
        html: outputPath,
        stages: {
          validation: "passed",
          canonical_static_charts: "passed",
          structural_verification: "passed",
          enhanced_reader_browser_qa: "blocked_by_html_block_loader",
        },
        counts: verification.counts,
      })}\n`,
    );
  } finally {
    await rm(temporaryDirectory, { recursive: true, force: true });
  }
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
