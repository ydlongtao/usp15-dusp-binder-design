nextflow.enable.dsl=2

/*
 * Minimal OVO-facing ESM3 descriptor workflow.
 * The container is intentionally external to the OVO Python environment.
 */
params.input_fasta = null
params.output_dir = 'esm3_eval_results'
params.adapter_script = 'scripts/run_esm3_eval.py'
params.esm3_site = null
params.esm3_cache = null
params.container = 'ovo-esm:latest'
params.model = 'esm3-sm-open-v1'
params.num_steps = 1

process ESM3_EVAL {
    tag "${candidate_id}"
    publishDir params.output_dir, mode: 'copy', saveAs: { filename -> "${candidate_id}/${filename}" }

    input:
    tuple val(candidate_id), path(fasta)
    path adapter

    output:
    path 'esm3_eval.json', emit: json
    path 'esm3_metrics.csv', emit: metrics
    path 'esm3_structure.pdb', emit: structure

    script:
    def site_mount = params.esm3_site ? "-v ${params.esm3_site}:/esm3_site:ro" : ''
    def cache_mount = params.esm3_cache ? "-v ${params.esm3_cache}:/esm3_cache" : ''
    """
    mkdir -p out
    cp -L ${fasta.getName()} input.fasta
    cp -L ${adapter.getName()} adapter.py
    docker run --rm --gpus all \\
      ${site_mount} ${cache_mount} \\
      -v \"\$PWD:/work\" -e PYTHONPATH=/esm3_site -e HF_HOME=/esm3_cache \\
      ${params.container} /usr/bin/python3 /work/adapter.py \\
      --fasta /work/input.fasta --candidate-id ${candidate_id} \\
      --output-dir /work/out --model ${params.model} --num-steps ${params.num_steps}
    cp out/esm3_eval.json esm3_eval.json
    cp out/esm3_metrics.csv esm3_metrics.csv
    cp out/esm3_structure.pdb esm3_structure.pdb
    """
}

workflow {
    if (!params.input_fasta) {
        error 'Set --input_fasta to a FASTA file'
    }
    fasta_ch = Channel.fromPath(params.input_fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    adapter_ch = Channel.fromPath(params.adapter_script, checkIfExists: true)
    ESM3_EVAL(fasta_ch, adapter_ch)
}
