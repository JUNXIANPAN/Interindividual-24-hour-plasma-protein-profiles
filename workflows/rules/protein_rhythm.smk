"""Deterministic protein-rhythm screening and internal-validation workflow."""

import sys


PR = config.get("protein_rhythm", {})
PR_MATRIX = PR.get("input_matrix", "data/report.pg_matrix.tsv")
PR_OUT = PR.get("output_dir", "outputs/snakemake_protein_rhythm")
PR_PYTHON = PR.get("python_executable", sys.executable)


rule protein_rhythm_all:
    input:
        f"{PR_OUT}/result_manifest.json",
        f"{PR_OUT}/plots/five_candidates_individual_cosinor_fits.png",
        f"{PR_OUT}/plots/five_candidates_centered_trajectories.png",


rule protein_rhythm_screen:
    input:
        matrix=PR_MATRIX,
        script="analysis/simplified_mixed_cosinor.py",
        qc_kernel="analysis/first_pass_mixed_cosinor.py",
    output:
        results=f"{PR_OUT}/screen/simplified_mixed_cosinor_results.tsv",
        manifest=f"{PR_OUT}/screen/sample_manifest.tsv",
        protein_qc=f"{PR_OUT}/screen/protein_qc.tsv",
        summary=f"{PR_OUT}/screen/summary.json",
    log:
        f"{PR_OUT}/logs/screen.log",
    params:
        outdir=f"{PR_OUT}/screen",
        python=PR_PYTHON,
        bootstrap_top=lambda wildcards: PR.get("screen_bootstrap_top", 10),
        bootstrap_reps=lambda wildcards: PR.get("screen_bootstrap_reps", 50),
        seed=lambda wildcards: PR.get("screen_seed", 20250902),
    shell:
        """
        {params.python:q} {input.script:q} \
          --input {input.matrix:q} \
          --output-dir {params.outdir:q} \
          --bootstrap-top {params.bootstrap_top} \
          --bootstrap-reps {params.bootstrap_reps} \
          --seed {params.seed} > {log:q} 2>&1
        """


rule protein_rhythm_validate_candidates:
    input:
        matrix=PR_MATRIX,
        screen=f"{PR_OUT}/screen/simplified_mixed_cosinor_results.tsv",
        script="analysis/validate_simplified_candidates.py",
        model_kernel="analysis/simplified_mixed_cosinor.py",
    output:
        summary_table=f"{PR_OUT}/validation/candidate_validation_summary.tsv",
        loo=f"{PR_OUT}/validation/leave_one_subject_out.tsv",
        origin=f"{PR_OUT}/validation/time_origin_sensitivity.tsv",
        summary=f"{PR_OUT}/validation/summary.json",
    log:
        f"{PR_OUT}/logs/validation.log",
    params:
        outdir=f"{PR_OUT}/validation",
        python=PR_PYTHON,
        bootstrap_reps=lambda wildcards: PR.get("validation_bootstrap_reps", 200),
        seed=lambda wildcards: PR.get("validation_seed", 20250903),
    shell:
        """
        {params.python:q} {input.script:q} \
          --matrix {input.matrix:q} \
          --screen-results {input.screen:q} \
          --output-dir {params.outdir:q} \
          --bootstrap-reps {params.bootstrap_reps} \
          --seed {params.seed} > {log:q} 2>&1
        """


rule protein_rhythm_candidate_plots:
    input:
        matrix=PR_MATRIX,
        manifest=f"{PR_OUT}/screen/sample_manifest.tsv",
        script="analysis/plot_candidate_trajectories.py",
    output:
        data=f"{PR_OUT}/plots/candidate_trajectory_data.tsv",
        observed=f"{PR_OUT}/plots/five_candidates_centered_trajectories.png",
        fitted=f"{PR_OUT}/plots/five_candidates_individual_cosinor_fits.png",
        dcd=f"{PR_OUT}/plots/DCD_individual_trajectories.png",
        apoc2=f"{PR_OUT}/plots/APOC2_individual_trajectories.png",
        adipoq=f"{PR_OUT}/plots/ADIPOQ_individual_trajectories.png",
        apoa4=f"{PR_OUT}/plots/APOA4_individual_trajectories.png",
        lman2=f"{PR_OUT}/plots/LMAN2_individual_trajectories.png",
    log:
        f"{PR_OUT}/logs/plots.log",
    params:
        outdir=f"{PR_OUT}/plots",
        python=PR_PYTHON,
    shell:
        """
        MPLCONFIGDIR={PR_OUT}/.matplotlib {params.python:q} {input.script:q} \
          --matrix {input.matrix:q} \
          --manifest {input.manifest:q} \
          --output-dir {params.outdir:q} > {log:q} 2>&1
        """


rule protein_rhythm_result_manifest:
    input:
        matrix=PR_MATRIX,
        screen_results=f"{PR_OUT}/screen/simplified_mixed_cosinor_results.tsv",
        screen_summary=f"{PR_OUT}/screen/summary.json",
        validation_results=f"{PR_OUT}/validation/candidate_validation_summary.tsv",
        validation_summary=f"{PR_OUT}/validation/summary.json",
        fitted_plot=f"{PR_OUT}/plots/five_candidates_individual_cosinor_fits.png",
        script="analysis/build_protein_rhythm_manifest.py",
    output:
        f"{PR_OUT}/result_manifest.json",
    params:
        run_id=lambda wildcards: PR.get("run_id", "protein-rhythm-run"),
        job_id=lambda wildcards: PR.get("job_id", "mixed-cosinor-validation"),
        protocol=lambda wildcards: PR.get("protocol_version", "1.0.0"),
        python=PR_PYTHON,
    shell:
        """
        {params.python:q} {input.script:q} \
          --run-id {params.run_id:q} \
          --job-id {params.job_id:q} \
          --protocol-version {params.protocol:q} \
          --input-matrix {input.matrix:q} \
          --screen-results {input.screen_results:q} \
          --screen-summary {input.screen_summary:q} \
          --validation-results {input.validation_results:q} \
          --validation-summary {input.validation_summary:q} \
          --fitted-plot {input.fitted_plot:q} \
          --output {output:q}
        """
