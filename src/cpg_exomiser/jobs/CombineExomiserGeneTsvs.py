from typing import TYPE_CHECKING

from cpg_utils import Path, hail_batch, config
from cpg_flow import targets

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_combine_exomiser_gene_tsvs_job(
    dataset: targets.Dataset,
    input_files: dict[str, Path],
    project_id: str,
    output: Path,
    job_attrs: dict[str, str],
) -> 'BashJob':
    """consolidate all the exomiser gene tsvs into one file."""
    batch = hail_batch.get_batch()

    localised_files = []
    for family, file in input_files.items():
        if family.endswith('_variants'):
            continue

        localised_files.append(batch.read_input(str(file)))

    job = batch.new_bash_job(
        f'CombineExomiserGeneTsvs: {dataset.name}, {project_id}',
        attributes=job_attrs | {'tool': 'python'},
    )
    job.storage('10Gi')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command(f"""
        python -m cpg_exomiser.scripts.combine_exomiser_gene_tsvs \\
            --project {project_id} \\
            --input {' '.join(localised_files)} \\
            --output {job.output}
    """)
    batch.write_output(job.output, output)

    return job
