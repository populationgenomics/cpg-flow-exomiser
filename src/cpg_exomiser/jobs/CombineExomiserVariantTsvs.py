from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_combine_exomiser_variant_tsvs_job(
    dataset: targets.Dataset,
    input_files: dict[str, Path],
    output: str,
    job_attrs: dict[str, str],
) -> 'BashJob':
    """consolidate all the exomiser gene tsvs into one file."""
    batch = hail_batch.get_batch()

    localised_files = []
    for family, file in input_files.items():
        if family.endswith('_variants'):
            localised_files.append(batch.read_input(str(file)))

    job = batch.new_bash_job(
        f'CombineExomiserVariantTsvs: {dataset.name}',
        attributes=job_attrs | {'tool': 'python'},
    )
    job.storage('10Gi')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.declare_resource_group(
        output={
            'json': '{root}.json',
            'ht.tar': '{root}.ht.tar',
        }
    )
    job.command(f"""
        python -m cpg_exomiser.scripts.combine_exomiser_variant_tsvs \\
        --input {' '.join(localised_files)} \\
        --output {job.output}
    """)
    # tar the Hail Table so we can remove as a single file
    job.command(f'tar --remove-files -cf {job.output}.ht.tar {job.output}.ht')
    batch.write_output(job.output, output.removesuffix('.json'))

    return job
