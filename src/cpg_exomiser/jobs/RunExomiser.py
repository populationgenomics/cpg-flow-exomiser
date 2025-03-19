from typing import TYPE_CHECKING
from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_utils.hail_batch import get_batch

from cpg_flow.utils import chunks


if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_exomiser(content_dict: dict) -> list['Job']:
    """
    run jobs through Exomiser 14

    Args:
        content_dict ():

    Returns:

    """

    exomiser_version = image_path('exomiser_14').split(':')[-1]
    exomiser_dir = f'/exomiser/exomiser-cli-{exomiser_version}'

    # localise the compressed exomiser references
    inputs = get_batch().read_input_group(
        core=reference_path('exomiser_2402_core'),
        pheno=reference_path('exomiser_2402_pheno'),
    )

    # now chunk the jobs - load resources, then run a bunch of families
    probands = sorted(content_dict.keys())
    all_jobs = []
    for chunk_number, proband_chunk in enumerate(
        chunks(probands, config_retrieve(['workflow', 'exomiser_chunk_size'], 8)),
    ):
        # see https://exomiser.readthedocs.io/en/latest/installation.html#linux-install
        job = get_batch().new_bash_job(f'Run Exomiser for chunk {chunk_number}')
        all_jobs.append(job)
        job.storage(config_retrieve(['workflow', 'exomiser_storage'], '100Gi'))
        job.memory(config_retrieve(['workflow', 'exomiser_memory'], '60Gi'))
        job.cpu(config_retrieve(['workflow', 'exomiser_cpu'], 4))
        job.image(image_path('exomiser_14'))

        # unpack references, see linux-install link above
        job.command(rf'unzip {inputs}/\* -d "{exomiser_dir}/data"')

        job.command(f'echo "This job contains families {" ".join(proband_chunk)}"')

        # number of chunks should match cpu, accessible in config
        # these will all run simultaneously using backgrounded tasks and a wait
        for parallel_chunk in chunks(
            proband_chunk,
            chunk_size=config_retrieve(['workflow', 'exomiser_parallel_chunks'], 4),
        ):
            for proband in parallel_chunk:
                # read in VCF & index
                vcf = get_batch().read_input_group(
                    **{
                        f'{proband}_vcf': str(content_dict[proband]['vcf']),
                        f'{proband}_vcf_index': f'{content_dict[proband]["vcf"]}.tbi',
                    },
                )[f'{proband}_vcf']

                # read in ped & phenotype JSON
                ped = get_batch().read_input(str(content_dict[proband]['ped']))
                ppk = get_batch().read_input(str(content_dict[proband]['pheno']))

                # # this was really satisfying syntax to work out
                job.declare_resource_group(
                    **{
                        proband: {
                            'json': '{root}.json',
                            'tsv': '{root}.tsv',
                            'variants.tsv': '{root}.variants.tsv',
                            'yaml': '{root}.yaml',
                        },
                    },
                )

                # generate a config file based on the batch tmp locations
                job.command(f'python3 {exomiser_dir}/config_shuffle.py {ppk} {job[proband]["yaml"]} {ped} {vcf} ')

                # now run it, as a backgrounded process
                job.command(
                    f'java -Xmx10g -Xms4g -jar {exomiser_dir}/exomiser-cli-{exomiser_version}.jar '
                    f'--analysis {job[proband]["yaml"]} --ped {ped} '
                    f'--spring.config.location={exomiser_dir}/application.properties &',
                )

            # wait for backgrounded processes to finish, show current state
            job.command('wait && ls results')

            # move the results, then copy out
            for proband in parallel_chunk:
                job.command(f'mv results/{proband}.json {job[proband]["json"]}')
                job.command(f'mv results/{proband}.genes.tsv {job[proband]["tsv"]}')
                job.command(f'mv results/{proband}.variants.tsv {job[proband]["variants.tsv"]}')

                get_batch().write_output(
                    job[proband],
                    str(content_dict[proband]['output']).removesuffix('.tsv'),
                )
    return all_jobs
