from typing import TYPE_CHECKING

from cpg_flow import utils
from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_an_exomiser_job(name: str) -> 'BashJob':
    """
    Sets up a single job with the given title, and returns it.
    """
    job = hail_batch.get_batch().new_bash_job(name=name)
    job.storage(config.config_retrieve(['workflow', 'exomiser_storage']))
    job.memory(config.config_retrieve(['workflow', 'exomiser_memory']))
    job.cpu(config.config_retrieve(['workflow', 'exomiser_cpu']))
    job.image(config.config_retrieve(['images', 'exomiser_14']))

    return job


def run_exomiser(content_dict: dict) -> list['BashJob']:
    """
    run jobs through Exomiser 14

    Args:
        content_dict ():

    Returns:

    """

    exomiser_version = config.config_retrieve(['images', 'exomiser_14']).split(':')[-1]
    exomiser_dir = f'/exomiser/exomiser-cli-{exomiser_version}'

    chunks_per_vm: int = config.config_retrieve(['workflow', 'exomiser_chunk_size'], 8)
    chunks_in_parallel: int = config.config_retrieve(['workflow', 'exomiser_parallel_chunks'], 4)

    # localise the compressed exomiser references
    inputs = hail_batch.get_batch().read_input_group(
        core=config.reference_path('exomiser_2402_core'),
        pheno=config.reference_path('exomiser_2402_pheno'),
    )

    # now chunk the jobs - load resources, then run a bunch of families
    probands = sorted(content_dict.keys())
    all_jobs = []
    for chunk_number, proband_chunk in enumerate(utils.chunks(probands, chunks_per_vm)):
        # see https://exomiser.readthedocs.io/en/latest/installation.html#linux-install
        job = make_an_exomiser_job(f'Run Exomiser for chunk {chunk_number}')
        all_jobs.append(job)

        # unpack references, see linux-install link above
        job.command(rf'unzip {inputs}/\* -d "{exomiser_dir}/data"')

        job.command(f'echo "This job contains families {" ".join(proband_chunk)}"')

        # number of chunks should match cpu, accessible in config
        # these will all run simultaneously using backgrounded tasks and a wait
        for parallel_chunk in utils.chunks(proband_chunk, chunk_size=chunks_in_parallel):
            for proband in parallel_chunk:
                # read in VCF & index
                vcf = hail_batch.get_batch().read_input_group(
                    **{
                        f'{proband}_vcf': content_dict[proband]['vcf'],
                        f'{proband}_vcf_index': f'{content_dict[proband]["vcf"]}.tbi',
                    },
                )[f'{proband}_vcf']

                # read in ped & phenotype JSON
                ped = hail_batch.get_batch().read_input(content_dict[proband]['ped'])
                ppk = hail_batch.get_batch().read_input(content_dict[proband]['pheno'])

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

                hail_batch.get_batch().write_output(
                    job[proband],
                    str(content_dict[proband]['output']).removesuffix('.tsv'),
                )
    return all_jobs
