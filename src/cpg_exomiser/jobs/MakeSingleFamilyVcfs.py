from typing import TYPE_CHECKING

from cpg_flow import targets, utils
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def family_vcf_from_gvcf(family_members: list[targets.SequencingGroup], out_path: str) -> 'BashJob':
    """
    Does a quick blast of gVCF -> VCF
    Strips out ref-only sites, and splits Alt, Non_Ref into just Alt
    If there are multiple members, this finishes with a merge

    Args:
        family_members (list[SequencingGroup]): member(s) to convert to VCF
        out_path (str): where to write the VCF (and implicitly, the index)

    Returns:
        the job, for dependency setting
    """

    family_ids = [sg.id for sg in family_members]
    job = hail_batch.get_batch().new_bash_job(f'Generate VCF {out_path} from {family_ids}')
    job.image(config.config_retrieve(['images', 'bcftools']))

    # for genomes, we need a bit more storage, exomes the 5GB default is fine
    if config.config_retrieve(['workflow', 'sequencing_type']) == 'genome':
        job.storage('10Gi')

    # read input
    family_vcfs = [hail_batch.get_batch().read_input(sg.gvcf) for sg in family_members]

    # declare a resource group
    job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # view -m 3 to strip out ref-only sites
    # norm -m -any to split Alt, Non_Ref into just Alt
    # grep -v NON_REF to remove the NON_REF sites
    # bgzip -c to write to a compressed file
    if len(family_vcfs) == 1:
        gvcf_input = family_vcfs[0]
        job.command(
            f"""
            bcftools view -m3 {gvcf_input} | \\
            bcftools norm -m -any | \\
            grep -v NON_REF | \\
            bgzip -c  > {job.output['vcf.bgz']}
            """
        )
        job.command(f'tabix {job.output["vcf.bgz"]}')
        hail_batch.get_batch().write_output(job.output, out_path.removesuffix('.vcf.bgz'))
        return job

    # if there are multiple members, convert and merge them
    job.cpu(4)
    paths = []
    for index, gvcf_input in enumerate(family_vcfs):
        job.command(
            f'bcftools view -m3 {gvcf_input} | bcftools norm -m -any | grep -v NON_REF | bgzip -c  > {index}.vcf.bgz',
        )
        job.command(f'tabix {index}.vcf.bgz')
        paths.append(f'{index}.vcf.bgz')
    # merge the VCFs
    # -m all to merge all sites
    # -0 to replace missing with WT (potentially inaccurate, but predictable parsing in exomiser)
    # --threads 4 to use 4 threads
    # -Oz to write a compressed VCF
    job.command(f'bcftools merge {" ".join(paths)} -Oz -o {job.output["vcf.bgz"]} --threads 4 -m all -0')
    job.command(f'tabix {job.output["vcf.bgz"]}')
    hail_batch.get_batch().write_output(job.output, out_path.removesuffix('.vcf.bgz'))
    return job


def create_gvcf_to_vcf_jobs(
    proband_dict: dict[str, list[targets.SequencingGroup]],
    previous_completions: set[str],
    out_paths: dict[str, Path],
) -> list['BashJob']:
    """
    Create Joint VCFs for families of SG IDs

    Args:
        proband_dict (): dict of proband ID to list of SG IDs
        previous_completions (set[str]): set of analyses we've already completed
        out_paths (): dict of family ID to output path
    Returns:
        list of Jobs
    """

    jobs = []

    # take each family
    for proband, members in proband_dict.items():
        # skip if already done
        if utils.exists(out_paths[proband]) or proband in previous_completions:
            continue

        jobs.append(family_vcf_from_gvcf(members, out_paths[proband]))
    return jobs
