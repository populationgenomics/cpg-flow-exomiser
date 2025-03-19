"""
As a dataset Stage

- Iterate over the families
- For each family export a VCF & subselected PED
- split the families into chunks
- for each chunk copy in the resources & set up a config file
- run exomiser, multiple families per setup
"""

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

# import job logic per-stage
from cpg_exomiser.jobs.CreateFamilyVcfs import create_gvcf_to_vcf_jobs
from cpg_exomiser.jobs.MakePhenopackets import make_phenopackets
from cpg_exomiser.jobs.MakePedExtracts import extract_mini_ped_files
from cpg_exomiser.jobs.RunExomiser import run_exomiser

from cpg_exomiser.utils import find_seqr_projects, find_probands, find_previous_analyses

from cpg_flow.utils import get_logger
from cpg_flow.stage import DatasetStage, SequencingGroupStage, stage
from cpg_flow.workflow import get_workflow

if TYPE_CHECKING:
    from cpg_utils import Path
    from cpg_flow.targets.dataset import Dataset
    from cpg_flow.targets.sequencing_group import SequencingGroup
    from cpg_flow.stage import StageInput, StageOutput


EXOMISER_ANALYSIS_TYPE: str = 'exomiser'


@stage
class MakeSingleFamilyVcfs(DatasetStage):
    """
    it's a gVCF combiner, densification, Mt -> VCF
    we do this densification per-proband, not per-family... If there are multiple probands in the family we'll
    produce multiple versions of the same VCF. It's simpler to just do it per-proband rather than implement some
    sharing logic here, which would be more complex for a marginal benefit
    """

    def expected_outputs(self, dataset: 'Dataset') -> dict[str, 'Path']:
        """
        this now writes to a temporary bucket, we don't need these VCFs again in the future
        they cost more to keep than to regenerate
        """
        prefix = dataset.tmp_prefix() / 'exomiser_inputs'
        return {proband: prefix / f'{proband}.vcf.bgz' for proband in find_probands(dataset)}

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(dataset)
        jobs = create_gvcf_to_vcf_jobs(
            proband_dict=find_probands(dataset),
            previous_completions=find_previous_analyses(dataset.name),
            out_paths=outputs,
        )
        return self.make_outputs(dataset, outputs, jobs=jobs)


@stage
class MakeSingleFamilyPhenopackets(DatasetStage):
    """
    for each relevant family, make some Phenopackets
    """

    def expected_outputs(self, dataset: 'Dataset'):
        dataset_prefix = dataset.analysis_prefix() / 'exomiser_inputs'
        return {proband: dataset_prefix / f'{proband}_phenopacket.json' for proband in find_probands(dataset)}

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        """
        this actually doesn't run as Jobs, but as a function
        """

        expected_out = self.expected_outputs(dataset)
        make_phenopackets(
            proband_dict=find_probands(dataset),
            out_paths=expected_out,
        )
        return self.make_outputs(dataset, data=expected_out)


@stage
class MakeSingleFamilyPedFiles(DatasetStage):
    """
    from the dataset MT, we make a PED per-proband
    """

    def expected_outputs(self, dataset: 'Dataset'):
        dataset_prefix = dataset.analysis_prefix() / 'exomiser_inputs'
        return {proband: dataset_prefix / f'{proband}.ped' for proband in find_probands(dataset)}

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        """
        this actually doesn't run as Jobs, but as a function...
        these outputs will be written by the time the pipeline starts
        """

        get_logger().debug(f'not using inputs: {inputs}')

        expected_out = self.expected_outputs(dataset)
        extract_mini_ped_files(
            proband_dict=find_probands(dataset),
            out_paths=expected_out,
        )
        return self.make_outputs(dataset, data=expected_out)


@stage(required_stages=[MakeSingleFamilyVcfs, MakeSingleFamilyPedFiles, MakeSingleFamilyPhenopackets])
class RunExomiser(DatasetStage):
    """
    Run exomiser on the family VCFs

    Only run this for families with at least one affected individual

    frustratingly due to the output format we can't really register these in metamist
    we could if it was a SequencingGroupStage, but it's not
    """

    def expected_outputs(self, dataset: 'Dataset'):
        """
        dict of outputs for this dataset, keyed on family ID
        """
        proband_dict = find_probands(dataset)

        dataset_prefix = dataset.analysis_prefix() / 'exomiser_results'

        # only the TSVs are required, but we need the gene and variant level TSVs
        # populate gene-level results
        return_dict = {proband: dataset_prefix / f'{proband}.tsv' for proband in proband_dict}
        # add more keys pointing to the variant-level TSVs
        return_dict.update(
            {f'{family}_variants': dataset_prefix / f'{family}.variants.tsv' for family in proband_dict},
        )

        return return_dict

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':

        output_dict = self.expected_outputs(dataset)

        vcf_inputs = inputs.as_dict(target=dataset, stage=MakeSingleFamilyVcfs)
        ped_inputs = inputs.as_dict(target=dataset, stage=MakeSingleFamilyPedFiles)
        ppk_files = inputs.as_dict(target=dataset, stage=MakeSingleFamilyPhenopackets)

        # combining all these stage outputs into one object
        single_dict = {
            family: {
                'output': output_dict[family],
                'vcf': vcf_inputs[family],
                'ped': ped_inputs[family],
                'pheno': ppk_files[family],
            }
            for family in output_dict
            if '_variants' not in family
        }

        jobs = run_exomiser(single_dict)
        return self.make_outputs(dataset, data=output_dict, jobs=jobs)


@stage(
    analysis_keys=['gene_level', 'variant_level'],
    required_stages=[RunExomiser],
    analysis_type=EXOMISER_ANALYSIS_TYPE,
)
class RegisterSingleSampleExomiserResults(SequencingGroupStage):
    """
    this is a tricky little fella'
    The previous Stage, RunExomiser, runs per-Dataset, but writes output per-Proband
    The reason is that Exomiser instances are expensive (time and compute) to start up, so we start a few instances,
    and really crush high numbers of samples through each VM. That's really cost effective, but it means that the
    expected_outputs dict is populated at runtime, so the individual output files can't be entered into analysis_keys
    That means that we can't write the per-proband result entries to metamist, from the previous stage, so we do it here
    """

    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> dict[str, 'Path']:
        """
        These output paths are identical to the previous Stage, but will be registered in metamist
        """
        # check if we're interested in this SG ID
        proband_dict = find_probands(sequencing_group.dataset)
        if sequencing_group.id not in proband_dict:
            return {}
        output_prefix = sequencing_group.dataset.analysis_prefix() / 'exomiser_results'
        return {
            'gene_level': output_prefix / f'{sequencing_group.id}.tsv',
            'variant_level': output_prefix / f'{sequencing_group.id}.variants.tsv',
        }

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':
        """
        get the dictionary of all probands to generate an analysis for
        if an entry in that dict exists for this SG ID, create a ghost job to get it registered in metamist
        metamist/status_reporter registrations won't run unless there's a job in the Stage
        """

        get_logger().debug(f'not using inputs: {inputs}')

        if not (outputs := self.expected_outputs(sequencing_group)):
            return self.make_outputs(sequencing_group, jobs=None)

        # if we're interested in this SG ID, create a spooky ghost job - exists, but no actions
        ghost_job = get_batch().new_job(f'Register {sequencing_group.id} Exomiser files in metamist')
        ghost_job.image(config_retrieve(['workflow', 'driver_image']))
        ghost_job.command(f'echo "I am the ghost of {sequencing_group.id}, oOOooOooOoOo"')
        return self.make_outputs(sequencing_group, jobs=ghost_job, data=outputs)


@stage(required_stages=[RunExomiser], analysis_type=EXOMISER_ANALYSIS_TYPE)
class CombineExomiserGeneTsvs(DatasetStage):
    """
    Parse the Exomiser results into a TSV for Seqr
    """

    def expected_outputs(self, dataset: 'Dataset') -> 'Path':
        return dataset.analysis_prefix() / get_workflow().output_version / 'exomiser_results.tsv'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':

        output = self.expected_outputs(dataset)

        # is there a seqr project?
        projects = find_seqr_projects()

        if dataset.name not in projects:
            get_logger(__file__).info(f'No Seqr project found for {dataset.name}, skipping')
            return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[], skipped=True)

        results = inputs.as_dict(target=dataset, stage=RunExomiser)

        # just collect the per-family gene-level TSVs
        local_files: list = []
        for family, file in results.items():
            if family.endswith('_variants'):
                continue

            local_files.append(get_batch().read_input(file))

        # assign the project name to a variable
        project = projects[dataset.name]

        job = get_batch().new_bash_job(f'CombineExomiserGeneTsvs: {dataset.name}, {project}')
        job.storage('10Gi')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'combine_exomiser_genes --project {project} --input {" ".join(local_files)} --output {job.output}')
        get_batch().write_output(job.output, str(output))

        return self.make_outputs(dataset, data=output, jobs=job)


@stage(required_stages=[RunExomiser], analysis_type=EXOMISER_ANALYSIS_TYPE, analysis_keys=['json', 'ht'])
class CombineExomiserVariantTsvs(DatasetStage):
    """
    Parse the Exomiser variant-level results into a JSON file and a Hail Table
    """

    def expected_outputs(self, dataset: 'Dataset'):

        prefix = dataset.analysis_prefix() / get_workflow().output_version

        return {
            'json': prefix / 'exomiser_variant_results.json',
            'ht': prefix / 'exomiser_variant_results.ht.tar',
        }

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':

        if not find_probands(dataset):
            get_logger().info('No families found, skipping exomiser')
            return self.make_outputs(dataset, jobs=None, skipped=True)

        outputs = self.expected_outputs(dataset)

        # collect the per-family variant-level TSVs
        family_files: list = []
        for family, file in inputs.as_dict(target=dataset, stage=RunExomiser).items():
            if family.endswith('_variants'):
                family_files.append(get_batch().read_input(file))

        job = get_batch().new_job(f'CombineExomiserVariantTsvs: {dataset.name}')
        job.declare_resource_group(
            output={
                'json': '{root}.json',
                'ht.tar': '{root}.ht.tar'
            }
        )

        job.image(config_retrieve(['workflow', 'driver_image']))

        # generate the outputs
        job.command(f'combine_exomiser_variants --input {" ".join(family_files)} --output {job.output}')

        # tar the Hail Table so we can remove as a single file
        job.command(f'tar --remove-files -cf {job.output}.ht.tar {job.output}.ht')

        # extract both output files together
        get_batch().write_output(job.output, str(outputs['json']).removesuffix('.json'))

        return self.make_outputs(dataset, data=outputs, jobs=job)
