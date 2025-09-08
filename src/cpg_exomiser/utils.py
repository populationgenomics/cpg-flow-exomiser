"""
utilities methods for interacting with metamist and detecting pedigrees in an analysis
"""

from functools import cache

from cpg_flow import targets
from cpg_utils import config, to_path
from loguru import logger
from metamist.apis import ProjectApi
from metamist.graphql import gql, query

# GraphQl query for all analysis entries in this Dataset of type: exomiser
# this is used to find the exomiser results for the current dataset
ANALYSIS_QUERY = gql(
    """
query MyQuery($dataset: String!, $analysis_type: String!) {
  project(name: $dataset) {
    analyses(type: {eq: $analysis_type}) {
      meta
      outputs
    }
  }
}
""",
)
HPO_KEY: str = 'HPO Terms (present)'
EXOMISER_ANALYSIS_TYPE: str = 'exomiser'
EXOMISER_VERSION = config.config_retrieve(['workflow', 'exomiser_version'])
EXOMISER_DATA_VERSION = config.config_retrieve(['workflow', 'exomiser_data_version'])


@cache
def find_seqr_projects() -> dict[str, str]:
    """
    query for all the seqr projects
    map datasets to corresponding seqr names
    """

    project_api = ProjectApi()
    seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
    seq_type_key = f'seqr-project-{seq_type}'
    return_dict: dict[str, str] = {}

    for element in project_api.get_seqr_projects():
        if element['meta'].get('is_seqr', False) and (meta_proj := element['meta'].get(seq_type_key, '')):
            return_dict[element['dataset']] = meta_proj

    return return_dict


@cache
def find_previous_analyses(dataset: str) -> set[str]:
    """
    Find all the analysis entries in this dataset of type: exomiser

    Args:
        dataset (str): the name of the dataset to find analysis entries for

    Returns:
        set[str]: analysis paths in metamist for this Dataset
    """

    metamist_results = query(ANALYSIS_QUERY, variables={'dataset': dataset, 'analysis_type': EXOMISER_ANALYSIS_TYPE})
    completed_runs: set[str] = set()
    for analysis in metamist_results['project']['analyses']:
        # uses the new metamist outputs formatted block
        if outputs := analysis['outputs']:
            # sub-select analysis entries to the ones with this matching Exomiser version (software and data)
            if (analysis['meta'].get('exomiser_version', None) != EXOMISER_VERSION) or (
                analysis['meta'].get('exomiser_data_version', None) != EXOMISER_DATA_VERSION
            ):
                continue
            if isinstance(outputs, dict) and (nameroot := outputs.get('nameroot')):
                completed_runs.add(nameroot)
            elif isinstance(outputs, str):
                completed_runs.add(to_path(outputs).stem)
    return completed_runs


@cache
def find_probands(dataset: targets.Dataset) -> dict[str, list[targets.SequencingGroup]]:
    """
    Find all the families in the project
    group on family ID and check for affected individuals & HPO terms
    Export a dictionary with every affected member in the family separately

    Do some management there to make sure we only retain affecteds at the leaf-end of the pedigree
    - once we get all affected in the pedigree, scan them to see if any have parents
    - if at least one has a parent, we don't include any affected without parents in the analysis
    - they can be affected parents of probands, but shouldn't be the centre of an exomiser analysis?

    Args:
        dataset ():

    Returns:
        dict[str, list[SequencingGroup]]: a dictionary of affected individuals to all family members
    """

    dict_by_family: dict[str, list[targets.SequencingGroup]] = {}
    for sg in dataset.get_sequencing_groups():
        # skip over members without a gVCF - unusable in this analysis
        if sg.gvcf is None:
            continue
        family_id = str(sg.pedigree.fam_id)
        dict_by_family.setdefault(family_id, []).append(sg)

    dict_of_affecteds: dict[str, list[targets.SequencingGroup]] = {}

    # now remove families with no affected individuals
    for family, members in dict_by_family.items():
        # check for at least one retained member
        affected = [sg for sg in members if str(sg.pedigree.phenotype) == '2']

        # remove families with no affected members
        if not affected:
            logger.info(f'Family {family} has no affected individuals, skipping')
            continue

        # check that the affected members have HPO terms - required for exomiser
        if any(sg.meta['phenotypes'].get(HPO_KEY, '') == '' for sg in affected):
            logger.info(f'Family {family} has affected individuals with no HPO terms, skipping')
            continue

        # check to see if any affected member has parents
        # (if so, we don't treat any affected-without-parents as relevant for Exomiser; cannot be a proband)
        # mom and dad here are SGs, not just empty pedigree members
        affected_with_parents: bool = False
        for member in affected:
            if member.pedigree.mom or member.pedigree.dad:
                affected_with_parents = True
                break

        for member in affected:
            if not (member.pedigree.mom or member.pedigree.dad) and affected_with_parents:
                # if there are affected members without parents, we don't want to include them
                continue

            dict_of_affecteds[member.id] = members

    return dict_of_affecteds
