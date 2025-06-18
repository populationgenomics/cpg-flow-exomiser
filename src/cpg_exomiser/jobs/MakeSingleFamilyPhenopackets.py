import json

from cpg_flow import targets, utils
from cpg_utils import Path

HPO_KEY: str = 'HPO Terms (present)'


def make_phenopackets(
    proband_dict: dict[str, list[targets.SequencingGroup]],
    out_paths: dict[str, Path],
) -> None:
    """
    find the minimal data to run an exomiser analysis
    n.b. these are not actually phenopackets - they are a simplified version

    Args:
        proband_dict (dict[str, list[SequencingGroup]]): proband CPG ID: [family members]
        out_paths (dict[str, Path]): corresponding output paths per family
    """

    for proband, members in proband_dict.items():
        # skip if already done
        if utils.exists(out_paths[proband]):
            continue

        # get all affected and unaffected
        affected = [sg for sg in members if str(sg.pedigree.phenotype) == '2']

        if not affected:
            raise ValueError(f'Family {proband} has no affected individuals, should not have reached here')

        # select the specific proband SG based on the ID (key in the dictionary)
        proband_sg = next(sg for sg in members if sg.id == proband)

        hpo_term_string = proband_sg.meta['phenotypes'].get(HPO_KEY, '')

        hpo_terms = hpo_term_string.split(',')

        # https://github.com/exomiser/Exomiser/blob/master/exomiser-cli/src/test/resources/pfeiffer-family.yml
        phenopacket: dict = {'family': proband, 'proband': proband_sg.id, 'hpoIds': hpo_terms}

        with out_paths[proband].open('w') as ppk_file:
            json.dump(phenopacket, ppk_file, indent=2)
