import pandas as pd

from cpg_utils import Path
from cpg_flow import utils, targets


def extract_mini_ped_files(proband_dict: dict[str, list[targets.SequencingGroup]], out_paths: dict[str, Path]):
    """
    write the mini-ped for each family

    Args:
        proband_dict (dict[str, list[SequencingGroup]]): proband SG ID to list of SG IDs
        out_paths (Path): temp dir to write the mini-peds to
    """

    # query for SG entities, group by family
    for proband_id, members in proband_dict.items():
        ped_path = out_paths[proband_id]
        # don't recreate if it exists
        if not utils.exists(ped_path):
            # make the pedigree for this family
            ped_df = pd.DataFrame([sg.pedigree.get_ped_dict() for sg in members])
            with ped_path.open('w') as ped_file:
                ped_df.to_csv(ped_file, sep='\t', index=False, header=False)
