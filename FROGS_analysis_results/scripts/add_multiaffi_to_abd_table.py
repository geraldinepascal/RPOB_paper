#!/usr/bin/env python3

"""
Description

:Example:
python template.py -v
"""

# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'support.bioinfo.genotoul@inra.fr'
__status__ = 'dev'


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import sys
from frogs_analysis_fct import process_frogs_affiliation, get_corresponding_species_mock, clean_mock_sp_relation


def load_mock_taxonomies(mock_taxonomies_file):
    
    with open(mock_taxonomies_file) as fl:
        return {l.rstrip().split("\t")[1]:l.split("\t")[0] for l in fl if l.rstrip()}
    

def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('--abundance_table', required=True, help='Affiliation abundance table, output of Frogs affiliations_stat.py script followed by biom_to_tsv.py. ')

    parser.add_argument('--multiaffi_table', required=True, help='Multi Affiliation table, output of Frogs affiliations_stat.py script followed by biom_to_tsv.py. ')

    parser.add_argument('--region', required=True, help='Name of the region analysed (ie 16S, 23S or 16Sv3v4).')

    parser.add_argument('--affi_db_name', required=True, help='Name of the taxonomic affiliation database. ie 16S silva. custom 16S-23S')
    
    parser.add_argument('--min_identity',  default=98, type=float, help='Identity threshold to consider an affilition as weak.')
    
    parser.add_argument('--min_coverage', default=99, type=float, help='Coverage threshold to consider an affilition as weak.')

    parser.add_argument('--taxonomic_ranks', default='Domain Phylum Class Order Family Genus Species Species', help='Taxonomic ranks of the affiliation')

    parser.add_argument('--mock_taxonomies', default=None, help='Specify this argument if the sample is from mock communities. '
                        'It allows clusters to be linked to their corresponding mock species. '
                        'Provide a tabulated file with two columns: the first column containing the name of the mock species, '
                        'and the second column containing the taxonomy of this species based on the affiliation database used. '
                        'The taxonomy should be represented from Phylum to Species, separated by semicolons.')

    parser.add_argument('-o', '--output', default='affi_tables_merged.tsv', help='output table name')

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()
    return args


def main():

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    affi_abundance_fl = args.abundance_table
    multihit_fl = args.multiaffi_table
    taxonomic_ranks = args.taxonomic_ranks.split(' ')
    
    min_identity = args.min_identity
    min_coverage = args.min_coverage

    region  = args.region
    affi_db_name = args.affi_db_name
    mock_taxonomies_file = args.mock_taxonomies


    analysis_df =  process_frogs_affiliation(affi_abundance_fl, multihit_fl, taxonomic_ranks,
                                                                min_ids = [min_identity], min_covs = [min_coverage])

    df_valid_affi = analysis_df.loc[analysis_df[f'id>{min_identity}_cov>{min_coverage}']]
    logging.info(f"{len(df_valid_affi)}/{len(analysis_df)} clusters have a valid affiliation with identity > {min_identity} and coverage > {min_coverage}.")

    sum_seq_valid = df_valid_affi['observation_sum'].sum()
    abd_of_valid_affi = 100 * sum_seq_valid/analysis_df['observation_sum'].sum()

    logging.info(f"They represent a total of {sum_seq_valid} sequences and a relative abundance of {abd_of_valid_affi:.3f}%")


    analysis_df['valid_affiliation'] = analysis_df[f'id>{min_identity}_cov>{min_coverage}']

    analysis_df['region'] = region
    analysis_df['db'] = affi_db_name

    if mock_taxonomies_file:
        logging.info(f'Linking cluster to their corresponding mock species using {mock_taxonomies_file}')
        taxonomies2mock_species = load_mock_taxonomies(mock_taxonomies_file)

        analysis_df["sp_mock_taxonomy"] = analysis_df[["blast_taxonomy", "mutliaffiliation"]].apply(
                get_corresponding_species_mock, args=(taxonomies2mock_species,), axis=1 )

        analysis_df["sp_mock_taxonomy"] = analysis_df["sp_mock_taxonomy"].apply(clean_mock_sp_relation)
        analysis_df["mock_species"] = analysis_df["sp_mock_taxonomy"].apply(lambda x: x if x not in taxonomies2mock_species else taxonomies2mock_species[x] )


    logging.info(f'Writting merged table in {args.output}')
    analysis_df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
