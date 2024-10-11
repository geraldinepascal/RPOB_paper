# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'


import pandas as pd
import logging
from itertools import takewhile

def clean_mock_sp_relation(mock_sp_relation):
    if len(mock_sp_relation) == 1:
        return mock_sp_relation.pop()
    elif len(mock_sp_relation) == 0:
        return "Other"
    elif len(mock_sp_relation) > 1:
        logging.warning(f"A cluster have more than one mock species in its affiliation: \n {mock_sp_relation}")
        return 'multiple mock species'
    
    
def get_corresponding_species_mock(tax_and_multiaffi, mock_taxonomies):
    taxonomy_str = tax_and_multiaffi["blast_taxonomy"]
                
    mutliaffiliation_str = tax_and_multiaffi["mutliaffiliation"]
    
    
    sp_hits = set()
    if taxonomy_str.endswith("Multi-affiliation"):
        # mutliaffiliation_str = mutliaffiliation_str.replace(
        #    alternative_fermentum_taxon, "Lactobacillus" )
        
        for affi in mutliaffiliation_str.split('|'):
            
            sp_hits |= get_corresponding_species_mock({'blast_taxonomy':affi, 'mutliaffiliation':None}, mock_taxonomies)
          

    else:
        for mock_tax in mock_taxonomies:
            if taxonomy_str.startswith(mock_tax):
                sp_hits.add(mock_tax)
    
    return sp_hits


def consider_only_selected_samples(df, all_sample, samples_to_keep):

    samples_to_remove = all_sample - samples_to_keep

    # df = df.drop(columns=samples_to_remove)
    df["observation_sum"] = df[samples_to_keep].sum(axis=1)
    df = df.loc[df["observation_sum"] > 0]
    df["abundance"] = 100 * df["observation_sum"] / df["observation_sum"].sum()

    return df


def process_frogs_affiliation(
    affi_abundance_file, multiaff_file, ranks, min_ids=[99], min_covs=[40, 99]
):

    df = pd.read_csv(affi_abundance_file, sep="\t")
    df_multiaff = pd.read_csv(multiaff_file, sep="\t")

    df_multiaff = df_multiaff.set_index("#observation_name", drop=False)

    add_multi_affi_to_df(df, df_multiaff)

    improve_affi_with_multiaffi(df, df_multiaff)

    # Add coverage and identity value to the cluster with a multiaffi by taking the value in multiaffi table
    cols_impacted_by_multihit = [
        "observation_name",
        "blast_perc_identity",
        "blast_perc_query_coverage",
    ]
    df[cols_impacted_by_multihit] = df[cols_impacted_by_multihit].apply(
        get_best_covid_from_multihit, args=(df_multiaff,), axis=1
    )

    # Add rank of the affiliation
    df["taxon_affi"] = df["blast_taxonomy"].apply(lambda x: get_taxon_affi(x, ranks))
    df["rank_affi"] = df["blast_taxonomy"].apply(lambda x: get_rank_affi(x, ranks))

    df["blast_taxonomy_cleaned"] = df["blast_taxonomy"].apply(
        lambda x: clean_taxonomy(x, ranks)
    )

    # df['blast_taxonomy'] = df['blast_taxonomy'].apply(rm_strain_from_lineage)

    df[ranks] = df["blast_taxonomy"].str.split(pat=";", expand=True)

    df = df.astype({"blast_perc_identity": float, "blast_perc_query_coverage": float})

    for min_id in min_ids:
        for min_cov in min_covs:

            df[f"coverage_>=_{min_cov}"] = df["blast_perc_query_coverage"] >= min_cov
            df[f"identity_>=_{min_id}"] = df["blast_perc_identity"] >= min_id
            df[f"id>{min_id}_cov>{min_cov}"] = (df["blast_perc_identity"] >= min_id) & (
                df["blast_perc_query_coverage"] >= min_cov
            )

    df["abundance"] = 100 * df["observation_sum"] / df["observation_sum"].sum()

    return df


# def get_taxon_affi(taxonomy_str):
#     taxonomy = [t for t in taxonomy_str.split(';') if t != 'Multi-affiliation']
#     if len(taxonomy) == 8:
#         taxonomy.pop() # ignore strain name
#     taxon_affi = taxonomy.pop()
#     if taxon_affi.count('_') > 1: # remove strain name in silva taxo
#         taxon_affi = '_'.join(taxon_affi.split('_')[:1])

#     return taxon_affi


def get_taxon_affi(taxonomy_str, ranks):

    rank, taxon = get_rank_and_taxon_affi(taxonomy_str, ranks)
    return taxon


def get_rank_affi(taxonomy_str, ranks):
    rank, taxon = get_rank_and_taxon_affi(taxonomy_str, ranks)
    return rank


def get_rank_and_taxon_affi(taxonomy_str, ranks):
    print(taxonomy_str)
    if taxonomy_str == "no data":
        return ("unknown", "unknown")

    # taxonomy = [t for t in taxonomy_str.split(';') if t != 'Multi-affiliation']
    # print(taxonomy)

    rank = "unknown"
    taxon = "unknown"
    for rank, taxon in zip(ranks[::-1], taxonomy_str.split(";")[::-1]):
        if (
            taxon.startswith("unknown")
            or "metagenome" in taxon
            or taxon == "Multi-affiliation"
        ):
            continue
        else:
            break

    # taxo_len = len(set(taxonomy))
    print(rank, taxon)
    return rank, taxon  # ranks[taxo_len-1]


def clean_taxonomy(taxonomy_str, ranks):
    
    if len(taxonomy_str.split(';')) > 7:
        # remove strain info when needed
        taxonomy_str = ';'.join(taxonomy_str.split(';')[:7])

    if len(ranks) > 7:
        ranks = ranks[:7]

    assert len(ranks) == len(taxonomy_str.split(';'))

    # affiliation
    uninformative_taxon_words = ["unknown", "metagenome"]
    clean_taxo = None
    # going trough
    taxo_split = taxonomy_str.split(";")
    for i, (rank, taxon) in enumerate(zip(ranks[::-1], taxo_split[::-1])):
        if any((w in taxon for w in uninformative_taxon_words)):
            logging.debug(
                f"Cleaning taxonomy: '{taxon}' is trimmed off because it is not informative enough from {taxonomy_str}"
            )
            continue
        else:
            clean_taxo = ";".join(taxo_split[: len(taxo_split) - i])
            break

    return clean_taxo


def get_sp_mock_relation(tax_and_multiaffi, mock_taxonomies):
    # print('='*50)

    taxonomy_str = tax_and_multiaffi["blast_taxonomy"]
    mutliaffiliation_str = tax_and_multiaffi["mutliaffiliation"]

    alternative_fermentum_taxon = "Limosilactobacillus"

    if alternative_fermentum_taxon in taxonomy_str:
        taxonomy_str = taxonomy_str.replace(
            alternative_fermentum_taxon, "Lactobacillus"
        )

    # taxonomy = [t for t in taxonomy_str.split(';') if t != 'Multi-affiliation']

    if taxonomy_str.endswith("Multi-affiliation"):
        mutliaffiliation_str = mutliaffiliation_str.replace(
            alternative_fermentum_taxon, "Lactobacillus"
        )
        sp_mock = get_sp_mock_relation_in_multiaffi(
            mutliaffiliation_str, mock_taxonomies
        )

        return sp_mock

    if taxonomy_str == "no data":
        return "unknown"

    # print(taxonomy_str)
    sp_hit = "other"
    for mock_tax in mock_taxonomies:
        # print('->', mock_tax)
        if taxonomy_str.startswith(mock_tax):
            # print(f'HITS')
            sp_hit = mock_tax.split(";")[-1]
            if sp_hit == "Escherichia_coli":
                break
        else:
            pass
            # print('      --> not in taxo')

    return sp_hit


def get_sp_mock_relation_in_multiaffi(mutliaffiliation_str, mock_taxonomies):

    mutliaffiliations = mutliaffiliation_str.split("|")
    for mock_tax in mock_taxonomies:

        for affi in mutliaffiliations:
            if len(affi.split(";")) == 8:
                affi = ";".join(affi.split(";")[:-1])  # remove strain
            # print(affi)
            if affi == mock_tax:
                # print('multi affi is found in mock_tax')
                return mock_tax.split(";")[-1]

    # print('multiaffi is not in mock')
    # [print('    ', a) for a in mutliaffiliations]
    return "other"


def improve_affiliation_df(df, mock_taxonomies, ranks):
    # df['is_in_mock'] = df['mutliaffiliation'].apply(is_multiaffi_in_mock, args=(mock_taxonomies,))

    df["taxon_affi"] = df["blast_taxonomy"].apply(lambda x: get_taxon_affi(x, ranks))
    df["rank_affi"] = df["blast_taxonomy"].apply(lambda x: get_rank_affi(x, ranks))

    df["sp_mock_related"] = df[["blast_taxonomy", "mutliaffiliation"]].apply(
        get_sp_mock_relation, args=(mock_taxonomies,), axis=1
    )

    df["species_mock"] = df["sp_mock_related"].str.replace("_", " ")


def add_multi_affi_to_df(df, df_multihit):
    for index, row in df.iterrows():
        if row["blast_subject"] == "multi-subject":

            cluster = row["observation_name"]

            multiaffi_list = list(df_multihit.loc[cluster, "blast_taxonomy"])

            df.loc[index, "mutliaffiliation"] = "|".join(multiaffi_list)




def get_best_covid_from_multihit(df_covid, df_multiaff):
    cluster = df_covid["observation_name"]
    blast_perc_identity = df_covid["blast_perc_identity"]
    blast_perc_query_coverage = df_covid["blast_perc_query_coverage"]

    if (
        blast_perc_identity == "multi-identity"
        or blast_perc_query_coverage == "multi-coverage"
    ):
        # select multiaffi of the cluster
        filt = df_multiaff["#observation_name"] == cluster
        cluster_multiaffi = df_multiaff.loc[filt]

        # sort multiaffi on id and cov
        cluster_multiaffi = cluster_multiaffi.sort_values(
            by=["blast_perc_identity", "blast_perc_query_coverage"], ascending=False
        )
        # set best hit in df
        df_covid["blast_perc_identity"] = cluster_multiaffi["blast_perc_identity"].iloc[
            0
        ]
        df_covid["blast_perc_query_coverage"] = cluster_multiaffi[
            "blast_perc_query_coverage"
        ].iloc[0]

    return df_covid

def improve_affi_with_multiaffi(df, df_multihit):
    
    df['blast_taxonomy_original'] = df['blast_taxonomy']

    for index, row in df.iterrows():
        # print(index, row['blast_taxonomy'])
        

        if row["blast_subject"] == "multi-subject":
            
            cluster = row["observation_name"]
            print('___', cluster)
            multiaffi_list = list(df_multihit.loc[cluster, "blast_taxonomy"])
            [print(t) for t in multiaffi_list]
            common_taxonomy = get_common_taxonomy(multiaffi_list, delete_strain_info=True)
            df.loc[index, 'blast_taxonomy'] = common_taxonomy
            print("COMMON",common_taxonomy)
            assert len(df.loc[index, 'blast_taxonomy'].split(';')) == 7, df.loc[index, 'blast_taxonomy']
            
def manage_strain_in_taxo(taxonomy, delete_strain_info=False):
    logging.debug(taxonomy)
    assert len(taxonomy) == 7
    if taxonomy[-1].count(' ') > 1 and 'metagenome' not in taxonomy[-1]:
        logging.debug('taxonomy species look like strain')
        sp = ' '.join(taxonomy[-1].split(' ')[:2])
        if delete_strain_info:
            taxonomy = taxonomy[:-1] + [sp]
        else:
            taxonomy = taxonomy[:-1] + [sp, taxonomy[-1]]
    logging.debug(taxonomy)
    return taxonomy

def get_common_taxonomy(taxonomies, delete_strain_info=False):
    for t in taxonomies:
        logging.debug(t)
        
    taxonomies_list = (manage_strain_in_taxo(t.split(";"), delete_strain_info) for t in taxonomies)
    taxonomy =  common_prefix(taxonomies_list)
    
    missing_taxon_count = 7 - len(taxonomy)
    
    return ';'.join(taxonomy + ['Multi-affiliation']*missing_taxon_count)

def common_prefix(its):
    return [items[0] for items in takewhile(all_equal, zip(*its))]

     
def all_equal(items) -> bool:
    '''
    A helper function to test if 
    all items in the given iterable 
    are identical. 

    Arguments:
    item -> the given iterable to be used

    eg.
    >>> all_same([1, 1, 1])
    True
    >>> all_same([1, 1, 2])
    False
    >>> all_same((1, 1, 1))
    True
    >> all_same((1, 1, 2))
    False
    >>> all_same("111")
    True
    >>> all_same("112")
    False
    '''
    return all(item == items[0] for item in items)

