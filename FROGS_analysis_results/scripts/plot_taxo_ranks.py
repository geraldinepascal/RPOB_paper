
# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import pandas as pd
import plotly.express as px
import os
import re


def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('--affi_tables', nargs="+", required=True, 
                        help='Affiliation abundance table with multiaffi debug of the 16S23S db affiliation, output of the script add_multiaffi_to_abd_table.py ')
    parser.add_argument('--labels', nargs="+", required=False, 
                        help='Name associated to table. order maters. Name of the table is used by default.')

    parser.add_argument('-s', '--samples', nargs="+", required=True, 
                        help='Sample names in the table. Can be a list of names or you can use 1:32 to express sample from 1 to 32 included. 1:5 8 19:21 would be valid and would plot samples 1 2 3 4 5 8 19 20 21')

    parser.add_argument('-o', '--outdir', default='.',
                        help='Path of the output dir.')

    parser.add_argument('-f', '--outformat', nargs="+", default={'html', 'svg'}, type=str, choices=['png', 'jpg', 'jpeg', "webp", 'svg', "pdf", "html"],
                        help='Format of the output plots.')

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    parser.add_argument("--debug", help="increase a lot output verbosity",
                        action="store_true")

    args = parser.parse_args()
    return args




def main():

    args = parse_arguments()

    if args.debug:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.debug('Mode debug ON')
    elif args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")


    # rank2color = {'superkingdom': 'rgb(95, 70, 144)',
    #                 "Domain":'rgb(95, 70, 144)',
    #                 'phylum': 'rgb(29, 105, 150)',
    #                 'class': 'rgb(56, 166, 165)',
    #                 'order': 'rgb(15, 133, 84)',
    #                 'family': 'rgb(115, 175, 72)',
    #                 'genus': 'rgb(237, 173, 8)',
    #                 'species': 'rgb(204, 80, 62)',
    #                 'Superkingdom': 'rgb(95, 70, 144)',
    #                 'Phylum': 'rgb(29, 105, 150)',
    #                 'Class': 'rgb(56, 166, 165)',
    #                 'Order': 'rgb(15, 133, 84)',
    #                 'Family': 'rgb(115, 175, 72)',
    #                 'Genus': 'rgb(237, 173, 8)',
    #                 'Species': 'rgb(204, 80, 62)', 
    #                 "weak affiliation":"grey"}
    col2code = {"grey":"#999999", 
                    "orange":"#E69F00",
                    "ligth_bleu":"#56B4E9",
                    "green":"#009E73",
                    "yellow":"#F0E442",
                    "darker_blue":"#0072B2",
                    "darker_orange":"#D55E00",
                    "pink":"#CC79A7"}

    rank2color = {"Domain":col2code['pink'],
                    'Phylum': col2code['yellow'],
                    'Class': col2code['ligth_bleu'],
                    'Order': col2code['darker_blue'],
                    'Family': col2code['green'],
                    'Genus': col2code['orange'],
                    'Species': col2code['darker_orange'],
                    "weak affiliation":col2code['grey'],}
    
    ranks = [ "Domain", "Phylum",  "Class", "Order",
                "Family", "Genus", "Species", "Species"]

# color  1-main 000000   0   0   0 black, grey, grey, grey, rich black, grey, cod grey, grey, almost black, grey
# color  2-main 2271B2  34 113 178 honolulu blue, bluish, strong cornflower blue, spanish blue, medium persian blue, sapphire blue, ocean boat blue, french blue, windows blue, tufts blue
# color  2-alt  AA0DB4 170  13 180 barney, strong magenta, heliotrope magenta, strong heliotrope, steel pink, barney purple, purple, violet, violet eggplant, deep magenta
# color  3-main 3DB7E9  61 183 233 summer sky, cyan, picton blue, vivid cerulean, deep sky blue, brilliant cornflower blue, malibu, bright cerulean, cerulean, cerulean
# color  3-alt  FF54ED 255  84 237 light magenta, violet pink, light brilliant magenta, pink flamingo, light brilliant orchid, brilliant magenta, purple pizzazz, candy pink, blush pink, shocking pink
# color  4-main F748A5 247  72 165 barbie pink, rose bonbon, wild strawberry, brilliant rose, brilliant rose, magenta, wild strawberry, light brilliant rose, frostbite, brilliant cerise
# color  4-alt  00B19F   0 177 159 strong opal, tealish, persian green, keppel, topaz, manganese blue, light sea green, sea green light, puerto rico, turquoise
# color  5-main 359B73  53 155 115 ocean green, sea green, viridian, mother earth, moderate spring green, moderate aquamarine, paolo veronese green, observatory, jungle green, ocean green
# color  5-alt  EB057A 235   5 122 vivid rose, red purple, mexican pink, bright pink, rose, strong pink, luminous vivid rose, deep pink, winter sky, hot pink
# color  6-main d55e00 213  94   0 bamboo, smoke tree, red stage, tawny, tenn, tenne, burnt orange, rusty orange, dark orange, mars yellow
# color  6-alt  F8071D 248   7  29 vivid red, luminous vivid amaranth, ruddy, ku crimson, vivid amaranth, light brilliant red, cherry red, red, red, bright red
# color  7-main e69f00 230 159   0 gamboge, squash, buttercup, marigold, dark goldenrod, medium goldenrod, fuel yellow, sun, harvest gold, orange
# color  7-alt  FF8D1A 255 141  26 dark orange, juicy, west side, tangerine, gold drop, pizazz, princeton orange, university of tennessee orange, tangerine, tahiti gold
# color  8-main f0e442 240 228  66 holiday, buzz, paris daisy, starship, golden fizz, dandelion, gorse, lemon yellow, bright lights, sunflower
# color  8-alt  9EFF37 158 255  55 french lime, lime, green yellow, green lizard, luminous vivid spring bud, spring frost, vivid spring bud, bright yellow green, spring bud, acid green


    affi_tables =  args.affi_tables
    labels =  args.labels 
    output_formats = set(args.outformat)
    outdir = args.outdir

    samples = []
    for s in args.samples:
        if re.match("^\d+\:\d+$", s):
            start, end = s.split(':')
            samples_range = [str(i) for i in range(int(start), int(end)+1)]
            samples += samples_range
            logging.info(f'from {s} to {samples_range}')
        else:
            samples.append(s)

    logging.info(f'Going to plot {len(samples)} samples : {samples}')    

    if not labels:
        labels = affi_tables

    logging.info(f'Going to plot {len(labels)} analysis : {"    ".join(labels)}')   
    
    ## Merge all table in one
    df_list = []
    for table, name in zip(affi_tables, labels):
        logging.info(f'Processing {table} labeled {name}')
        df = pd.read_csv(table, sep='\t')
        df['name'] = name
        df_list.append(df)

    df = pd.concat(df_list)
    

    # give a 'weak affiliation' label to the cluster that are have not a valid affi
    filt = df['valid_affiliation']
    df.loc[:, 'Taxonomic rank'] = df['rank_affi']
    df.loc[~filt, 'Taxonomic rank'] = f"weak affiliation"



    ### PLOT ALL SAMPLES IN ONE BAR
    # Group by rank 

    df_rank = df.groupby(['Taxonomic rank', "region", "db", "name"
                        ]).agg({"observation_sum":"sum", 'abundance':"sum",
                                                    "observation_name":"count", }).reset_index()

                                        
    df_rank['abundance_all_sample_round']  = df_rank['abundance'].round(2).astype(str) +"%"
    df_rank['abundance_all_sample_round'] = df_rank['abundance_all_sample_round'] + '<br>' + df_rank['observation_name'].astype(str) + ' clusters'

    fig = px.bar(df_rank, x="name", y="abundance", color="Taxonomic rank", # facet_col="name",
                 template="seaborn",
                text="abundance_all_sample_round",
                category_orders={"rank_affi":ranks, "name":labels, #"region":regions_analysed,
                                'Taxonomic rank':ranks + ["ident or cov < 99%"] }, color_discrete_map=rank2color, )
    
    fig.update_layout( # customize font and legend orientation & position
        legend={'traceorder':'reversed'}
    )
    
    #fig.update_layout(
    #    width=700,
    #    height=600,)

    output_base_name = os.path.join(outdir, "taxonomic_ranks_per_target")

    for extension in output_formats:
        logging.info(f'Writting {output_base_name}.{extension}')
        if extension == "html":
            fig.write_html(f"{output_base_name}.{extension}")
        else:
            fig.write_image(f"{output_base_name}.{extension}")
        

    ### PLOT ALL SAMPLES WITH ONE BAR PER SAMPLE
        

    df_rank_by_sample_list =[]
    for label in labels:
        filt_name = df['name'] == label
        

        for sample in samples:

            if str(sample) not in df.columns:
                logging.warning(f'sample {sample} is not found in the table.')
                continue

            df.loc[filt_name, 'abundance'] = 100 * df.loc[filt_name, f'{sample}']/df.loc[filt_name, f'{sample}'].sum()


            df.loc[filt_name, f'sample_sum'] = df.loc[filt_name, f'{sample}']

            df_rank_by_sample = df.loc[filt_name].groupby(['Taxonomic rank', "region", 'name']).agg({f"sample_sum":"sum", 
                                                    f'abundance':"sum",
                                                    "observation_name":"count",}).reset_index()
            df_rank_by_sample['sample'] = str(sample)
            df_rank_by_sample_list.append(df_rank_by_sample)
            
        
    df_rank = pd.concat(df_rank_by_sample_list)
    
    # Output tsv
    table_output_name = os.path.join(outdir, "taxonomic_ranks_per_target_and_per_sample.tsv")

    logging.info(f'Writting ranks table in {table_output_name}')
    
    df_rank.to_csv(table_output_name,  sep='\t', index=False)#, columns=[])


    df_rank["n"] = df_rank["name"]

    fig = px.bar(df_rank, x="sample", y="abundance", color="Taxonomic rank", facet_row="n",
             category_orders={"rank_affi":ranks, "n":labels,  # "dataset":dataset_analysed, "region":regions_analysed,
                              'Taxonomic rank':ranks + ["ident or cov < 99%"] },
                              template="seaborn",
             color_discrete_map=rank2color, )
    
    fig.update_layout( # customize font and legend orientation & position
        legend={'traceorder':'reversed'}
    )

    output_base_name = os.path.join(outdir, "taxonomic_ranks_per_target_and_per_sample")

    for extension in output_formats:
        logging.info(f'Writting {output_base_name}.{extension}')
        if extension == "html":
            fig.write_html(f"{output_base_name}.{extension}")
        else:
            fig.write_image(f"{output_base_name}.{extension}")


    # Raw 
    fig = px.bar(df_rank, x="sample", y="sample_sum", color="Taxonomic rank", facet_row="n",
            category_orders={"rank_affi":ranks, "n":labels,  # "dataset":dataset_analysed, "region":regions_analysed,
                            'Taxonomic rank':ranks + ["ident or cov < 99%"] },
                            template="seaborn",
            color_discrete_map=rank2color, )
    
    fig.update_layout( # customize font and legend orientation & position
        legend={'traceorder':'reversed'}
    )

    output_base_name = os.path.join(outdir, "taxonomic_ranks_per_target_and_per_sample_raw_count")

    for extension in output_formats:
        logging.info(f'Writting {output_base_name}.{extension}')
        if extension == "html":
            fig.write_html(f"{output_base_name}.{extension}")
        else:
            fig.write_image(f"{output_base_name}.{extension}")




if __name__ == '__main__':
    main()
