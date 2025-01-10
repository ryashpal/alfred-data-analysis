import logging

log = logging.getLogger("EHR-ML")

import warnings

warnings.filterwarnings("ignore")


def run(cutoff, minTokenLength, minRatioDifference):

    import os
    from pathlib import Path
    from scipy import stats
    import pandas as pd

    demoTrainDf = pd.read_csv(os.environ['EHR_DATA_BASE'] + '/blood_pos_cohort_20240614/data/wb_365_wa_1/splits_v1/demographics/mortality_normal_train.csv', sep='\t')
    demoTestDf = pd.read_csv(os.environ['EHR_DATA_BASE'] + '/blood_pos_cohort_20240614/data/wb_365_wa_1/splits_v1/demographics/mortality_normal_test.csv', sep='\t')
    demoValidateDf = pd.read_csv(os.environ['EHR_DATA_BASE'] + '/blood_pos_cohort_20240614/data/wb_365_wa_1/splits_v1/demographics/mortality_normal_validate.csv', sep='\t')
    demoDf = pd.concat([demoTrainDf, demoTestDf, demoValidateDf], ignore_index=True)
    mappingDf = pd.read_csv(os.environ['GENOMICS_DATA_BASE'] + '/patient_tube_id_mapping_full.tsv', sep='\t')
    mappedJourneyTubeIdsDf = mappingDf[['tube_code', 'PATIENT_ID', 'EPISODE_ID']].drop_duplicates().merge(
        demoDf[['person_id', 'visit_occurrence_id', 'JOURNEY_ID']],
        how='inner',
        left_on=['PATIENT_ID', 'EPISODE_ID'],
        right_on=['person_id', 'visit_occurrence_id']
    ).drop(
        columns=['PATIENT_ID', 'EPISODE_ID', 'visit_occurrence_id']
    )[['person_id', 'JOURNEY_ID', 'tube_code']].drop_duplicates()

    overlappingFilesDir = Path(os.environ['GENOMICS_DATA_BASE'], 'genome_nlp_tokens', 'overlapping_with_annotations')

    overlappingDfList = []
    for overlappingFile in os.listdir(overlappingFilesDir):
        df = pd.read_csv(Path(overlappingFilesDir, overlappingFile), sep='\t', names=['contig_id', 'start_position', 'end_position', 'tokens', 'score', 'feature_type', 'id', 'name', 'gene', 'atributes'])
        df['tube_code'] = [overlappingFile.split('_')[0]]*df.shape[0]
        overlappingDfList.append(df)
    overlappingDf = pd.concat(overlappingDfList, ignore_index=True)
    overlappingDf = overlappingDf.merge(
        mappedJourneyTubeIdsDf,
        how='inner',
        on=['tube_code']
    )
    overlappingDf = overlappingDf.drop_duplicates()
    overlappingDf = overlappingDf[overlappingDf.tokens.apply(lambda x: (len(x) >= minTokenLength))].reset_index()
    overlappingDf['gene'] = overlappingDf.gene.str.lower()
    overlappingDf = overlappingDf[overlappingDf.feature_type.isin(['CDS', 'ncRNA', 'oriC', 'regulatory_region', 'oriT'])]

    cutoffValue = overlappingDf.score.mean() + cutoff * overlappingDf.score.std()

    highScoreOverlappingDf = overlappingDf[(overlappingDf.score > cutoffValue)].tokens.value_counts().reset_index()
    lowScoreOverlappingDf = overlappingDf[(overlappingDf.score < cutoffValue)].tokens.value_counts().reset_index()
    mergedOverlappingDf = highScoreOverlappingDf.add_suffix('_hs').merge(
        lowScoreOverlappingDf.add_suffix('_ls'),
        how='inner',
        left_on=['tokens_hs'],
        right_on=['tokens_ls']
    )[['tokens_hs', 'count_hs', 'count_ls']].rename(columns={'tokens_hs': 'tokens'})

    mergedOverlappingDf['proportion_ls'] = mergedOverlappingDf.count_ls/mergedOverlappingDf.count_ls.sum()

    mergedOverlappingDf['count_expected'] = mergedOverlappingDf.proportion_ls * mergedOverlappingDf.count_hs.sum()

    filteredOverlappingDf = mergedOverlappingDf[(mergedOverlappingDf.count_hs >= 5) & (mergedOverlappingDf.count_ls >= 5)]

    chi2, p, dof, expected = stats.chi2_contingency(pd.crosstab(filteredOverlappingDf.count_hs, filteredOverlappingDf.count_expected), correction=True)
    significant = p < 0.05  # 5% significance level
    print(chi2, p, significant)

    filteredOverlappingDf.loc[:, 'ratio_difference'] = (filteredOverlappingDf.count_hs / filteredOverlappingDf.count_expected)

    overrepresentedTokensDf = filteredOverlappingDf.sort_values(by=['ratio_difference'])[filteredOverlappingDf.ratio_difference > minRatioDifference]

    mergedDf = overrepresentedTokensDf.merge(overlappingDf[['tokens', 'feature_type', 'id', 'name', 'gene']], how='inner', on=['tokens'])

    print(mergedDf[['tokens', 'feature_type', 'gene', 'name']].groupby(by=['feature_type', 'gene', 'name']).agg('count').reset_index().sort_values(by=['tokens'], ascending=False)[:15])

def main():

    import logging
    import sys
    import argparse

    log = logging.getLogger("EHR-ML")
    log.setLevel(logging.INFO)
    format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(format)
    log.addHandler(ch)

    log.info("Parsing command line arguments")

    parser = argparse.ArgumentParser(description='An utility to perform overlap of annotations to the tokens file')

    parser.add_argument('-c', '--cutoff', nargs=1, type=int, default=2, help='Cutoff value to define the high-attribution score. By default: [cutoff=2 SD]')

    parser.add_argument('-c', '--min_token_length', nargs=1, type=int, default=5, help='Minimum token length to retain for the analysis. By default: [min_token_length=5]')

    parser.add_argument('-c', '--min_ratio_difference', nargs=1, type=int, default=4, help='Minimum ratio difference to display the output. By default: [min_ratio_difference=4]')

    args = parser.parse_args()

    log.info('args.cutoff: ' + str(args.cutoff[0]))
    log.info('args.min_token_length: ' + str(args.min_token_length[0]))
    log.info('args.min_ratio_difference: ' + str(args.min_ratio_difference[0]))

    run(args.cutoff[0], args.min_token_length[0], args.min_ratio_difference[0])


if __name__ == "__main__":
    main()     
