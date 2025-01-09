import logging

log = logging.getLogger("EHR-ML")

import warnings

warnings.filterwarnings("ignore")


def run(tubeCode, tokensDir, annotationsDir, outputDir):

    from pathlib import Path
    import pandas as pd


    annotationFile = Path(annotationsDir, tubeCode + '.gff3')

    log.info('Reading annotation from the file: ' + str(annotationFile))

    annotationDf = pd.read_csv(
        annotationFile,
        sep='\t',
        comment='#',
        names=['sequence_id', 'source', 'feature_type', 'feature_start', 'feature_end', 'score', 'strand', 'phase', 'atributes']
    )
    annotationDf['id'] = annotationDf.atributes.apply(lambda x: [att.split('=')[1] for att in x.split(';') if att.split('=')[0] == 'ID']).apply(lambda x: None if (len(x) == 0) else x[0])
    annotationDf['name'] = annotationDf.atributes.apply(lambda x: [att.split('=')[1] for att in x.split(';') if att.split('=')[0] == 'Name']).apply(lambda x: None if (len(x) == 0) else x[0])
    annotationDf['gene'] = annotationDf.atributes.apply(lambda x: [att.split('=')[1] for att in x.split(';') if att.split('=')[0] == 'gene']).apply(lambda x: None if (len(x) == 0) else x[0])
    annotationDf = annotationDf[annotationDf.feature_type != 'region']

    tokensFile = Path(tokensDir, tubeCode + '_short_token_attributions.bed')

    log.info('Reading tokens from the file: ' + str(tokensFile))

    tokensDf = pd.read_csv(tokensFile, sep='\t', names=['contig_id', 'start_position', 'end_position', 'tokens', 'score'])

    log.info('Obtaining the overlaps')

    mergedDf = tokensDf.merge(
        annotationDf,
        how='inner',
        left_on=['contig_id'],
        right_on=['sequence_id']
    )

    overlappingDf = mergedDf[(mergedDf.start_position >= mergedDf.feature_start) & (mergedDf.end_position <= mergedDf.feature_end)]

    overlappingDf = overlappingDf[['contig_id', 'start_position', 'end_position', 'tokens', 'score_x', 'id', 'name', 'gene', 'atributes']]

    Path(outputDir).mkdir(parents=True, exist_ok=True)

    overlapFile = Path(outputDir, tubeCode + '_annotation_overlaps.bed')

    log.info('Saving the overlaps to the file: ' + str(overlapFile))

    overlappingDf.to_csv(overlapFile, sep='\t', index=False, header=None)


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

    parser.add_argument('tube_code', nargs=1, help='Tube code')

    parser.add_argument('tokens_dir', nargs=1, help='Path to a directory containing tokens in bed format')

    parser.add_argument('annotations_dir', nargs=1, help='Path to a directory containing annotation files')

    parser.add_argument('output_dir', nargs=1, help='Path to a directory to save the output containing overlapping annotations for the tokens')

    args = parser.parse_args()

    log.info('args.tube_code: ' + str(args.tube_code[0]))
    log.info('args.tokens_dir: ' + str(args.tokens_dir[0]))
    log.info('args.annotations_dir: ' + str(args.annotations_dir[0]))
    log.info('args.output_dir: ' + str(args.output_dir[0]))

    run(args.tube_code[0], args.tokens_dir[0], args.annotations_dir[0], args.output_dir[0])

if __name__ == "__main__":
    main()     
