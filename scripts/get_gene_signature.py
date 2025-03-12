print('Start!!')

import os
import pandas as pd
from pathlib import Path


bedDir = Path(os.environ['GENOMICS_DATA_BASE'], 'genome_nlp_tokens', 'bed_files', 'ECOLI')

print('Reading bed files')

bedDfList = []
for bedFile in os.listdir(bedDir):
    tubeid = bedFile.split('.')[0].split('_')[0]
    tokensDf = pd.read_csv(Path(bedDir, bedFile), sep='\t', names=['contig_id', 'start_position', 'end_position', 'tokens', 'score'])
    tokensDf['tube_id'] = tubeid
    bedDfList.append(tokensDf)

bedDf = pd.concat(bedDfList, ignore_index=True)

print('Calculating token lengths')

bedDf['token_length'] = bedDf.tokens.apply(lambda x: len(x))

highscoreCutoff = bedDf.score.mean() + 5 * bedDf.score.std()

highscoreDf = bedDf[bedDf.score > highscoreCutoff]

print('High score token lengths: ', highscoreDf.token_length.value_counts())

print('Reading overlapping files')

overlappingFilesDir = Path(os.environ['GENOMICS_DATA_BASE'], 'genome_nlp_tokens', 'overlapping_with_annotations', 'ECOLI')

overlappingDfList = []
for overlappingFile in os.listdir(overlappingFilesDir):
    df = pd.read_csv(Path(overlappingFilesDir, overlappingFile), sep='\t', names=['contig_id', 'start_position', 'end_position', 'tokens', 'score', 'feature_type', 'id', 'name', 'gene', 'atributes'])
    df['tube_code'] = [overlappingFile.split('_')[0]]*df.shape[0]
    overlappingDfList.append(df)
overlappingDf = pd.concat(overlappingDfList, ignore_index=True)

print('Obtain annotation overlapping with top tokens')

for tokenLengthCutoff in [5, 6, 7, 12, 13, 14]:
    print('tokenLengthCutoff: ', str(tokenLengthCutoff))
    tokensDf = highscoreDf[highscoreDf.token_length == tokenLengthCutoff]
    signatureTokens = tokensDf.tokens.unique()
    tokenOverlapsDfs = []
    for signatureToken in signatureTokens:
        tokenOverlapsDf = overlappingDf[overlappingDf.tokens == signatureToken]
        tokenOverlapsDfs.append(tokenOverlapsDf)
    if tokenOverlapsDfs:
        finalTokenOverlapsDf = pd.concat(tokenOverlapsDfs, ignore_index=True)
        finalTokenOverlapsDf.drop_duplicates().to_csv(Path(os.environ['GENOMICS_DATA_BASE'], 'survival_analysis', 'ECOLI', 'constant_score_overlaping_annotations', 'sd_5_tl_' + str(tokenLengthCutoff) + '.csv'), index=False)

print('Done!!')
