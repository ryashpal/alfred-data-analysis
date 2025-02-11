if __name__ == "__main__":
    import os
    inDir = os.environ['GENOMICS_DATA_BASE'] + '/annotations/e_coli/'
    outDir = os.environ['GENOMICS_DATA_BASE'] + '/annotations/e_coli_gff3'
    for file in os.listdir(inDir):
        print('file: ', file)
        if file.endswith('.gff3'):
            with open((inDir + '/' + file), mode='r') as in_file, open((outDir + '/' + file), mode='w') as out_file:
                for line in in_file.readlines():
                    if '##FASTA' in line:
                        break
                    out_file.write(line)
