#!/usr/bin/env python

import re
import sys
import pandas as pd

# 7 letter abreviations of Pathogens as found in the FASTQ sequence headers
names = {
    'Alt.alt':'Alternaria alternata',
    'Alt.sol':'Alternaria solani',
    'Ca.Lib.':'Ca Liberibacter solanacearum',
    'Cla.mic':'Clavibacter michiganensis sepedonicus',
    'Colleto':'Colletotrichum coccodes',
    'Dic.dad':'Dickeya dadantii',
    'Dic.dia':'Dickeya dianthicola',
    'Dic.sol':'Dickeya solani',
    'Fus.oxy':'Fusarium oxysporum', #Fusarium solani
    'Fus.sam':'Fusarium sambucinum',
    'Glo.ell':'Globodera ellingtonae',
    'Glo.pal':'Globodera pallida',
    'Glo.ros':'Globodera rostachiensis',
    'Hel.sol':'Helminthosporium solani',
    'Mel.chi':'Meloidogyne chitwoodi',
    'Mel.fal':'Meloidogyne fallax',
    'Mel.hap':'Meloidogyne hapla',
    'Mel.inc':'Meloidogyne incognita',
    'Pec.atr':'Pectobacterium atrosepticum',
    'Pec.car':'Pectobacterium carotovora',
    'Pyt.ery':'Phytophthora erythroseptica',
    'Phy.inf':'Phytophthora infestans',
    'PLRV.GC':'Potato leaf roll virus, PLRV',
    'PMTV.GC':'Potato mop top virus, PMTV',
    'PVY.GCA':'Potato virus Y, PVY',
    'Pra.neg':'Pratylenchus neglectus',
    'Pra.pen':'Pratylenchus penetrans',
    'Pyt.ult':'Pythium ultimum',
    'Ral.sol':'Ralstonia solanacearum',
    'Rhiz.so':'Rhizoctonia solani',
    'Scl.scl':'Sclerotinia sclerotiorum',
    'Ath.rol':'Athelia rolfsii',
    'Spo.sub':'Spongospora subterranea',
    'Str.sca':'Streptomyces spp.',
    'Syn.end':'Synchytrium endobioticum',
    'TRV.GCA':'Tobacco Rattle Virus, TRV',
    'Ver.dah':'Verticillium dahliae',
    'Stubero':'Stuberosum',
}

# Pathogen Dataframe
dfPathogens = pd.DataFrame(index=[
    'total sequences (count)',
    'Total mapped to potato (counts)',
    'Sequences post potato removal (counts)',
    'Sequences that mapped to pathogens (counts)',
    "Sequences that didn't match to anything",
    '% total reads to potato',
    '% total reads to pathogens (post potato removal)',
    "% that didn't map to potato or pathogen"])

# Method to add data to the global Dataframe
def add_to_df(pathogens:dict, count:int, count_potato:int, unmapped_potato:int, count_pathogen:int, unmapped:int):
    global dfPathogens
    df = pd.DataFrame.from_dict(pathogens, orient='index', columns=[name])
    df.loc['total sequences (count)'] = count
    df.loc['Total mapped to potato (counts)'] = count_potato
    df.loc['Sequences post potato removal (counts)'] = unmapped_potato
    df.loc['Sequences that mapped to pathogens (counts)'] = count_pathogen
    df.loc["Sequences that didn't match to anything"] = unmapped
    df.loc['% total reads to potato'] = round( 100.0 * count_potato / count, 2 )
    df.loc['% total reads to pathogens (post potato removal)'] = round( 100.0 * count_pathogen / unmapped_potato, 2 )
    df.loc["% that didn't map to potato or pathogen"] = round( 100.0 * unmapped / count, 2 )
    dfPathogens = pd.concat([dfPathogens, df], axis=1)


# MAIN Loop
with open(sys.argv[1]) as reader:
    name = ''
    pathogens = dict()
    for line in reader:
        # Sample start
        match = re.search('Sequences in ([A-Za-z0-9_]*): ([0-9]*)', line)
        if match:
            if name:
                add_to_df(pathogens, count, count_potato, unmapped_potato, count_pathogen, unmapped)
                pathogens = dict()
            name = match.group(1)
            count = int(match.group(2))
        # Mapped potato count
        match = re.search('Mapped potato reads: ([0-9]*)', line)
        if match:
            count_potato = int(match.group(1))
        # Unmapped potato count
        match = re.search('Unmapped potato reads: ([0-9]*)', line)
        if match:
            unmapped_potato = int(match.group(1))
        # Mapped Pathogen count
        match = re.search('Mapped pathogen reads: ([0-9]*)', line)
        if match:
            count_pathogen = int(match.group(1))
        # Unmapped reads
        match = re.search('Unmapped reads: ([0-9]*)', line)
        if match:
            unmapped = int(match.group(1))
        # Mapped pathogens
        match = re.search('([A-Za-z.]*):([0-9]{1,})', line)
        if match:
            pathogens[match.group(1)] = int(match.group(2))
    # add last group
    add_to_df(pathogens, count, count_potato, unmapped_potato, count_pathogen, unmapped)

# rename index based on names dictionary
dfPathogens.rename(index=names, inplace=True)
# fill all unknowns to zero
dfPathogens.fillna(0, inplace=True)
print(dfPathogens)
# save combined TSV file
dfPathogens.to_csv(sys.argv[2], sep='\t', index=True)
