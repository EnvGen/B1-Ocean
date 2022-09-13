#!/usr/bin/env python

import pandas as pd
import os


def orf2feat(d, val_name="vals", regex=""):
    import re
    keys = []
    vals = []
    for key, val in d.items():
        for item in val.split(","):
            if regex != "":
                m = re.search(regex, item)
                if m == None:
                    continue
                item = m.group()
            keys.append(key)
            vals.append(item)
    df = pd.DataFrame(data={"orf": keys, val_name: vals})
    df.set_index("orf", inplace=True)
    return df


def parse_emapper(sm):
    from pandas.errors import EmptyDataError
    db = sm.wildcards.db
    db_att = {'kos': {'val_name': 'ko', 'regex': "K\d{5}"},
              'pathways': {'val_name': 'pathway', 'regex': "map\d{5}"},
              'modules': {'val_name': 'module', 'regex': ""},
              'enzymes': {'val_name': 'enzyme', 'regex': ""}}
    df = pd.read_csv(sm.input.annotations, sep="\t", index_col=0)
    df.rename(columns={'KEGG_ko': 'kos', 'KEGG_Pathway': 'pathways',
                          'KEGG_Module': 'modules', 'EC': 'enzymes'},
               inplace=True)
    df.fillna("-", inplace=True)
    d = df.loc[df[db] != "-", db].to_dict()
    how = "inner"
    try:
        info_df = pd.read_csv(sm.input.info, sep="\t", index_col=0)
    except EmptyDataError:
        info_df = pd.DataFrame()
        how = "right"
    annot = orf2feat(d, val_name=db_att[db]["val_name"],
                     regex=db_att[db]["regex"])
    annot = pd.merge(info_df, annot, left_index=True,
                     right_on=db_att[db]["val_name"], how=how)
    annot.to_csv(sm.output[0], sep="\t", index=True, header=True)


def parse_rgi(sm):
    annot = pd.read_csv(sm.input.txt, sep="\t", index_col=0)
    annot = annot.loc[:, ["Model_ID", "AMR Gene Family", "Resistance Mechanism"]]
    annot.loc[:, "Model_ID"] = ["RGI_{}".format(x) for x in annot.Model_ID]
    annot.rename(index=lambda x: x.split(" ")[0], inplace=True)
    annot.to_csv(sm.output.tsv, sep="\t", index=True)


def parse_pfam(sm):
    annot = pd.read_csv(sm.input[0], comment="#", header=None, sep=" +",
                        usecols=[0, 5, 7, 14], engine="python",
                        names=["orf", "pfam", "pfam_type", "pfam_clan"])
    if os.path.exists(sm.input[1]) and os.path.exists(sm.input[2]):
        clan_files = True
        clans = pd.read_csv(sm.input[1], header=None, names=["clan", "clan_name"],
                            usecols=[0, 3], sep="\t")
        info = pd.read_csv(sm.input[2], header=None,
                            names=["pfam", "clan", "pfam_name"], usecols=[0, 1, 4],
                            sep="\t")
    else:
        clan_files = False
    # Strip suffix for pfams
    annot.loc[:, "pfam"] = [x.split(".")[0] for x in annot.pfam]
    # Select unique orf->pfam mappings
    # TODO: This masks multiple occurrences of domains on the same orf. Figure out if this is wanted or not.
    annot = annot.groupby(["orf", "pfam"]).first().reset_index()
    # Merge with pfam info and clan info
    cols = ["orf", "pfam", "pfam_name"]
    if clan_files:
        cols += ["clan", "clan_name"]
        annot = pd.merge(annot, info, left_on="pfam", right_on="pfam")
        annot = pd.merge(annot, clans, left_on="clan", right_on="clan", how="left")
        annot.fillna("No_clan", inplace=True)
    annot = annot.loc[:, cols]
    annot.sort_values("orf", inplace=True)
    # Write to file
    annot.to_csv(sm.output[0], sep="\t", index=False)


def read_hmmout(f, evalue_threshold=0.001):
    orf_hits = {}
    d = {}
    with open(f, 'r') as fhin:
        i = 0
        for line in fhin:
            if line.rstrip().startswith("#"):
                continue
            items = line.rstrip().rsplit()
            orf, name, acc, evalue, score, start, stop = (items[0], items[3], items[4], float(items[6]), float(items[7]),  int(items[15]), int(items[16]))
            # Filter by evalue
            if evalue > evalue_threshold:
                continue
            # For each hmm with '-' as hmm id and only a name, set hmm id to name
            if acc == "-":
                acc = name
            if not orf in orf_hits.keys():
                orf_hits[orf] = []
            if acc in orf_hits[orf]:
                continue
            orf_hits[orf].append(acc)
            d[i] = {'orf': orf, 'hmm': acc, "hmm_name": name, 'evalue': evalue,
                    'score': score, 'start': start, 'stop': stop}
            i+=1
        annot = pd.DataFrame(d).T
        annot.set_index("orf", inplace=True)
    return annot


def get_non_overlapping(annot):
    """
    Iterates a dataframe with start, stop coordinates and only
    selects highest scoring hits that do not overlap
    """
    d = {}
    orf_hits = {}
    i = 0
    # Sort the dataframe by score from highest to lowest
    for row in annot.sort_values("score", ascending=False).iterrows():
        orf = row[0]
        # orf_hits is a dictionary storing all positions occupied by
        # a hit of some kind
        try:
            orf_hits[orf]
        except KeyError:
            orf_hits[orf] = {'residues': []}
        r = row[1]
        # Generate a list of all residues occupied by this hit
        stretch = list(range(r["start"], r["stop"]+1))
        # Check overlap
        overlap = False
        if len(orf_hits[orf].keys()) > 0:
            residues = orf_hits[orf]["residues"]
            # Check intersection between current stretch and stored
            # occupied residues. If intersection is > 0 there's an overlap
            if len(set(stretch).intersection(residues)) > 0:
                overlap = True
        if not overlap:
            d[i] = {'orf': orf, 'hmm': r["hmm"], 'hmm_name': r["hmm_name"], 'evalue': r["evalue"], 'score': r["score"],
                    'start': r["start"], 'stop': r["stop"]}
            orf_hits[orf]["residues"] = sorted(orf_hits[orf]["residues"]+stretch)
            i += 1
    filtered = pd.DataFrame(d).T
    return filtered.set_index("orf")


def parse_hmmsearch(sm):
    f = sm.input[0]
    annot = read_hmmout(f, evalue_threshold=sm.params.evalue)
    filtered = get_non_overlapping(annot)
    df = filtered.loc[:, ["hmm", "hmm_name"]]
    df.to_csv(sm.output[0], sep="\t", index=True)


def main(sm):
    toolbox = {"parse_pfam": parse_pfam,
               "parse_hmmsearch": parse_hmmsearch,
               "parse_emapper": parse_emapper,
               "parse_rgi": parse_rgi}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
