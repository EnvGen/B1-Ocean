#!/usr/bin/env python

import pandas as pd
from ete3 import NCBITaxa
import os
import sys
from pathlib import Path
from argparse import ArgumentParser


def init_sqlite_taxdb(taxdir, force=False):
    """Creates ete3 sqlite database"""
    taxdump_file = f"{taxdir}/taxdump.tar.gz"
    dbfile = f"{taxdir}/taxonomy.sqlite"
    if not os.path.exists(taxdump_file) or os.path.exists(dbfile) and not force:
        taxdump_file = None
    if not os.path.exists(taxdir):
        os.makedirs(taxdir)
    Path(dbfile).touch()
    ncbi_taxa = NCBITaxa(dbfile=dbfile, taxdump_file=taxdump_file)
    return ncbi_taxa


def translate_taxids_to_names(res_df, reportranks, name_dict):
    """
    Takes a pandas dataframe with ranks as columns and contigs as rows and
    taxids as values and translates taxids to names column by column using a
    taxid->name dictionary

    Parameters
    ----------
    res_df: pandas.DataFrame
        Results with taxids
    reportranks: list
        List of taxonomic ranks to report results for
    name_dict: dictionary
        Dictionary mapping taxids -> names
    Returns
    -------
    res: pandas.DataFrame
        Dataframe with names instead of taxids
    """

    res = {}
    for rank in reportranks:
        res[rank] = [name_dict[taxid] for taxid in res_df.loc[:, rank]]
    res = pd.DataFrame(res)
    res.index = res_df.index
    res = res.loc[:, reportranks]
    return res


def make_name_dict(df, ranks):
    """
    Creates a dictionary of taxids to taxonomy names, including Unclassified
    ranks

    Parameters
    ----------
    df: pandas.DataFrame
        Lineage dataframe
    ranks: list
        Ranks to store names information for

    Returns
    -------
    name_dict: dict
        Name dictionary mapping taxonomy ids to names
    """

    name_dict = {}
    for rank in ranks:
        name_dict.update(
            dict(zip(df[rank].values, df["{}.name".format(rank)].values)))
        name_dict.update(dict(
            zip(-abs(df[rank]), "Unclassified." + df["{}.name".format(rank)])))
    name_dict[-1] = "Unclassified"
    return name_dict


def unique_cols(x):
    """
    Checks the columns for a taxid entry and only returns unique indices
    This is meant to fix the problem of some taxids having multiple entries
    for the same taxonomic rank.

    Example: Alouatta palliata mexicana (Howler monkey)
    taxid = 182248
            superkingdom phylum order genus species  class  family    class
            2759         7711   9443  9499  30589    40674  378855    1338369

    This function would return indices for columns to use with 'iloc' in order
    to get:
            superkingdom phylum order genus species  class  family
            2759         7711   9443  9499  30589    40674  378855
    :param x:
    :return:
    """
    col_index = []
    cols = []
    for i, c in enumerate(x.columns):
        if c not in cols:
            cols.append(c)
            col_index.append(i)
    return col_index


def add_names(x, taxid, ncbi_taxa):
    """
    This function translates taxonomy ids to names. It operates per-row in the
    lineage dataframe.

    Parameters
    ----------
    x: pandas.DataFrame
        DataFrame of one taxid and its taxonomic ranks
    taxid: int
        Taxid being evaluated
    ncbi_taxa: ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        The ete3 sqlite database connection

    Returns
    -------
        The original DataFrame merged with the taxa names
    """

    # Get a names dictionary for all taxids in the row
    names = ncbi_taxa.get_taxid_translator(list(x.loc[taxid].values) + [taxid])
    n = {}
    # Iterate ranks
    for rank in list(x.columns):
        # Get taxid for the current rank
        t = x.loc[taxid, rank]
        # If taxid is negative it means that there is no classified taxonomy
        # at this rank. Instead we get the last known name in the hierarchy.
        # We can then use the negative values to translate into the name with
        # the "Unclassified." prefix. If the name is 'root' we just use
        # 'Unclassified'
        if t < 0:
            known_name = names[-t]
            if known_name == "root":
                name = "Unclassified"
            else:
                name = known_name
        # If taxid is positive we just use the name from the dictionary
        else:
            name = names[t]
        # Add name to a dictionary with keys in the form of {rank}.name
        n["{}.name".format(rank)] = name
    name_df = pd.DataFrame(n, index=[taxid])
    return pd.merge(x, name_df, left_index=True, right_index=True)



def propagate_lower(x, taxid, ranks):
    """
    Shift known ranks down through the taxonomic hierarchy.

    Parameters
    ----------
    x: pandas.DataFrame
        DataFrame of one taxid and its taxonomic ranks
    taxid:  int
        Taxid being evaluated
    ranks: list
        Ranks used for assigning

    Returns
    -------
        pandas.DataFrame updated with missing ranks

    Some proteins in the database may map to a taxonomic rank above the
    lowest taxonomic rank that we are trying to assign. For instance,
    if we use the ranks 'superkingdom phylum genus species' and a protein
    maps to a taxid at rank phylum then we want to add the taxonomic
    information at the genus and species levels. This is done here by adding
    the negative taxid of the lowest known rank to the lower ranks.

    Example: In the Uniref90 database the entry 'E1GVX1' maps to taxonomy id
    838 (rank: genus, name: Prevotella). When creating the lineage for taxid
    838 we add '-838' to rank species.
    """

    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    missing = {}
    known = taxid
    for rank in rev_ranks[0:]:
        if rank not in x.columns:
            missing[rank] = -known
        else:
            known = x.loc[taxid, rank]
    return pd.merge(x, pd.DataFrame(missing, index=[taxid]), left_index=True,
                    right_index=True)


def process_lineage(lineage, taxid, ncbi_taxa, ranks):
    """
    Looks up lineage information from taxids.

    The lineage object is a list of taxonomic ids corresponding to the full
    lineage of a single taxid.
    """
    # Get ranks for each taxid in the lineage
    lineage_ranks = ncbi_taxa.get_rank(lineage)
    x = pd.DataFrame(lineage_ranks, index=["rank"]).T
    x = x.loc[x["rank"].isin(ranks)].reset_index().T
    x.columns = x.loc["rank"]
    x.drop("rank", inplace=True)
    x.index = [taxid]
    # Only select unique columns
    x = x.iloc[:, unique_cols(x)]
    # Add taxids for lower ranks in the hierarchy
    x = propagate_lower(x, taxid, ranks)
    # Add names for taxids
    x = add_names(x, taxid, ncbi_taxa)
    return x


def make_lineage_df(taxids, ncbi_taxa, ranks):
    """
    Creates a lineage dataframe with full taxonomic information for a list of
    taxids.

    Example:
    taxid   species phylum  genus   genus.name      phylum.name     species.name
    859655  305     1224    48736   Ralstonia       Proteobacteria  Ralstonia solanacearum
    387344  1580    1239    1578    Lactobacillus   Firmicutes      Lactobacillus brevis
    358681  1393    1239    55080   Brevibacillus   Firmicutes      Brevibacillus brevis

    Parameters
    ----------
    taxids: list
        List of taxonomic ids to obtain information for
    taxdir: str
        Path to directory holding taxonomic info
    dbname: str
        Name of ete3 sqlite database within taxdir
    ranks: list
        Ranks to store information for
    Returns
    -------
    lineage_df: pandas.DataFrame
        Data Frame with full taxonomic info
    """
    lineages = ncbi_taxa.get_lineage_translator(taxids)
    # Store potential missing taxids and warn user
    missing_taxids = set([int(x) for x in taxids]).difference(lineages.keys())
    # Get possible translations for taxids that have been changed
    _, translate_dict = ncbi_taxa._translate_merged(
        list(set(taxids).difference(lineages.keys())))
    rename = {y: x for x, y in translate_dict.items()}
    # Update lineages with missing taxids
    lineages.update(ncbi_taxa.get_lineage_translator(translate_dict.values()))
    res = []
    for taxid, lineage in lineages.items():
        res.append(process_lineage(lineage, taxid, ncbi_taxa, ranks))
    lineage_df = pd.concat(res, sort=False)
    lineage_df.rename(index=rename, inplace=True)
    lineage_df.rename(index=lambda x: int(x), inplace=True)
    for rank in ranks:
        lineage_df[rank] = pd.to_numeric(lineage_df[rank])
    name_dict = make_name_dict(lineage_df, ranks)
    if len(missing_taxids) > 0:
        sys.stderr.write("#WARNING: Missing taxids found:\n")
        sys.stderr.write(
            "#{}\n".format(",".join([str(x) for x in missing_taxids])))
        sys.stderr.write(
            "#To fix this, you can try to update the "
            "taxonomy database.")
    return lineage_df.loc[:, lineage_df.dtypes == int], name_dict


def main(args):
    ncbi_taxa = init_sqlite_taxdb(args.taxdir, args.force)
    df = pd.read_csv(args.infile,
                     header=None, sep="\t", index_col=1, usecols=[0, 1, 2, 3],
                     names=["result", "contig_id", "taxid", "length"])
    lineage_df, name_dict = make_lineage_df(
        df.loc[df.result == "C"].taxid.unique(), ncbi_taxa, args.ranks)
    df_taxid = pd.merge(df, lineage_df, left_on="taxid", right_index=True,
                        how="outer")
    df_taxid.fillna(-1, inplace=True)
    results = translate_taxids_to_names(df_taxid, reportranks=args.ranks, name_dict=name_dict)
    results.to_csv(args.outfile, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", help="Kraken results file")
    parser.add_argument("outfile", help="Outfile")
    parser.add_argument("taxdir", help="Directory for NCBI taxonomy info")
    parser.add_argument("--ranks", help="Ranks to report",
                        nargs="+", default=["superkingdom", "phylum", "class", "order", "family", "genus", "species"])
    parser.add_argument("--dbname", help="Name of taxonomy file",
                        default="taxonomy")
    parser.add_argument("--force", action="store_true",
                        help="Force download of taxonomy")
    args = parser.parse_args()
    main(args)