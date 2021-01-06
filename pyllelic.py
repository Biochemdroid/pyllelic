#!/usr/bin/env python3
"""pyllelic: a tool for detection of allelic-specific varation in DNA sequencing.

## Example usage in ipython / jupyter notebook:

    import pyllelic

    pyllelic.set_up_env_variables(
        base_path="/Users/abonham/documents/test_allelic/",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
    )

    pyllelic.main("output.xlsx")  # runs every step all at once

----------------------------------

## Example exploratory / step-by-step use in ipython / jupyter notebook:

    import pyllelic

    pyllelic.set_up_env_variables(
        base_path="/Users/abonham/documents/test_allelic/",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
    )

    files_set = pyllelic.make_list_of_bam_files()  # finds bam files

    positions = pyllelic.index_and_fetch(files_set)  # index bam and creates bam_output folders/files

    pyllelic.genome_parsing()  # writes out genome strings in bam_output folders

    cell_types = extract_cell_types(files_set)  # pulls out the cell types available for analysis

    df_list = run_quma_and_compile_list_of_df(cell_types, filename)  # run quma, get dfs

    means_df = process_means(df_list, positions, files_set)  # process means data from dataframes

    modes_df = process_modes(df_list, positions, cell_types)  # process modes data from dataframes

    diff_df = find_diffs(means_df, modes_df)  # find difference between mean and mode

    write_means_modes_diffs(means_df, modes_df, diffs_df, filename)  # write output to excel files
"""

# Imports
import pandas as pd
import numpy as np
import pysam
import os
import sys
from skbio.alignment import StripedSmithWaterman
import plotly.express as px  # noqa
import subprocess
from pathlib import Path
from scipy import stats
from tqdm.notebook import tqdm
import config


def set_up_env_variables(base_path, prom_file, prom_start, prom_end, chrom):
    """Helper method to set up all our environmental variables, such as for testing.

    Args:
        base_path (str): directory where all processing will occur, put .bam files in "test"
                         sub-directory in this folder
        prom_file (str): filename of genmic sequence of promoter region of interest
        prom_start (str): start position to analyze in promoter region
        prom_end (str): final position to analyze in promoter region
        chrom (str): chromosome promoter is located on
    """

    config.base_directory = Path(base_path)
    config.promoter_file = Path(base_path) / prom_file
    config.results_directory = Path(base_path) / "results"
    config.bam_directory = Path(base_path) / "bam_output"
    config.analysis_directory = Path(base_path) / "test"
    config.promoter_start = prom_start
    config.promoter_end = prom_end
    config.chromosome = chrom


##################################################################################
##################################################################################
def main(filename):
    """Run a given set of Pyllelic analysis, using values from supplied environmental variables.

    Args:
        filename: filename to write output of analysis to.
    """

    files_set = make_list_of_bam_files()
    positions = index_and_fetch(files_set)
    genome_parsing()
    cell_types = extract_cell_types(files_set)
    df_list = run_quma_and_compile_list_of_df(cell_types, filename)
    means_df = process_means(df_list, positions, cell_types)
    modes_df = process_modes(df_list, positions, cell_types)
    diffs_df = find_diffs(means_df, modes_df)
    write_means_modes_diffs(means_df, modes_df, diffs_df, filename)


##################################################################################
##################################################################################


def genome_range(position, genome_string):  # TODO: This may be specific only for TERT
    """Helper to return a genome string (e.g., "ATCGACTAG")
    given a position and an entire string.

    Args:
        position: genomic position on chromesome
        genome_string: string representation of genomic promoter known sequence

    Returns:
        str: genomic bases for indicated read / position
    """

    start = 1298163 - (int(position) + 30)
    end = 1298163 - (int(position) + 1)

    return genome_string[start:end]


def run_quma(directory, genomic_seq_file, reads_seq_file):
    """Helper function to run external QUMA tool.

    Args:
        directory (str): directory path to analyze
        genomic_seq_file (str): text file with known genomic sequence
        reads_seq_file (str): text file with experimental reads

    Returns:
        bytes: shell output from quma command
    """

    quma_path = os.fspath(config.base_directory.joinpath("quma_cui"))
    command = [
        "perl",
        f"{quma_path}/quma.pl",
        "-g",
        f"{directory}/{genomic_seq_file}",
        "-q",
        f"{directory}/{reads_seq_file}",
    ]

    out = subprocess.run(command, text=True, capture_output=True).stdout
    return out


def make_list_of_bam_files():
    """Check analysis directory for all valid .bam files.

    Returns:
        list[str]: list of files
    """

    indv_bam = [b for b in config.analysis_directory.iterdir()]

    init_files_set = [bam.name for bam in indv_bam if bam.suffix == ".bam"]

    bai_set = [bai.name for bai in indv_bam if bai.suffix == ".bai"]

    f_set = [str(bam).split(".")[0] for bam in init_files_set]

    baii = [str(bai).split(".")[0] for bai in bai_set]

    files_set = [file + (".TERT.bam") for file in f_set if file in baii]

    final_file_sets = []

    # Code truncates to aleviate excel file_name limit of <=31 characters:
    final_file_sets = [data[:29] for data in files_set]

    return final_file_sets


def index_and_fetch(files_set):
    """Wrapper to call processing of each sam file.

    Args:
        files_set (list[str]): list of bam/sam files

    Returns:
        list[str]: list of genomic positions analyzed
    """

    sam_path = [config.base_directory / "test" / f for f in files_set]

    all_pos = set()
    for sams in sam_path:
        pos = run_sam_and_extract_df(sams)
        all_pos.update(pos)

    return sorted(list(all_pos))


def run_sam_and_extract_df(sams):
    """Process samfiles, pulling out sequence and position data and writing to folders/files.

    Args:
        sams (str): path to a samfile

    Returns:
        list: list of unique positions in the samfile
    """

    # Make sure each sam file has an index by calling external samtools index function
    _ = samtools_index(sams)  # we don't care what the output is

    # Grab the promoter region of interest
    samm = pysam.AlignmentFile(sams, "rb")
    itern = samm.fetch(
        config.chromosome, int(config.promoter_start), int(config.promoter_end)
    )

    position = []
    sequence = []

    for x in itern:
        cols = str(x).split()
        position.append(cols[3])
        sequence.append(cols[9])

    # Transfer into dataframe for processing
    df = pd.DataFrame(list(zip(position, sequence)), columns=["positions", "sequence"])

    df2 = df.set_index("positions")
    # will set the inital index (on the leftmost column) to be position
    df3 = df2.stack()
    # if confused, see: https://www.w3resource.com/pandas/dataframe/dataframe-stack.php

    # Now, iterate through the dataframe
    for each1 in df2.index.unique():
        alignments = []

        # Set up query using alignment algorithm
        query_sequence = df3.loc[each1].head(1).tolist()[0]
        query = StripedSmithWaterman(query_sequence)

        # Set up sequences to check for alignment
        target_sequences = df3.loc[each1].tolist()
        for target_sequence in target_sequences:
            alignment = query(target_sequence)
            alignments.append(alignment)

        read_file = []
        for index, each in enumerate(alignments):
            read_file.append(str(">read" + str(index)))
            read_file.append(alignments[0].aligned_target_sequence)
            # returns aligned target sequence

        # Make sure bam_output directory and sam subdirectories exist
        config.base_directory.joinpath("bam_output", sams.name).mkdir(
            parents=True, exist_ok=True
        )
        directory = config.base_directory.joinpath("bam_output", sams.name)

        with open(directory.joinpath(str(each1) + ".txt"), "w") as file_handler:
            for item in read_file:
                file_handler.write("{}\n".format(item))

    return df2.index.unique()  # return all unique positions in the data


def samtools_index(sams):
    """Helper function to run external samtools index tool.

    Args:
        sams (str): filepath to samfile

    Returns:
        str: output from samtools index shell command, usually discarded
    """

    command = ["samtools", "index", os.fspath(sams)]

    out = subprocess.run(command, capture_output=True, text=True)  # remove shell=True,
    return out


def genome_parsing():
    """Writes out a list of genomic sequence strings for comparison to read data."""

    # Grab list of directories
    subfolders = [x for x in config.bam_directory.iterdir() if x.is_dir()]

    # Grab genomic sequence
    with open(config.promoter_file, "r") as f:
        genome_base = f.readlines()
        genome_base_lines = [s.rstrip("\n") for s in genome_base]
        genome_string = "".join(map(str, genome_base_lines))

    # Wrap everything in processing them one at a time
    for folder in subfolders:

        # Grab list of read files in that directory:
        raw_read_files = os.listdir(folder)
        read_files = [
            os.path.splitext(i)[0]
            for i in raw_read_files
            if not i.lstrip().startswith("g")
            and not i.lstrip().startswith(".ip")
            and not i.lstrip().startswith(".DS")
        ]

        # Now, process each file:
        for read_name in read_files:
            file_lines = []
            # Grab the genomic sequence and write it
            file_lines.append(str(">genome" + str(read_name)))
            file_lines.append(str(genome_range(read_name, genome_string)))

            # Save the reads as a text file for each position
            with open(
                folder.joinpath("g_" + str(read_name) + ".txt"), "w"
            ) as file_handler:
                for item in file_lines:
                    file_handler.write("{}\n".format(item))


def quma_full(cell_types, filename):
    """Run external QUMA methylation analysis on all specified cell lines.

    Args:
        cell_types (list[str]): list of cell lines in our dataset
        filename (str): desired output filename for xlsx output
    """

    # Grab list of directories
    subfolders = [f.path for f in os.scandir(config.bam_directory) if f.is_dir()]

    writer = pd.ExcelWriter(config.base_directory.joinpath(filename))

    # Wrap everything in processing them one at a time
    for folder in tqdm(subfolders, desc="Cell Lines"):

        if any(substring in folder for substring in cell_types):
            # Set up a holding data frame from all the data
            holding_df = pd.DataFrame()

            # Grab list of read files in that directory:
            raw_read_files = os.listdir(folder)
            read_files = [
                os.path.splitext(i)[0]
                for i in raw_read_files
                if i.endswith(".txt") and not i.lstrip().startswith("g")
            ]

            # Now, process each file:
            for read_name in tqdm(read_files, desc="Positions"):
                # file_lines = []

                quma_result = run_quma(
                    folder, "g_{}.txt".format(read_name), "{}.txt".format(read_name)
                )

                dots = []
                for line in quma_result.splitlines():
                    if not line.lstrip().startswith("g"):
                        fields = line.split("\t")
                        if float(fields[7]) < 80:
                            dots.append("FAIL")
                        else:
                            dots.append(fields[13])
                # Sort dots output by number of "1"
                dots2 = sorted(dots, key=lambda t: t.count("1"))

                # Next, add this readname to the holding data frame
                int_df = pd.DataFrame({read_name: dots2})
                holding_df = pd.concat([holding_df, int_df], axis=1)
                # holdng_df[read_name]=dots

            # Now, save it to an excel file
            holding_df.to_excel(writer, sheet_name=Path(folder).name)

            del holding_df

    writer.save()


def extract_cell_types(file_sets):
    """Returns a list[str] of cell lines in the dataset."""

    return [file.split("_")[1] for file in file_sets]


def run_quma_and_compile_list_of_df(cell_types, filename, run_quma=True):
    """Wrapper to run QUMA on all cell lines in the dataset and write output files.

    Args:
        cell_types (list[str]): list of cell lines in the dataset
        filename (str): desired output filename
        run_quma (bool): whether we should invoke the external quma run, default True

    Returns:
        Dict[pd.DataFrame]: dict of dataframes of quma results
    """

    if run_quma:
        quma_full(cell_types, filename)

    df = pd.read_excel(
        config.base_directory.joinpath(filename),
        dtype=str,
        sheet_name=None,
        index_col=0,
    )

    return df


def process_means(dict_of_dfs, positions, cell_types):
    """Process the mean values at each position for each cell line.

    Args:
        dict_of_dfs (Dict[pd.DataFrame]): dict of dataframes of quma results
        positions (list[str]): list of genomic positions to analyze
        cell_types (list[str]): list of cell lines in the dataset

    Returns:
        pd.DataFrame: dataframes of mean values for each position in each cell line
    """

    bad_values = ["N", "F"]  # for interpreting quma returns

    df_full = dict_of_dfs

    # Gives the means of each individual positions-- NOT the mean of the entire dataframe!
    working_df = pd.DataFrame()
    for pos in positions:
        working_df[pos] = ""
        for key, each in df_full.items():
            values_list = []
            if pos in df_full[key].columns:
                if not (
                    len(df_full[key].loc[:, pos].dropna().astype(str)) < 5
                    and not len(df_full[key].loc[:, pos].dropna().astype(str)[0]) < 3
                ):
                    for value in df_full[key].loc[:, pos].dropna().astype(str):
                        if not any(substring in value for substring in bad_values):
                            fraction_val = float(value.count("1")) / float(len(value))
                            values_list.append(fraction_val)

            if values_list:
                pos_means = np.mean(values_list)
            else:
                pos_means = np.nan

            working_df.loc[key, pos] = pos_means

    means_df = working_df

    return means_df


def process_modes(dict_of_dfs, positions, cell_types):
    """Process the mode values at each position for each cell line.

    Args:
        dict_of_dfs (Dict[pd.DataFrame]): dict of dataframes of quma results
        positions (list[str]): list of genomic positions to analyze
        cell_types (list[str]): list of cell lines in the dataset

    Returns:
        pd.DataFrame: dataframes of mode values for each position in each cell line
    """

    bad_values = ["N", "F"]  # for interpreting quma returns

    df_full = dict_of_dfs

    # Gives the modes of each individual positions-- NOT the mean of the entire dataframe!
    working_df = pd.DataFrame()
    for pos in positions:
        working_df[pos] = ""
        for key, each in df_full.items():
            values_list = []
            if pos in df_full[key].columns:
                if not (
                    len(df_full[key].loc[:, pos].dropna().astype(str)) < 5
                    and not len(df_full[key].loc[:, pos].dropna().astype(str)[0]) < 3
                ):
                    for value in df_full[key].loc[:, pos].dropna().astype(str):
                        if not any(substring in value for substring in bad_values):
                            fraction_val = float(value.count("1")) / float(len(value))
                            values_list.append(fraction_val)

            if values_list:
                pos_modes = stats.mode(values_list)[0][0]
            else:
                pos_modes = np.nan

            working_df.loc[key, pos] = pos_modes

    modes_df = working_df

    return modes_df


def return_individual_data(dict_of_dfs, positions, cell_types):
    """Return a dataframe for methylation values at each position for each cell line.

    Args:
        dict_of_dfs (Dict[pd.DataFrame]): dict of dataframes of quma results
        positions (list[str]): list of genomic positions to analyze
        cell_types (list[str]): list of cell lines in the dataset

    Returns:
        pd.DataFrame: dataframes of methylation values for each position in each cell line
    """

    bad_values = ["N", "F"]  # for interpreting quma returns

    df_full = dict_of_dfs

    working_df = pd.DataFrame()
    for pos in positions:
        working_df[pos] = ""
        for key, each in df_full.items():
            values_list = []
            if pos in df_full[key].columns:
                if not (
                    len(df_full[key].loc[:, pos].dropna().astype(str)) < 5
                    and not len(df_full[key].loc[:, pos].dropna().astype(str)[0]) < 3
                ):
                    for value in df_full[key].loc[:, pos].dropna().astype(str):
                        if not any(substring in value for substring in bad_values):
                            fraction_val = float(value.count("1")) / float(len(value))
                            values_list.append(fraction_val)

            if values_list:
                data_for_df = values_list
            else:
                data_for_df = np.nan

            working_df.loc[key, pos] = data_for_df

    return working_df


def find_diffs(means_df, modes_df):
    """Find the differences between means and modes for each cell line at each position
       in means and modes data.

    Args:
        means_df (pd.DataFrame): dataframe of means values
        modes_df (pd.DataFrame): dataframe of modes values

    Returns:
        pd.DataFrame: dataframe of difference values
    """

    return means_df.subtract(modes_df)


def write_means_modes_diffs(means_df, modes_df, diff_df, filename):
    """Wite out files of means, modes, and diffs for future analysis.

    Args:
        means_df (pd.DataFrame): dataframe of means values
        modes_df (pd.DataFrame): dataframe of modes values
        diff_df (pd.DataFrame): dataframe of diff values
        filename (str): desired root filename
    """

    means_df.to_excel(
        config.base_directory / str(filename + "_means.xlsx"),
        sheet_name="Means",
        index=True,
    )
    modes_df.to_excel(
        config.base_directory / str(filename + "_modes.xlsx"),
        sheet_name="Modes",
        index=True,
    )
    diff_df.to_excel(
        config.base_directory / str(filename + "_diff.xlsx"),
        sheet_name="Diff",
        index=True,
    )


def histogram(data_df, cell_line, pos):
    """Return a histogram of the methylation data at a given cell line and position,
       using built-in pandas histogram capacilities.

    Args:
        data_df (pd.DataFrame): a dataframe of individual methylation data
        cell_line (str): Name of a cell line in our data_set
        pos (str): genomic position

    Returns:
        histogram: matplotlib histogram figure"""

    data_to_graph = data_df.loc[cell_line, pos]

    return pd.DataFrame(data_to_graph).hist()


if __name__ == "__main__":
    main(sys.argv[1])  # run the whole shebang using defaults and env variables!
