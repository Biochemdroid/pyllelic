
import Bio
from Bio.Seq import Seq
from typing import List, Dict, Any, Tuple
from collections import namedtuple
import re

def fasta_output(seq: str, seq_name: str, line: int = None) -> str:
    """Write a sequence string to a fasta-formatted text file contents.

    Args:
        seq (str): sequence string
        seq_name (str): sequence name
        line (int, optional): Max line length to process. Defaults to None.

    Returns:
        str: fasta-formatted text file contents.
    """
    # line = ling or MAX_LINE_LENGTH

    seq = re.sub(r"[0-9]| |\t|\n|\r|\f", "", seq)

    reg = r"(.{1," + str(line) + "})"
    seq = re.sub(reg, r"\1\n", seq)

    return f">{seq_name}\n{seq}"




def check_char_in_allowed(seq: str, pattern: str) -> str:
    """Return only charcters in string present in pattern.

    Args:
        seq (str): sequence string
        patt (str): string of allowed characters

    Returns:
        str: string with unallowed characters removed.
    """
    new = ""
    for each in seq:
        if each in pattern:
            new += each
    return new

    
def curate_seq(seq: str) -> str:
    """Curate a sequence to only have allowed characters.

    Args:
        seq (str): sequence to check.

    Returns:
        str: sequence without dis-allowed characters.
    """
    return check_char_in_allowed(seq, "ACGTURYMWSKDHBVNacgturymwskdhbvn")

def parse_seq(seq: str) -> str:
    """Extract sequence strings from the string of a text file.

    Args:
        seq (str): stringified text file.

    Returns:
        str: sequence string
    """
    _ = ""
    seq = re.sub(r"\r\n", "\n", seq)
    seq = re.sub(r"\r", "\n", seq)
    seq = seq.upper()

    reg = re.compile(r"^\s*>.*?\n", re.MULTILINE)
    if re.findall(reg, seq):
        seq = re.sub(r"^\s*>(.*)?\n", "", seq)

    elif re.findall(
        r"^ORIGIN\s*\n((\s+(\d+(\s+\w+)+))+)\s*\n//", seq, re.MULTILINE
    ):  # pragma: no cover
        seq = re.findall(reg, seq, re.MULTILINE)[0]

    elif re.findall(
        r"^SQ\s+SEQUENCE.*\n((\s+(\w+\s+)+\d+)+)\n\/\/", seq, re.MULTILINE
    ):  # pragma: no cover
        seq = re.findall(
            r"^SQ\s+SEQUENCE.*\n((\s+(\w+\s+)+\d+)+)\n\/\/", seq, re.MULTILINE
        )[0]

    elif re.findall(
        r"\.\.\s*\n((\s+(\d+(\s+\w+)+))+)\s*", seq, re.MULTILINE
    ):  # pragma: no cover
        seq = re.findall(r"\.\.\s*\n((\s+(\d+(\s+\w+)+))+)\s*", seq, re.MULTILINE)[0]
    elif re.findall(r"^\s*>.+\s.+", seq, re.MULTILINE):  # pragma: no cover
        seq = re.findall(r"^\s*>(.+?)\s(?=.+)", seq, re.MULTILINE)[0]
        _ = seq

    return curate_seq(seq)


def multi_fasta_parse(multi: Any) -> List[Dict[str, str]]:
    """Find all bisufite sequencing reads in a fasta file,
       and return as a dictionary.

    Args:
        multi (Any): multi-line sequence string

    Returns:
        List[Dict[str, str]]: list of dictionaries of sequence reads
    """
    multi = re.sub(r"\r\n", "\n", multi)
    multi = re.sub(r"\r", "\n", multi)
    biseq: List[Dict[str, str]] = []
    fa: Dict[str, str] = {}

    multi = re.findall(r"(.*)$", multi, re.MULTILINE)

    for line in multi:

        if ">" in line:
            if fa and not fa["seq"]:  # pragma: no cover
                biseq.pop()
            fa = {"com": line}
            fa["com"] = re.sub(r"^>", "", fa["com"])
            fa["com"] = re.sub(r"\s*$", "", fa["com"])
            biseq.append(fa)
        else:
            line = curate_seq(line)
            if line == "":
                continue
            if not fa:  # pragma: no cover
                return None  # does this ever happen?
            try:
                fa["seq"] += line.upper()
            except KeyError:
                fa["seq"] = line.upper()

    if fa:
        if not fa.get("seq"):  # pragma: no cover
            biseq.pop()

    return biseq



def parse_genome(file: str) -> str:
    """Parse genome file, removing white spaces and extra returns.

    Args:
        file (str): text file path for fasta file.

    Returns:
        str: parsed and curated string of genome sequence.
    """
    with open(file, "r") as f:
        seq = f.read()

    seq = re.sub(r"^[\r\s]+", "", seq)
    seq = re.sub(r"[\r\s]+$", "", seq)
    seq = re.sub(r"(\r|\n|\r\n){2}", "\r|\n|\r\n", seq)

    return parse_seq(seq)


def parse_biseq(file: str) -> List[Dict[str, str]]:
    """Parse bisulfite sequencing fasta file.

    Args:
        file (str): file path.

    Returns:
        List[Dict[str, str]]: list of dictionaries of sequence reads
    """
    multi = _multi_parser(file)

    return multi_fasta_parse(multi)

    
def parse_multi(file: str) -> Tuple[None, List[Dict[str, str]]]:
    """Parse fast sequencing file with multiple reads.

    Args:
        file (str): file path

    Returns:
        Tuple[None, List[Dict[str, str]]]: None and list of dicts of sequence reads
    """
    multi = _multi_parser(file)

    biseq = multi_fasta_parse(multi)
    return None, biseq


def _multi_parser(file: str) -> Any:
    """Helper to do substitution in fasta file contents."""
    with open(file, "r") as f:
        multi = f.read()

    multi = re.sub(r"^[\r\s]+", "", multi)
    multi = re.sub(r"[\r\s]+$", "", multi)
    multi = re.sub(r"(\r\n){2}", "\r\n", multi)
    multi = re.sub(r"(\n){2}", "\n", multi)
    multi = re.sub(r"(\r){2}", "\r", multi)

    return multi



def cpg_main(gfile: str, qfile: str) -> str:

    gseq = parse_genome(gfile)
    qseq = parse_biseq(qfile)

    return find_cpg(gseq,qseq)










# define some organizational Named Tuples (allow access by name instead of number)
SeqPair = namedtuple("SeqPair", ["genomic_Seq", "query_Seq"])
SeqAlignments = namedtuple("SeqAlignments", ["Both_Fwd", "G_fwd_Q_rev", "G_rev_Q_fwd"])
CpG_Counts = namedtuple("CpG_Counts", ["counts", "cpg_dict", "quma"])



def generate_seq_combinations(g_seq: Seq, q_seq: Seq) -> SeqAlignments:
    both_fwd = SeqPair(g_seq, q_seq)
    g_fwd_q_rev = SeqPair(g_seq, q_seq.reverse_complement())
    g_rev_q_fwd = SeqPair(g_seq.reverse_complement(), q_seq)
    return SeqAlignments(both_fwd, g_fwd_q_rev, g_rev_q_fwd)



def find_cpgs(seq):
    return [x.start() for x in list(re.finditer("CG", str(seq)))]   

def make_cpg_dict(cpg_locations, seq):
    cpg_dict = {key: None for key in cpg_locations}
    for pos in cpg_locations:
        if str(seq)[pos] == "C":
            cpg_dict[pos] = 1
        elif str(seq)[pos] == "T":
            cpg_dict[pos] = 0
        else:
            cpg_dict[pos] = None

    return cpg_dict

def cpg_dict_add(cpg_dict):
    # From https://stackoverflow.com/questions/47327646/handling-none-when-adding-numbers
    return sum(filter(None, [x for x in cpg_dict.values()]))


def make_quma_repr(cpg_dict):
    return "".join([str(x) if x is not None else "-" for x in cpg_dict.values()])


def find_cpg(all_combos):
    find_best_dict = {}
    for name, seqpair in zip(all_combos._fields, all_combos):
        cpg_locations = find_cpgs(seqpair.genomic_Seq)
        cpg_dict = make_cpg_dict(cpg_locations, seqpair.query_Seq)
        cpg_count = cpg_dict_add(cpg_dict)
        quma = make_quma_repr(cpg_dict)
        find_best_dict[name] = CpG_Counts(cpg_count, cpg_dict, quma)

    best = max(find_best_dict, key=lambda x: find_best_dict[x][0])

    return find_best_dict[best]

#