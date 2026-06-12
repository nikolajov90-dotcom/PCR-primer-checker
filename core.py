from Bio.Seq import Seq
from Bio import SeqIO
from io import TextIOWrapper

def load_genomes(fasta_file):

    # Streamlit UploadedFile
    if hasattr(fasta_file, "read"):

        fasta_file.seek(0)

        handle = TextIOWrapper(
            fasta_file,
            encoding="utf-8"
        )

    # običan filepath
    else:

        handle = fasta_file

    return {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(
            handle,
            "fasta"
        )
    }


def reverse_complement(seq: str) -> str:

    return str(
        Seq(seq).reverse_complement()
    )


def hamming_distance(
        seq1: str,
        seq2: str) -> int:

    distance = 0

    for a, b in zip(seq1, seq2):

        if a != b:
            distance += 1

    return distance


def find_binding_sites(
        genome,
        primer,
        max_mismatches=0):

    primer_hits = []

    primer_length = len(primer)

    # Brzi režim (perfect match)
    if max_mismatches == 0:

        pos = genome.find(primer)

        while pos != -1:

            primer_hits.append(
                (pos, 0)
            )

            pos = genome.find(
                primer,
                pos + 1
            )

        return primer_hits

    # Režim sa mismatch-evima
    for i in range(
            len(genome)
            - primer_length
            + 1):

        if genome[i] != primer[0]:
            continue

        window = genome[
            i:i + primer_length
        ]

        mismatches = hamming_distance(
            window,
            primer
        )

        if mismatches <= max_mismatches:

            primer_hits.append(
                (i, mismatches)
            )

    return primer_hits


def find_primer_hits(
        genomes,
        primer,
        max_mismatches=0):

    hits = {}

    for chromosome, genome in genomes.items():

        hits[chromosome] = (
            find_binding_sites(
                genome,
                primer,
                max_mismatches
            )
        )

    return hits


def find_amplicons(
        forward_hits,
        reverse_hits,
        reverse_primer_length,
        min_size=50,
        max_size=5000):

    amplicons = {}

    for chromosome in forward_hits:

        chromosome_amplicons = []

        f_hits = forward_hits[
            chromosome
        ]

        r_hits = reverse_hits.get(
            chromosome,
            []
        )

        for f_pos, f_mm in f_hits:

            for r_pos, r_mm in r_hits:

                # određivanje orijentacije

                if r_pos > f_pos:

                    strand = "+"

                elif f_pos > r_pos:

                    strand = "-"

                else:
                    continue

                size = (
                    abs(r_pos - f_pos)
                    + reverse_primer_length
                )

                if (
                        min_size
                        <= size
                        <= max_size):

                    chromosome_amplicons.append(
                        (
                            f_pos,
                            r_pos,
                            size,
                            strand,
                            f_mm,
                            r_mm,
                            f_mm + r_mm
                        )
                    )

        if chromosome_amplicons:

            amplicons[
                chromosome
            ] = chromosome_amplicons

    return amplicons


def filter_amplicons(
        amplicons,
        max_length):

    filtered = {}

    for chromosome, products in amplicons.items():

        valid_products = []

        for product in products:

            size = product[2]

            if size <= max_length:

                valid_products.append(
                    product
                )

        if valid_products:

            filtered[
                chromosome
            ] = valid_products

    return filtered


def get_amplicon_sequence(
        genomes,
        chromosome,
        start,
        end,
        reverse_primer_length):

    left = min(start, end)

    right = (
        max(start, end)
        + reverse_primer_length
    )

    return genomes[chromosome][
           left:right
           ]


def analyze_primers(
        fasta_file,
        forward,
        reverse,
        max_mismatches=0,
        max_length=1000,
        min_amplicon_size=50,
        max_amplicon_size=5000):

    forward = forward.upper()
    reverse = reverse.upper()

    genomes = load_genomes(
        fasta_file
    )

    reverse_rc = (
        reverse_complement(
            reverse
        )
    )

    forward_hits = (
        find_primer_hits(
            genomes,
            forward,
            max_mismatches
        )
    )

    reverse_hits = (
        find_primer_hits(
            genomes,
            reverse_rc,
            max_mismatches
        )
    )

    amplicons = (
        find_amplicons(
            forward_hits,
            reverse_hits,
            len(reverse),
            min_amplicon_size,
            max_amplicon_size
        )
    )

    amplicons = (
        filter_amplicons(
            amplicons,
            max_length
        )
    )


    return (
    genomes,
    amplicons,
)
