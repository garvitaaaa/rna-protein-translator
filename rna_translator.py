"""
RNA → Protein Translator
Translates an RNA sequence into amino acids.
Supports plain RNA input and FASTA file input.
"""

CODON_TABLE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

ONE_LETTER = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
    "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Stop": "*",
}


def read_fasta(file_path):
    """Read a FASTA file and return the sequence as a string."""
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()


def validate_rna(seq):
    """Validate that the sequence contains only valid RNA bases."""
    seq = seq.upper().strip().replace(" ", "").replace("\n", "")
    # auto-convert DNA to RNA if T is found
    if "T" in seq and "U" not in seq:
        print("  ⚠️  DNA sequence detected — auto-converting T → U")
        seq = seq.replace("T", "U")
    invalid = set(seq) - {"A", "U", "G", "C"}
    if invalid:
        return False, f"Invalid bases found: {', '.join(sorted(invalid))}"
    if len(seq) == 0:
        return False, "Sequence cannot be empty."
    return True, seq


def translate(rna):
    """
    Translate RNA from the first AUG start codon to a stop codon.
    Returns (protein_3letter, protein_1letter, start_pos, codons_list)
    """
    # find the start codon
    start = rna.find("AUG")
    if start == -1:
        return None, None, -1, []

    protein_3 = []
    protein_1 = []
    codons = []

    for i in range(start, len(rna) - 2, 3):
        codon = rna[i:i+3]
        if len(codon) < 3:
            break
        amino = CODON_TABLE.get(codon, "???")
        codons.append((codon, amino))
        if amino == "Stop":
            protein_3.append("Stop")
            protein_1.append("*")
            break
        else:
            protein_3.append(amino)
            protein_1.append(ONE_LETTER.get(amino, "?"))

    return protein_3, protein_1, start + 1, codons


def print_section(title):
    print(f"\n{'─' * 56}")
    print(f"  {title}")
    print(f"{'─' * 56}")


def display_results(rna, protein_3, protein_1, start_pos, codons):
    print("\n" + "═" * 58)
    print("       🧬  RNA → PROTEIN TRANSLATOR  🧬")
    print("═" * 58)

    print_section("INPUT RNA SEQUENCE")
    print(f"  Length   : {len(rna)} bases")
    print(f"  Sequence : {rna[:60]}{'...' if len(rna) > 60 else ''}")

    if protein_3 is None:
        print("\n  ❌ No AUG start codon found. Cannot translate.")
        return

    print_section("TRANSLATION DETAILS")
    print(f"  Start codon (AUG) found at position : {start_pos}")
    print(f"  Codons translated                   : {len(codons)}")
    has_stop = protein_3[-1] == "Stop" if protein_3 else False
    print(f"  Stop codon reached                  : {'Yes ✅' if has_stop else 'No ❌ (incomplete)'}")

    print_section("CODON TABLE")
    print(f"  {'Codon':<8} {'Amino Acid':<12} {'1-Letter'}")
    print(f"  {'─'*5:<8} {'─'*10:<12} {'─'*6}")
    for codon, amino in codons:
        one = ONE_LETTER.get(amino, "?")
        marker = " ◀ Start" if codon == "AUG" and codons.index((codon, amino)) == 0 else ""
        marker = " ◀ Stop" if amino == "Stop" else marker
        print(f"  {codon:<8} {amino:<12} {one}{marker}")

    print_section("PROTEIN SEQUENCE (3-letter)")
    chunk = 10
    aa_list = [a for a in protein_3 if a != "Stop"]
    for i in range(0, len(aa_list), chunk):
        group = aa_list[i:i+chunk]
        print(f"  {i+1:>4}  {' - '.join(group)}")

    print_section("PROTEIN SEQUENCE (1-letter)")
    one_str = "".join(protein_1).replace("*", "")
    for i in range(0, len(one_str), 60):
        print(f"  {i+1:>4}  {one_str[i:i+60]}")

    print_section("SUMMARY")
    print(f"  Total amino acids : {len(aa_list)}")
    print(f"  Protein length    : {len(aa_list)} aa")

    # amino acid frequency
    freq = {}
    for aa in aa_list:
        freq[aa] = freq.get(aa, 0) + 1
    top5 = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:5]
    print(f"  Most frequent AAs : {', '.join(f'{aa}({n})' for aa, n in top5)}")

    print("\n" + "═" * 58 + "\n")


def run_from_sequence(rna_sequence):
    """Translate a plain RNA or DNA string."""
    valid, rna = validate_rna(rna_sequence)
    if not valid:
        print(f"\n  ❌ Error: {rna}")
        return
    protein_3, protein_1, start_pos, codons = translate(rna)
    display_results(rna, protein_3, protein_1, start_pos, codons)


def run_from_fasta(file_path):
    """Read a FASTA file and translate its sequence."""
    sequence = read_fasta(file_path)
    run_from_sequence(sequence)


# ── Run ────────────────────────────────────────────────────
if __name__ == "__main__":

    # Option 1: translate a plain RNA sequence directly
    run_from_sequence("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")

    # Option 2: translate from a FASTA file — uncomment to use
    # run_from_fasta("fshb_sequence.fasta")