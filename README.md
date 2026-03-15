# rna-protein-translator
# 🧬 RNA → Protein Translator

A Python tool that translates RNA sequences into proteins.

## Files
| File | Description |
|------|-------------|
| `rna_translator.py` | RNA → Protein translation script |

## Features
- Finds AUG start codon automatically
- Full codon-by-codon breakdown
- Protein in both 3-letter and 1-letter format
- Auto-converts DNA to RNA if needed
- Amino acid frequency summary
- FASTA file support

## How to Run
```
python3 rna_translator.py
```

## Example Output
```
Start codon (AUG) found at position : 1
Codons translated                   : 17
Stop codon reached                  : Yes ✅

Protein : Met - Ala - Met - Ala - Pro - Arg...
1-letter: MAMAPRTEINSTRING
```

## Requirements
Python 3.x — no extra libraries needed!
