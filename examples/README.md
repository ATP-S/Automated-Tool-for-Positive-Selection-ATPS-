# ATPS Example Data

This folder contains example sequence files for testing the ATPS pipeline.

## Files

- `species_list.txt` - List of species to analyze
- `gene_list.txt` - List of genes to analyze
- `BRCA1/` - Example BRCA1 gene sequences
  - `CodingSequences.fasta` - Coding DNA sequences
  - `ProteinSequences.fasta` - Protein sequences
- `TP53/` - Example TP53 gene sequences
  - `CodingSequences.fasta` - Coding DNA sequences
  - `ProteinSequences.fasta` - Protein sequences

## Running the Example

```bash
# From the project root directory
python -m ATPS.main -G BRCA1,TP53 -IF examples -O output -I "Homo sapiens"
```

## Notes

- These are simplified example sequences for testing purposes
- Real analyses should use full-length sequences from NCBI
- The pipeline expects FASTA headers to contain species names
