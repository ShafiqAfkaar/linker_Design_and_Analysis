# Linker Design and Analysis

**Project website:** https://shafiqafkaar.github.io/linker_Design_and_Analysis/

Python scripts for designing protein linkers with RFdiffusion and ProteinMPNN, validating designs with ColabFold, and analyzing linker geometry and confidence from predicted PDB structures.

## Files

- `design_linker_pipeline.py` — runs RFdiffusion, ProteinMPNN, ColabFold, and the analysis step.
- `analyze_linker.py` — calculates linker sequence, end-to-end distance, mean solvent-accessible surface area (SASA), and mean pLDDT from PDB B-factors.

## Requirements

The analysis script requires Python 3 plus the packages in `requirements.txt`.

```bash
python -m pip install -r requirements.txt
```

The full design pipeline additionally expects working installations of:

- RFdiffusion
- ProteinMPNN
- ColabFold (`colabfold_batch`)

Before running the pipeline, update the tool paths near the top of `LinkerDesignPipeline.__init__` if your installations are not under `/root/autodl-fs`.

## Analyze predicted structures

```bash
python analyze_linker.py \
  --input /path/to/pdb_results \
  --phage_len 100 \
  --pep_len 20 \
  --output ./linker_results
```

Use `--chain A` to select a specific chain. Without it, the first chain is analyzed.

## Run the complete design pipeline

```bash
python design_linker_pipeline.py \
  --pdb_path /absolute/path/to/input.pdb \
  --phage_start 1 \
  --phage_end 100 \
  --pep_start 111 \
  --pep_end 130 \
  --output_dir ./unified_output
```

Run either script with `--help` for all options.

## Notes

- Residue lengths are used to identify the linker: the N-terminal phage length and C-terminal peptide length are excluded.
- pLDDT is read from PDB B-factors, as used in common AlphaFold/ColabFold PDB outputs.
- Generated outputs are ignored by Git through `.gitignore`.
