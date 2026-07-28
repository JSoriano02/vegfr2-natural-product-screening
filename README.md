# Nature-Inspired VEGFR-2 Inhibitor Discovery Pipeline

A Nextflow pipeline for large-scale virtual screening of natural compounds (LOTUS database) against the human Vascular Endothelial Growth Factor Receptor 2 (**VEGFR-2**, PDB [3VO3](https://www.rcsb.org/structure/3VO3)). Candidates are filtered by druglikeness (Lipinski) and predicted ADMET safety, prepared in 3D, and docked with **GNINA** (an AutoDock Vina-derived engine with CNN rescoring). Top-ranked poses are then inspected visually and carried forward, outside this repository, into molecular dynamics (MD) system preparation.

> **Scope note:** this git repository implements the *virtual screening* stage only (data acquisition → filtering → docking → ranking → visual inspection). Two top candidates were subsequently carried into MD system preparation and simulation (CHARMM-GUI + GROMACS), run manually in a **sibling directory (`../md_simulations/`) that is not part of this git repository**. That stage is documented in detail in [MD System Preparation and Simulation](#md-system-preparation-and-simulation-external-not-scripted-in-this-repository) from the artifacts and logs found there, but it is not reproducible by running anything in *this* repo — it required manual steps through the CHARMM-GUI web interface. TODO: decide whether `md_simulations/` should be published alongside this repo (e.g. as a companion repo or data archive) so reviewers can access the underlying trajectories/inputs.

## Citation

If you use this pipeline, please cite:

```bibtex
@article{TODO_author_year,
  title   = {TODO: Paper title},
  author  = {TODO: Author list},
  journal = {TODO: Journal},
  year    = {TODO},
  doi     = {TODO}
}
```

<!-- TODO: replace with the final BibTeX entry / DOI once the paper is published or a preprint is available. -->

## Requirements and Installation

### Core requirements

| Tool | Version used in this repo | Notes |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | not pinned (TODO) | Workflow manager; requires Java. `conda.enabled = true` in `nextflow.config` — Nextflow builds one Conda environment per process automatically. |
| Conda / Miniconda | not pinned (TODO) | Required for automatic per-process environment creation (no manual env activation needed). |
| [GNINA](https://github.com/gnina/gnina) | **v1.0.3** | Downloaded automatically as a prebuilt Linux binary at runtime by `modules/docking_gnina.nf` (`gnina/gnina` GitHub release), not installed via Conda/pip. Combines a Vina-based scoring function with CNN pose rescoring. |
| [MGLTools](https://ccsb.scripps.edu/mgltools/) | not pinned (TODO), installed via `bioconda::mgltools` | Provides `prepare_receptor4.py` for receptor PDBQT preparation. |
| [RDKit](https://www.rdkit.org/) | not pinned (TODO), installed via `conda-forge::rdkit` | Lipinski filtering, 3D embedding (ETKDGv3), MMFF94 minimization, SMILES recovery from docked poses. |
| [Meeko](https://github.com/forlilab/Meeko) | not pinned (TODO), installed via `conda-forge::meeko` | Ligand PDBQT writing. |
| [OpenBabel](https://openbabel.org/) | not pinned (TODO), installed via `conda-forge::openbabel` | Used only for pH-based protonation (`obabel -p 7.4`) before 3D embedding. |
| [ADMET-AI](https://github.com/swansonk14/admet_ai) | not pinned (TODO), installed via `pip install admet-ai` at process runtime | ADMET property prediction. |
| GPU (CUDA) | — | Not strictly required, but recommended for `FILTER_ADMET` and `DOCKING_GNINA` (both labeled `gpu_intensive`, `maxForks = 1`). |

Nextflow does **not** use a single top-level environment file for the whole pipeline: each process in `modules/*.nf` declares its own inline `conda "..."` package spec, and Nextflow builds/caches a separate Conda environment per process on first run.

`envs/admet_env.yml` exists in the repo but is **not referenced by `main.nf`, `nextflow.config`, or any module** — it appears to be a standalone environment definition for running/testing `admet-ai` manually, not part of the automated pipeline. TODO: confirm its intended use or remove it if stale.

### Setting up

```bash
# 1. Install Nextflow (requires Java 11+)
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# 2. Install Conda/Miniconda if not already available
# https://docs.conda.io/en/latest/miniconda.html

# No manual `conda env create` step is required for the docking/filtering
# processes: Nextflow creates them automatically on first `nextflow run`
# because `conda.enabled = true` is set in nextflow.config.
```

### External dependencies not installable via pip/conda alone

- **GNINA**: fetched automatically at runtime by `DOCKING_GNINA` via `wget` from the pinned GitHub release (`v1.0.3`). Requires internet access at pipeline run time; no manual install needed, but also not cached outside `work/`.
- **MGLTools**: installed via the `bioconda` channel (`bioconda::mgltools`); no manual download needed given `conda.enabled = true`.
- **CHARMM-GUI** (v3.7, <https://www.charmm-gui.org/>): used manually for the downstream MD system-building step, via its **Solution Builder** / **Ligand Reader & Modeler** modules and its GROMACS **FF-Converter**. No API/CLI integration exists in this repo — inputs were uploaded and downloaded by hand. See [MD System Preparation and Simulation](#md-system-preparation-and-simulation-external-not-scripted-in-this-repository).
- **CGenFF**: ligand force-field parameters were generated through CHARMM-GUI's Ligand Reader & Modeler (which wraps CGenFF) for each ligand, producing `UNL.itp`. The exact CGenFF version and atom-typing/penalty-score log are not exported/committed anywhere — TODO if the paper needs to report them.
- **GROMACS**: CHARMM-GUI exported input files are declared compatible with **GROMACS ≥ 2019.2**; downstream trajectory analysis (RMSD/RMSF/RoG/SASA/H-bonds/PCA) was run with **GROMACS 2026.1**. Force field: **CHARMM36** (via CHARMM-GUI), TIP3P water. Not installed by pip/conda in this repo — install separately (<https://manual.gromacs.org/current/install-guide/index.html>).
- **gmx_MMPBSA** (+ AmberTools, used only as its internal backend): used for MM-GBSA binding free-energy estimates on the MD trajectories. Installed separately from this pipeline's environments, e.g. `conda create -n gmxMMPBSA python=3.10 -c conda-forge && conda install -c conda-forge ambertools=23 && pip install gmx_MMPBSA` (as documented in `../md_simulations/04b_run_mmpbsa.sh`).

## Repository Structure

```
.
├── main.nf                  # Workflow entry point; wires all processes together
├── nextflow.config          # Executor, resource, and default parameter settings
├── modules/
│   ├── obtain_data.nf        # OBTAIN_DATA_RAW: downloads LOTUS SMILES + receptor PDB
│   ├── process_ligands.nf    # CHUNK_LOTUS, FILTER_RDKIT: batching + Lipinski filter
│   ├── admet_filter.nf       # FILTER_ADMET: ADMET-AI safety filtering
│   ├── prepare_meeko.nf      # PREPARE_MEEKO_3D: SMILES -> 3D -> PDBQT ligand prep
│   ├── prepare_receptor.nf   # PREPARE_RECEPTOR: receptor -> PDBQT via MGLTools
│   ├── docking_gnina.nf      # DOCKING_GNINA: GNINA docking + CNN rescoring
│   └── paste_results.nf      # CONSOLIDATE_RESULTS: ranking, SMILES recovery, master CSV
├── envs/
│   └── admet_env.yml         # Standalone Conda env for admet-ai (not wired into the pipeline — TODO)
├── notebooks/
│   └── 01_VGFR2_Top_Candidates_Visualization.ipynb  # Phase 2: 2D grid + 3D pocket inspection
├── raw_data/                 # Pipeline downloads (gitignored): lotus_full.smi, 3VO3.pdb
├── results/
│   └── top_final/
│       ├── top_candidates_master.csv  # Ranked top candidates with ADMET + docking scores
│       └── top_poses/*.sdf            # 3D poses of top candidates (gitignored intermediates aside)
├── images/
│   └── Flux_diagram.png      # Workflow diagram referenced below
└── LICENCE                   # MIT License
```

`raw_data/` and `results/` are listed in `.gitignore`; the copies present in a working checkout are pipeline outputs from a prior local run, not committed example data.

## Pipeline

Diagram of the automated stages:

<img src="images/Flux_diagram.png" alt="Workflow Diagram" width="50%">

Run the full pipeline with:

```bash
nextflow run main.nf -resume
```

Key parameters (`nextflow.config`, overridable with `--param value` or by editing the file):

```groovy
params.lotus_url  = "https://lotus.naturalproducts.net/download/smiles"
params.pdb_id      = "3VO3"
params.chunk_size  = 10000
params.raw_dir     = "${projectDir}/../raw_data"
params.outdir      = "${projectDir}/../results"
```

> **TODO (path check):** `raw_dir`/`outdir` are defined relative to `${projectDir}/..` (one level above the directory containing `main.nf`), but a working checkout shows `raw_data/` and `results/` populated *inside* the repository root (sibling to `main.nf`, not one level above it). Also, `notebooks/01_VGFR2_Top_Candidates_Visualization.ipynb` hardcodes `../nextflow_pipeline/results/...`, which does not match this repository's current flat layout. Verify the intended directory layout (a nested `nextflow_pipeline/` subfolder may have existed before a prior refactor) and correct either the config, the notebook paths, or this note before relying on it for reproduction.

### Step-by-step

**1. Data acquisition — `OBTAIN_DATA_RAW`**
Downloads the full LOTUS natural-product SMILES database and the VEGFR-2 receptor structure from the PDB.
- Command (as run inside the process): `wget -O lotus_full.smi <lotus_url>`, `wget -O 3VO3.pdb https://files.rcsb.org/download/3VO3.pdb`
- Input: `params.lotus_url`, `params.pdb_id`
- Output: `lotus_full.smi` (tab-separated `SMILES<TAB>ID`, 276,518 entries as of the last run), `3VO3.pdb`

**2. Batching and druglikeness filter — `CHUNK_LOTUS`, `FILTER_RDKIT`**
Splits the SMILES database into chunks and applies Lipinski's Rule of Five with RDKit.
- Command: `split -l 10000 lotus_full.smi batch_` then, per batch, a Python/RDKit script filtering on `MolWt <= 500`, `MolLogP <= 5`, `NumHDonors <= 5`, `NumHAcceptors <= 10`, `NumRotatableBonds < 10`
- Input: `lotus_full.smi`
- Output: `batch_*_viable.smi` (SMILES + ID, canonicalized by RDKit)

**3. ADMET safety filter — `FILTER_ADMET`**
Predicts ADMET properties with ADMET-AI and applies clinical-safety thresholds tuned for intravitreal/ophthalmic use.
- Command: `admet_predict --data_path input_with_header.csv --save_path predictions.csv`, followed by an in-process filter requiring `Eye_Irritation`, `Neurotoxicity`, `Immunotoxicity`, `AMES`, `Hepatotoxicity`/`DILI`, and `BBB_Martins` all `< 0.5`
- Input: `batch_*_viable.smi`
- Output: `clinically_safe_batch_*_viable.smi` (filtered SMILES for docking), `admet_data_*.csv` (full ADMET property table, carried through to the final report)

**4. 3D ligand preparation — `PREPARE_MEEKO_3D`**
Converts each safe 2D SMILES into a 3D structure and writes AutoDock-compatible PDBQT files.
- Command (conceptually): `obabel <input> -O protonated.smi -p 7.4`, then per molecule: `Chem.AddHs` → `AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())` → `AllChem.MMFFOptimizeMolecule(mol)` → `Meeko.MoleculePreparation().prepare(mol)` → `write_pdbqt_string()`
- Input: `clinically_safe_*.smi`
- Output: one `.pdbqt` file per successfully embedded ligand

**5. Receptor preparation — `PREPARE_RECEPTOR`**
Splits the crystal structure into protein and co-crystallized ligand, then prepares the receptor with MGLTools.
- Command:
  ```bash
  awk '/^HETATM/ && !/HOH/' 3VO3.pdb > control_ligand.pdb
  awk '/^ATOM/' 3VO3.pdb > receptor_clean.pdb
  prepare_receptor4.py -r receptor_clean.pdb -o ready_receptor.pdbqt -A hydrogens -U nphs_lps_waters
  ```
- Input: `3VO3.pdb`
- Output: `ready_receptor.pdbqt` (polar hydrogens added, Kollman charges via MGLTools defaults), `control_ligand.pdb` (co-crystallized ligand, used to define the docking box)

**6. Docking — `DOCKING_GNINA`**
Docks every ligand PDBQT against the prepared receptor with GNINA, using the co-crystallized ligand to auto-define the search box.
- Command:
  ```bash
  ./gnina -r ready_receptor.pdbqt \
          -l <ligand>.pdbqt \
          --autobox_ligand control_ligand.pdb \
          --exhaustiveness 8 \
          --cnn_scoring rescore \
          --out <ligand>_docked.sdf
  ```
- Input: ligand PDBQTs, `ready_receptor.pdbqt`, `control_ligand.pdb`
- Output: `results_batch_*.csv` (`SMILES_ID, Vina_Affinity_kcal_mol, CNN_Pose_Score, CNN_Affinity_pKd` for the top pose of every ligand attempted); `successful_poses/*.sdf` — poses kept only if `Vina_Affinity <= -7.0 kcal/mol` **and** `CNN_Pose_Score >= 0.5` (all other poses/logs are deleted to save disk space)

**7. Consolidation and ranking — `CONSOLIDATE_RESULTS`**
Aggregates all batch results, applies a stricter final cutoff, recovers SMILES from the 3D poses, and merges in ADMET data.
- Command: in-process Python/RDKit/pandas script (no external CLI)
- Filter: `Vina_Affinity <= -8.0 kcal/mol` **and** `CNN_Pose_Score >= 0.6`, sorted by affinity, top 50 kept
- Input: all `results_batch_*.csv`, all `successful_poses/*.sdf`, all `admet_data_*.csv`
- Output: `results/top_final/top_candidates_master.csv` (ranked candidates with docking scores, ADMET fields `BBB_Martins`/`AMES_Toxicity`/`DILI_Risk`/`Carcinogens`, and `Extracted_SMILES` recovered from the docked 3D pose), `results/top_final/top_poses/*.sdf`

**8. Visual inspection — `notebooks/01_VGFR2_Top_Candidates_Visualization.ipynb`** *(manual, not part of `nextflow run main.nf`)*
Loads `top_candidates_master.csv`, renders a 2D grid of the top 9 candidates with RDKit, and overlays the best pose on the receptor in 3D with py3Dmol.
- Dependencies: `pandas`, `rdkit`, `py3Dmol` (installed inline in the notebook via `!pip install`)
- Input: `results/top_final/top_candidates_master.csv`, `raw_data/3VO3.pdb`, `results/top_final/top_poses/*.sdf`
- Output: inline 2D/3D visualizations (not saved to disk by default)

### MD System Preparation and Simulation (external, not scripted in this repository)

Two docking candidates (**LTS0070549**, **LTS0174212**) plus a structural reference (**Imidazol**, the heterocyclic fragment co-crystallized in 3VO3) were carried forward into MD. This stage lives entirely in a sibling directory, `../md_simulations/` (outside this git repository — not committed, not tracked), and is **manual/semi-scripted, not runnable via `nextflow run main.nf`**. It is documented here, transcribed from the shell scripts, `.mdp` files, and `resumen_ejecutivo_VEGFR2.md` found in that directory, so that the full paper methodology is available in one place. TODO: decide whether `md_simulations/` should be published as supplementary data (its own repo/archive) so this protocol is independently reproducible by reviewers.

**System preparation (per system, before CHARMM-GUI)**
- Receptor: VEGFR-2 kinase domain (residues 812–1168), chain A of `3VO3.pdb`.
- Ligand pose: for LTS0070549/LTS0174212, the GNINA-docked pose (`*_docked.sdf`, output of this repo's `DOCKING_GNINA`) is remapped onto the original co-crystallized ligand template with RDKit (`fix_ligant.py`): maximum-common-substructure atom matching (`rdFMCS`, `BondCompare.CompareAny`) transplants the docked heavy-atom coordinates onto the reference scaffold, then explicit hydrogens are rebuilt (`Chem.AddHs(..., addCoords=True)`) before upload.

**System building (CHARMM-GUI v3.7, <https://www.charmm-gui.org/>)**
- Modules: **Solution Builder** (protein + ligand + water/ion box) and **Ligand Reader & Modeler** (ligand parametrization via CGenFF), with the GROMACS **FF-Converter** for output.
- Force field: **CHARMM36**. Water model: **TIP3P** (`toppar/TIP3.itp`). Ions: K⁺/Cl⁻ (`POT`/`CLA`), consistent with CHARMM-GUI's default neutralization + physiological salt (0.15 M, matching `saltcon=0.150` used downstream for MM-GBSA) — exact box padding/shape and the salt-concentration value entered in the wizard are not exported/committed (TODO).
- Ligand residue name: `UNL` in all three systems (index group 13 in every case — 53 atoms Imidazol, 59 LTS0174212, 72 LTS0070549).
- Example composition (`LTS0070549`, replica 1): 1 protein (`PROA`), 1 ligand (`UNL`), 53 K⁺, 55 Cl⁻, 18,936 TIP3P waters. Ion counts differ slightly per system/replica (e.g. `LTS0174212`: 45 K⁺ / 47 Cl⁻) as each was solvated/ionized independently by CHARMM-GUI.
- Exported files per system: `step3_input.{gro,pdb,psf}`, `topol.top`, `toppar/{forcefield,PROA,UNL,POT,CLA,TIP3}.itp`, `index.ndx`.

**Simulation protocol (GROMACS, CHARMM-GUI-generated `.mdp` files)**

| Stage | Ensemble | Length | Key settings |
|---|---|---|---|
| Minimization | — | 5,000 steps steepest descent | `emtol=1000 kJ·mol⁻¹·nm⁻¹`; backbone/side-chain position restraints 400/40 kJ·mol⁻¹·nm⁻² |
| Equilibration | NVT | 125 ps (`dt=0.001`, 125,000 steps) | Same restraints; v-rescale thermostat, 310.15 K, groups SOLU/SOLV, τ=1.0 ps; velocities regenerated (`gen-vel=yes`, `gen-seed=-1`) |
| Production | NPT | 100 ns/replica (`dt=0.002`, 5×10⁷ steps) | Restraints released; v-rescale thermostat (310.15 K) + C-rescale barostat (1.0 bar, isotropic, τ=5.0 ps, compressibility 4.5×10⁻⁵ bar⁻¹) |

Common to all stages: Verlet cutoff scheme, LINCS constraints on h-bonds, PME electrostatics (`rcoulomb=1.2 nm`), van der Waals force-switch (1.0–1.2 nm). Command sequence (`gmx grompp` → `gmx mdrun`, chained minimization → equilibration → production) is in `../md_simulations/LTS0070549/LTS0070549_gromacs.sh`; production used GPU acceleration (`gmx mdrun -nb gpu -pme gpu`).

- **Replicas**: 2 independent replicas per system (velocities re-randomized at each `grompp` call — `gen-seed=-1` — which is how replica independence is achieved; no fixed seed). 3 systems × 2 replicas × 100 ns = **6 trajectories, 600 ns total**.
- **Analysis** (GROMACS 2026.1, run from `../md_simulations/01_run_analyses.sh`): PBC correction via `gmx trjconv -pbc mol -center`, then `gmx rms` (backbone RMSD), `gmx rmsf -res`, `gmx gyrate`, `gmx sasa`, `gmx hbond` (protein vs. ligand `UNL`), and PCA (`gmx covar` + `gmx anaeig` on the protein backbone) with a free-energy landscape `ΔG = −kT·ln(P)` over PC1/PC2. First 20 ns of each 100 ns production trajectory discarded as additional relaxation before computing "stable-block" statistics. Plots produced with Python (numpy, pandas, matplotlib) in `02_plot_replicas.py`/`03_embed_figures.py`.
- **Binding free energy**: `gmx_MMPBSA` (MM-GBSA, `igb=5`, `saltcon=0.150`, CHARMM36 topology fed directly as `-cp topol.top` — no Amber force field specified), 81 frames of the stable block per replica, averaged over the 2 replicas per system, `T=298.15 K` (`../md_simulations/04b_run_mmpbsa.sh`). MM-PBSA (`&pb` block) is prepared but commented out/not run by default.

**Headline results** (from `resumen_ejecutivo_VEGFR2.md`; reported here for methods completeness, not reproducible from this repo):

| System | RMSD (nm) | RoG (nm) | SASA (nm²) | H-bonds | ΔG_bind MM-GBSA (kcal/mol) |
|---|---|---|---|---|---|
| Imidazol (reference fragment) | 0.194 ± 0.007 | 2.037 ± 0.008 | 163.13 ± 1.17 | 2.50 ± 0.02 | −48.93 ± 0.04 |
| LTS0070549 | 0.212 ± 0.005 | 2.044 ± 0.002 | 165.42 ± 0.78 | 0.43 ± 0.21 | −21.16 ± 2.58 |
| LTS0174212 | 0.219 ± 0.004 | 2.053 ± 0.001 | 166.24 ± 0.55 | 2.75 ± 0.15 | −35.76 ± 2.22 |

Values are mean ± SD across the 2 replicas per system, computed over the >20 ns stable block. LTS0174212 sustains protein–ligand H-bonds comparable to the co-crystallized fragment and has the more favorable ΔG_bind of the two candidates; LTS0070549 loses stable H-bond contacts in simulation despite its docking score.

Remaining TODOs for this stage: exact CHARMM-GUI wizard settings not exported (box padding/shape, salt-concentration field, protonation-state assignment method — e.g. PROPKA — if any); CGenFF version and atom-typing/penalty-score log not retained; explicit, code-level criterion for why exactly these two of the pipeline's top-50 candidates were promoted to MD (the summary states "best dockings" qualitatively, without a documented cutoff).

## Minimal Reproducible Example

The repository ships with the outputs of a prior full run rather than a small curated example:
- `raw_data/lotus_full.smi` — 276,518 SMILES entries (the full LOTUS download, not a subset)
- `raw_data/3VO3.pdb` — the VEGFR-2 receptor used
- `results/top_final/top_candidates_master.csv` — 50 ranked top candidates from that run
- `results/top_final/top_poses/*.sdf` — corresponding docked poses

To reproduce end-to-end from scratch:

```bash
nextflow run main.nf -resume
```

TODO: no small/sample SMILES subset is provided for a quick smoke test; running from scratch re-downloads and re-processes the full LOTUS database (`params.chunk_size = 10000` per batch). Consider adding a `--lotus_url` override pointing to a small local `.smi` file for fast reviewer verification.

## Reproducibility Notes

- **Docking box**: defined dynamically per run via `--autobox_ligand control_ligand.pdb` (GNINA auto-boxes around the co-crystallized ligand extracted from `3VO3.pdb`); no fixed center/size coordinates are hardcoded.
- **Search exhaustiveness**: `--exhaustiveness 8` (`modules/docking_gnina.nf`).
- **Scoring**: GNINA default Vina-based scoring function plus `--cnn_scoring rescore` (default CNN ensemble bundled with GNINA v1.0.3); no custom/trained CNN model is specified.
- **Random seed**: TODO — no `--seed` is passed to GNINA/AutoDock Vina anywhere in `modules/docking_gnina.nf`, so docking runs are not seeded and results may vary slightly between reruns.
- **3D embedding**: RDKit `ETKDGv3` with default RDKit seed behavior (also unseeded — TODO if determinism is required), followed by MMFF94 energy minimization (`AllChem.MMFFOptimizeMolecule`, default iteration count).
- **Ligand protonation**: OpenBabel `-p 7.4` (physiological pH) is attempted before 3D embedding; falls back silently to the unprotonated SMILES if OpenBabel fails (`modules/prepare_meeko.nf`).
- **Receptor preparation**: MGLTools `prepare_receptor4.py -A hydrogens -U nphs_lps_waters` — non-polar hydrogens, lone pairs, and waters are removed; polar hydrogens and Kollman charges are added (MGLTools defaults, no custom charge set).
- **Selection thresholds** (manual/scientific decisions baked into the code, not derived from data):
  - Per-batch pose retention: Vina affinity ≤ −7.0 kcal/mol **and** CNN pose score ≥ 0.5 (`modules/docking_gnina.nf`)
  - Final top-candidate cutoff: Vina affinity ≤ −8.0 kcal/mol **and** CNN pose score ≥ 0.6, top 50 by affinity (`modules/paste_results.nf`)
  - ADMET safety cutoff: all six liabilities (`Eye_Irritation`, `Neurotoxicity`, `Immunotoxicity`, `AMES`, `Hepatotoxicity`/`DILI`, `BBB_Martins`) < 0.5 (`modules/admet_filter.nf`)
- **Software pinning gaps** (TODO): GNINA is the only precisely pinned tool (`v1.0.3`, fixed download URL). RDKit, Meeko, OpenBabel, MGLTools, pandas, and admet-ai are installed as unpinned `conda`/`pip` latest-at-install-time packages — for exact reproducibility, pin versions in each module's `conda` directive and record the resolved environment (e.g., via `conda list --export` after a real run).
- **MD stage**: see [MD System Preparation and Simulation](#md-system-preparation-and-simulation-external-not-scripted-in-this-repository) — protocol is fully documented (CHARMM36/TIP3P via CHARMM-GUI, GROMACS minimization/NVT/NPT, 2×100 ns replicas per system, `gmx_MMPBSA` MM-GBSA), but it lives outside this git repository (`../md_simulations/`) and is not runnable from here; replica independence comes from GROMACS regenerating velocities per run (`gen-seed=-1`), not from a fixed seed.

## License

This project is licensed under the MIT License — see [`LICENCE`](LICENCE) for details.

## Contact

Jaime Soriano Ivorra — jsorianoi1402@gmail.com
