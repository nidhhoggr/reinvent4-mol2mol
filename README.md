

# TL Example use-case: 16 compounds

Before the mechanics: mol2mol TL doesn't train on your 16 molecules directly — it trains on _pairs_ of them that are similar enough to each other (Tanimoto ≥ `pairs.lower_threshold`). With 16 compounds you have at most 16×15 = 240 ordered pairs, and the threshold knocks that down fast. So the make-or-break question isn't "is 16 enough rows" — it's "how many pairs survive the threshold." That's why the first thing you run is the pair-checker, not the trainer. If your 16 are one tight congeneric series you'll be fine; if they're scattered chemotypes you may get near-zero pairs at the default 0.7 and the run will barely learn. Everything below is built around managing that.

## Step 1 — Convert SDF → SMILES (`scripts/tl/01_sdf_to_smiles.py`)

mol2mol reads plain SMILES, one per line, first column used and a second column ignored. This desalts, neutralizes, canonicalizes, and dedupes (tested: a Na-salt and a duplicate collapsed correctly to one canonical parent):

Run it: `python scripts/tl/01_sdf_to_smiles.py data/raw_sdfs/ -o configs/compounds.smi`

One caveat worth knowing: the mol2mol prior has a fixed token vocabulary (drug-like organics). If a compound has exotic atoms or unusual valences, TL silently drops it. Standardizing here reduces that, but eyeball the output — if 16 went in and 11 came out, find out why before training.

## Step 2 — Check pair counts _before_ training (`scripts/tl/02_check_pairs.py`)

This is the diagnostic that decides your `lower_threshold`. It prints how many pairs survive each threshold:

Run: `python scripts/tl/02_check_pairs.py configs/compounds.smi`

How to read it: you want a threshold where you get a healthy number of pairs _and_ most of your 16 appear as sources (the right-hand column). Pick the **highest** threshold that still gives you a few hundred pairs and near-full source coverage. For a tight series that's often 0.6–0.7; for a looser set you may have to drop to 0.5 or 0.4. If even 0.3 gives you almost nothing, your set is too diverse for mol2mol TL to find structure in — that's a signal to either add compounds or accept a very light bias.

### Example output (Poor Pairing):

```
17 valid molecules -> 272 ordered pairs possible (excl. self) threshold | pairs >= thr | sources w/>=1 target 
------------------------------------------------ 
0.3 | 12 | 7 
0.4 | 8  | 5 
0.5 | 2  | 2 
0.6 | 2  | 2 
0.7 | 0  | 0 
0.8 | 0  | 0
```


Here's how to read it, column by column, and then the verdict — because this table is telling you something important.

**The columns:**

- "ordered pairs possible" — 17 molecules can form 17×16 = 272 directed pairs (A→B counts separately from B→A).
- "pairs >= thr" — how many of those 272 actually clear that Tanimoto cutoff. Because they're ordered, divide by 2 for the number of _distinct_ similar duos: 12 ordered at 0.3 ≈ 6 real pairs; 2 ordered at 0.5/0.6 = exactly **one** real pair.
- "sources w/>=1 target" — how many of your 17 have at least one similar partner. The rest are orphans that contribute nothing to TL.

**The verdict: this set is too diverse for mol2mol TL.** Three things jump out:

At the tutorial's default threshold of 0.7, you have **zero pairs**. With that config the TL run has literally nothing to train on. Same at 0.8. So the defaults are off the table entirely.

At 0.5–0.6 you have one similar pair out of 17 molecules. Training a model on a single pair isn't training — it's memorizing one edit.

Even at 0.3 — which is a very loose cutoff, "vaguely related" on ECFP4, not really analogs — you only scrape ~6 pairs, and **10 of your 17 compounds still have no partner at all** (only 7 sources). Those 10 sit out of TL no matter what threshold you pick.

So the table is saying your 17 are scattered chemotypes, not a congeneric series. mol2mol TL learns similarity-preserving edits _from pairs of mutually similar molecules_, and you don't have them. Lowering the threshold to force pairs is a trap here: you'd be training on a handful of barely-related pairs while two-thirds of your set is ignored, and you'd get an unfocused or degenerate agent for your trouble.

**What to actually do depends on your goal, and there are two clean paths:**

If what you want is _analogs of each of these 17 compounds_, you may not need TL at all. mol2mol is conditional — at sampling time you hand it any seed molecule and it generates analogs around it using the base prior. So you can skip the broken pairing step entirely: sample from `mol2mol_medium_similarity.prior` (or the high-similarity one) seeded with each of your 17 individually. That sidesteps the whole problem and is probably the right move for a diverse hit set.

If what you want is _an agent biased toward this chemical space as a series_, you need congeneric data the set doesn't currently have. That means augmenting (Step 3): for each compound (or each small cluster), pull in known analogs — from a SAR series, a similarity search against ChEMBL/your in-house DB, or R-group enumeration — until clusters form and pairs appear at a respectable threshold (aim for a few hundred pairs at ≥0.5–0.6 with most molecules covered). Then re-run `check_pairs.py` and you'll see the table fill in.


## Step 3 — ChEMBL Augmentation (`scripts/tl/03_augment_from_chembl.py`)

When you need congeneric data that your compound pool doesn't have you need to augment your data using a similarity search. This script used ChEMBL to find similar compounds to de-diversify your data set.

**How it works.** For each seed it hits ChEMBL's `similarity/{smiles}/{percent}` endpoint, standardizes every hit with the same desalt/neutralize/canonicalize logic as your converter, drops anything outside the heavy-atom window, caps how many it keeps per seed, and dedupes against everything (seeds win, so they keep their names). Originals come out first, then analogs tagged in the second column — `CHEMBLID|from:<seed>|sim:<tanimoto>` — which mol2mol ignores but lets you audit where each molecule came from.

**The two knobs that matter.** `--similarity` is your relevance/volume tradeoff: 70 pulls tight analogs; drop to 60 if check_pairs still looks thin and you need bigger, more overlapping clusters (at the cost of looser neighbors). `--max-per-seed` caps cluster size so one popular scaffold doesn't swamp the set — it's the data-prep mirror of `pairs.max_cardinality`. At 17 seeds × 25 that's up to ~425 analogs before dedup, comfortably into the pair count you want.

**A few honest caveats.** if `molecules` or `molecule_structures.canonical_smiles` ever come back shaped differently, that's the spot to check. Second, ChEMBL hits are _structural_ neighbors, not compounds known to be active against your pocket — which is exactly right for mol2mol TL (you're teaching edit style, not activity), but don't mistake the analogs for validated actives. Third, some seeds may return zero analogs if they're novel or proprietary chemotypes ChEMBL doesn't cover; the script lists those at the end, and they'll stay orphans until you add in-house analogs or R-group enumeration for them.

Once check_pairs shows a few hundred pairs at ≥0.5–0.6 with most molecules covered, set `pairs.lower_threshold` to that value in your TOML, point `smiles_file`/`validation_smiles_file` at the augmented split, and you're ready to actually run TL.
### How to run it:

```bash
python scripts/tl/03_augment_from_chembl.py configs/compounds.smi -o configs/compounds_augmented.smi --similarity 70 --max-per-seed 25

python scripts/tl/02_check_pairs.py configs/compounds_augmented.smi   # confirm the table fills in
```

### Example Output (Unbalanced)

```bash
python scripts/tl/03_augment_from_chembl.py configs/compounds.smi -o configs/compounds_augmented.smi --similarity 60 --max-per-seed 25 
seed1: +1 analogs 
seed2: +2 analogs 
seed3: +4 analogs 
seed4: +0 analogs 
seed5: +13 analogs 
seed6: +7 analogs 
seed7: +4 analogs 
seed8: +1 analogs 
seed9: +2 analogs 
seed10: +20 analogs 
seed11: +2 analogs 
seed12: +5 analogs 
seed13: +22 analogs 
seed14: +6 analogs 
seed15: +3 analogs 
seed16: +5 analogs 
seed17: +2 analogs 
seeds=17 analogs_added=99 total_unique=116 -> configs/compounds_augmented.smi seeds with 0 analogs (1): seed4 N
ow re-run check_pairs.py on the output to see the threshold table fill 

scripts/tl/02_check_pairs.py configs/compounds_augmented.smi 
116 valid molecules -> 13340 ordered pairs possible (excl. self) threshold | pairs >= thr | sources w/>=1 target 
------------------------------------------------ 
0.3 | 1568 | 115 
0.4 | 1348 | 115 
0.5 | 1058 | 115 
0.6 | 832  | 115 
0.7 | 374  | 78 
0.8 | 124  | 28
```

Compared to the 0.70 run it roughly tripled everything: 116 molecules, and at 0.6 it went from 286 pairs / 47 covered to **832 pairs / 115 covered**. Only one true orphan left (seed4), and full coverage holds all the way down to 0.5. 

- **0.5** → 1058 pairs, all 116 covered, but looser edits (you're training on pairs only ~50% similar).
- **0.6** → 832 pairs, ~full coverage, moderate edit size. 
- **0.7** → 374 pairs, tightest/most conservative edits, but ~38 molecules form no pairs and drop out (78 covered). Use only if you want strictly close-in analogs and don't mind the smaller effective set.

So I'd set `lower_threshold = 0.6` on this set and move on. 

**The cluster imbalance is now more pronounced** — seed13 (+22), seed10 (+20), and seed5 (+13) account for 55 of your 99 analogs. Those three chemotypes will dominate the learned bias. For a first run that's acceptable, and it may just reflect that those scaffolds are well-studied in ChEMBL. If you want a more even bias across your 17, re-run with `--max-per-seed 10` or `12` — that trims the big three without much affecting the thin clusters, and you'll still clear several hundred pairs. Worth a quick comparison, but not mandatory.

seed4 is a genuine orphan — zero analogs at both 60 and 70, so ChEMBL simply has nothing near it. That's a novel or proprietary scaffold; it'll sit out of TL unless you generate analogs for it in-house. Nothing the augmentation script can do about that one.

### Example Output (Well Balanced after rerunning with --max-per-seed 10)

```bash
python scripts/tl/03_augment_from_chembl.py configs/compounds.smi -o configs/compounds_augmented.smi --similarity 60 --max-per-seed 10 
seed1: +1 analogs 
seed2: +2 analogs 
seed3: +4 analogs 
seed4: +0 analogs 
seed5: +9 analogs 
seed6: +7 analogs 
seed7: +4 analogs 
seed8: +1 analogs 
seed9: +2 analogs 
seed10: +8 analogs 
seed11: +2 analogs 
seed12: +5 analogs 
seed13: +8 analogs 
seed14: +6 analogs 
seed15: +3 analogs 
seed16: +5 analogs 
seed17: +2 analogs 
seeds=17 analogs_added=69 total_unique=86 -> configs/compounds_augmented.smi seeds with 0 analogs (1): seed4 Now re-run check_pairs.py on the output to see the threshold table fill in. 

python scripts/tl/02_check_pairs.py configs/compounds_augmented.smi 
86 valid molecules -> 7310 ordered pairs possible (excl. self) 
threshold | pairs >= thr | sources w/>=1 target 
------------------------------------------------ 
0.3 | 694 | 85 
0.4 | 488 | 85 
0.5 | 376 | 85 
0.6 | 290 | 85 
0.7 | 164 | 52 
0.8 | 66 | 13
```

The balancing worked exactly as intended — the analog counts now spread evenly from 1 to 9 instead of three runaway clusters of 13–22. No chemotype dominates anymore, coverage is still 85 of 86 (only seed4, your permanent orphan), and you've got 290–376 pairs across the 0.5–0.6 band. This is a clean, well-balanced training set.

Here's the call between your two good options, because you've now got both on disk:

This balanced set (86 mol, cap 10) trades raw pair count for even representation. The 116-molecule set had more pairs (832 at 0.6) but half of them came from three chemotypes. Since the whole point is to bias toward your _series of 17_ — not toward whichever three scaffolds happen to be well-studied in ChEMBL — balance is worth more than volume here. 290–376 pairs with full coverage is plenty for a focused TL run, especially as a warm start. 

For the threshold, use `lower_threshold = 0.5` on this set: 376 pairs, all 85 covered, and since this set is a bit leaner on pairs you want to keep the count up. Use 0.6 (290 pairs, still full coverage) if you specifically want tighter edits. Skip 0.7 here — it drops you to 52 covered, undoing the balancing.

## Step 4 — Split (`scripts/tl/04_split_smiles.py`)

Run: `python scripts/tl/04_split_smiles.py configs/compounds.smi --val-frac 0.15`

### Example Output (Trying different seeds: none, 7, 15, 25)

```bash
# No seed paramater
python scripts/tl/04_split_smiles.py configs/compounds_augmented.smi --val-frac 0.15
total=86 train=73 val=13 
#check resultung pairs
python scripts/tl/02_check_pairs.py configs/compounds_augmented_val.smi
13 valid molecules -> 156 ordered pairs possible (excl. self) 
threshold | pairs >= thr | sources w/>=1 target 
------------------------------------------------ 
0.3 | 10 | 7 
0.4 | 10 | 7 
0.5 | 10 | 7 
0.6 | 6  | 6 
0.7 | 4  | 4 
0.8 | 2  | 2

# Trying --seed 7
python scripts/tl/04_split_smiles.py configs/compounds_augmented.smi --val-frac 0.15 --seed 7 
total=86 train=73 val=13 
#check resultung pairs
python scripts/tl/02_check_pairs.py configs/compounds_augmented_val.smi 13 valid molecules -> 156 ordered pairs possible (excl. self) 
threshold | pairs >= thr | sources w/>=1 target 
------------------------------------------------ 
0.3 | 12 | 6 
0.4 | 10 | 6 
0.5 | 6  | 3 
0.6 | 6  | 3 
0.7 | 6  | 3 
0.8 | 6  | 3


## Try --seed 15
python scripts/tl/04_split_smiles.py configs/compounds_augmented.smi --val-frac 0.15 --seed 15 
total=86 train=73 val=13
#check resultung pairs
python scripts/tl/02_check_pairs.py configs/compounds_augmented_val.smi 13 valid molecules -> 156 ordered pairs possible (excl. self) 
threshold | pairs >= thr | sources w/>=1 target 
------------------------------------------------ 
0.3 | 26 | 11 
0.4 | 20 | 11 
0.5 | 8  | 8 
0.6 | 4  | 4 
0.7 | 2  | 2 
0.8 | 2  | 2 

## Try --seed 25
python scripts/tl/04_split_smiles.py configs/compounds_augmented.smi --val-frac 0.15 --seed 25 
total=86 train=73 val=13
#check resultung pairs
python scripts/tl/02_check_pairs.py configs/compounds_augmented_val.smi 13 valid molecules -> 156 ordered pairs possible (excl. self) 
threshold | pairs >= thr | sources w/>=1 target 
------------------------------------------------ 
0.3 | 10 | 7 
0.4 | 8  | 7 
0.5 | 4  | 4 
0.6 | 4  | 4 
0.7 | 2  | 2 
0.8 | 0  | 0
```

### Interpreting Results
Across the default split (10 pairs) and seeds 7, 15, 25 (6, 8, 4 pairs at 0.5), you're stuck in the same thin band. That's not bad luck you can seed your way out of; it's the structural floor of holding 13 molecules out of an 86-molecule set. No seed is going to hand you 20 robust val pairs.

- **Default split (no `--seed`)**: 10 pairs / 7 covered — the most pairs.
- **Seed 15**: 8 pairs / 8 covered — the best coverage.

For a validation set, the thing that actually matters is having enough pairs to estimate the loss, so I'd lean to the default (10 pairs) — coverage is a _training_ concern, not a val one. But they're close enough that it doesn't matter. 

Then judge overfitting primarily by sampling from each checkpoint — valid-SMILES rate staying high, outputs resembling your series without collapsing into memorized copies — and use the val curve only as a coarse "is the gap widening" sanity check. That's all a 10-pair val set can honestly tell you, and it's enough.

One alternative if the weak val genuinely bothers you: with data this scarce, you can also just train on all 86 molecules and skip the held-out val entirely, monitoring purely by sampling at each checkpoint. That puts every pair into training, which is where the value is when you're pair-limited. It's a legitimate small-data approach — you trade the (already weak) overfitting curve for ~15% more training signal.

## Step 5 — TL config (`configs/tl_mol2mol.toml`)


```toml
run_type = "transfer_learning"
device = "cuda:0" #change to CPU for GPU-less machines
tb_logdir = "tb_TL"

[parameters]
input_model_file = "/workspace/configs/priors/mol2mol_medium_similarity.prior"  # point at the prior in your image
smiles_file = "/workspace/configs/compounds_train.smi"
validation_smiles_file = "/workspace/configs/compounds_val.smi"
output_model_file = "/workspace/models/tl_mol2mol.model"

num_epochs = 25            # small set overfits fast; watch the curve, take an early checkpoint
save_every_n_epochs = 2    # frequent, so you can pick the best epoch rather than the last
batch_size = 64            # lower toward 16-32 if check_pairs showed only a couple hundred pairs

pairs.type = "tanimoto"
pairs.lower_threshold = 0.5   # SET THIS from check_pairs.py output
pairs.upper_threshold = 1.0
pairs.min_cardinality = 1     # keep every source with >=1 valid target
pairs.max_cardinality = 199   # irrelevant at this size; harmless to leave
```

Two notes. `input_model_file`: the tutorial's `priors/mol2mol_medium_similarity.prior` is the right _default_ (medium-similarity biases toward analogs without collapsing to near-duplicates); if you want tighter, closer-in analogs of your 16, the high-similarity prior is the alternative. Confirm where the priors actually live inside `htvs-pipeline:latest` and fix the path accordingly. And `num_epochs = 25` is a starting point — the real stop signal is the loss curves and the sampling check, not the number.

**To run it**, simplest is interactively in the container:

```bash
reinvent -l tl.log configs/tl_mol2mol.toml

#using docker compose instead
docker compose up tl
```

### Example results

At the final step (25), your three curves ended at roughly: 
* training loss: 3.89
* validation loss: 4.71
* sample loss: 3.47

The validation loss landing a bit above training (a gap of ~0.8) is normal and expected — a model always does somewhat better on the molecules it studied than on held-out ones, so a small persistent gap is the healthy default, not a problem. What you'd worry about is the validation line _climbing steeply away_ from training, and that's not what happened here. In fact the "smoothed" value sitting above the actual "value" for all three curves tells you none of them had turned upward at the end — they were all still gently drifting down. So this is mild, well-behaved learning with only slight overfitting, which is about the best you can hope for from a small set. (The "sample loss" being lowest is just the model being confident about molecules it generates itself; it's not a generalization metric)

The final epoch-25 model is fine to use. If, when you sample from it, the outputs look like near-copies of your training molecules, that's the sign to fall back to the checkpoint around epoch 16–18 (near where validation flattened) — but try the final model first.

## Step 6 — sampling config (`configs/tl_mol2mol.toml`)


```toml
run_type = "sampling"
device = "cpu"
tb_logdir = "tb_sampling"

[parameters]
model_file = "/workspace/models/tl_mol2mol.model"      # your TL output (or an epoch-16/18 .chkpt)
smiles_file = "/workspace/configs/compounds.smi"        # SEEDS: the 17 molecules you want analogs OF
output_file = "/workspace/results/sampled.csv"
sample_strategy = "multinomial"              # diverse pool; "beamsearch" = deterministic top analogs
temperature = 1.0                            # multinomial only; ↑ = more diverse, ↓ = more conservative
num_smiles = 1000                            # total molecules to generate (tune up for a bigger pool)
unique_molecules = true                      # dedupe + canonicalize
randomize_smiles = true
```

Three choices worth understanding:

_What to seed with._ Use your original 17 (`configs/compounds.smi`), not the augmented 86. The augmentation was scaffolding to train the model; now you want novel analogs of the compounds you actually care about for your pocket. The seeds define which neighborhoods get explored, so seed with your real targets.

_Sampling strategy._ Mol2Mol supports either multinomial sampling with temperature or beam search. Multinomial is stochastic — it gives you a varied pool, which is what you want for feeding a downstream screen. Beamsearch is deterministic and returns the few most-probable analogs of each seed: a tighter, reproducible, smaller set. For hit-finding I'd start with multinomial; switch to beamsearch if you specifically want "the best handful of analogs per compound."

_How many._ `num_smiles` is roughly the total number generated across your seeds, so ~1000 over 17 seeds is ~60 attempts per compound before dedup. Bump it to several thousand once you've confirmed the output looks good.

**To run it**, simplest is interactively in the container:

```bash
reinvent -l sampling.log configs/sample_mol2mol.toml

#using docker compose instead
docker compose up sample
```


**Then triage `results/sampled.csv`.** This is the validity check I mentioned earlier, now with real output. You want to confirm: the SMILES are valid and varied (not degenerate repeats), they're genuinely _novel_ (not just your 17 seeds echoed back), and they sit in a sensible similarity band to their seeds (close enough to keep the properties, different enough to be new). That filtered, novel pool is what eventually goes into docking against your pocket.

## Step 7 — scripts/tl/05_qc_sampled.py
 

```bash
python scripts/tl/05_qc_sampled.py results/sampled.csv \
    --train configs/compounds_augmented_train.smi \
    -o results/sampled_qc.csv
```

It auto-detects the SMILES and input columns. If REINVENT's actual column names differ from what it expects, it'll print the real column names and you just pass `--smiles-col`/`--input-col`.

**How to read the summary** (plain-language, since the numbers are the point):

_Valid %_ — should be high, ideally >90%. The sampler canonicalizes its output so this is usually fine; a low number means the model drifted off into nonsense and you'd want an earlier checkpoint.

_Unique %_ — how much variety it produced. Very low means it's repeating itself and you should sample more or raise the temperature.

_Memorized (in train)_ and _seed returned unchanged_ — these are your two "wasted output" numbers, and you want both **low**. Memorized means it spat back a molecule that was literally in its training set; seed-unchanged means it handed your input back without editing it. Both are useless to you because you already had those compounds. If "memorized" is high, that's the concrete sign of the overfitting we talked about — fall back to the epoch-16/18 checkpoint and re-sample.

_Novel & unique_ — **this is the number that matters.** It's the count of genuinely new molecules you didn't already have. This is your usable pool.

_Tanimoto-to-seed distribution_ — the sweet spot is mass in the **moderate-to-close bands (0.4–0.8)**: close enough to your seeds to likely keep their properties, different enough to be worth making. A big pile in "identical (~1.0)" means it's barely changing anything; a big pile in "distant <0.4" means it's wandering off into unrelated chemistry and may have lost the thread of your series. A healthy run is a hump in the middle.

So a good result looks like: high validity, low memorized + low seed-copy, a solid novel-&-unique count, and Tanimoto sitting mostly in 0.4–0.8.

The annotated CSV (`results/sampled_qc.csv`) is per-molecule with all the flags, so your docking input is just the rows where `novel=1` and `duplicate=0` — that filtered set is exactly what feeds the `dock` stage we set up. If you want, once you've run real sampling and have actual numbers, paste the summary back and I'll tell you whether it looks healthy or whether to adjust temperature/checkpoint before you spend time docking.


### Example Output

```bash
python ./scripts/tl/05_qc_sampled.py results/sampled.csv --train configs/compounds_augmented_train.smi -o results/sampled_qc.csv
=== QC: results/sampled.csv ===
rows (generated):        2499
valid:                   2499  (100.0%)
unique (valid):          2499  (100.0% of valid)
memorized (in train):    67  (2.7% of valid)
seed returned unchanged: 17  (0.7% of valid)
novel & unique:          2432  <- your usable pool

Tanimoto-to-seed (ECFP4):
  min 0.12 | 25% 0.46 | median 0.55 | 75% 0.63 | max 1.00 | mean 0.55
  distribution:
    identical (~1.0) :    18  (0.7%)
    very close .8-1  :    50  (2.0%)
    close .6-.8      :   811  (32.5%)
    moderate .4-.6   :  1325  (53.0%)
    distant <.4      :   295  (11.8%)
```


Validity 100% and uniqueness 100% are about as good as it gets — the model stayed completely on the rails and didn't repeat itself once across 2499 molecules. Memorized at 2.7% and seed-unchanged at 0.7% are both low, which is the concrete confirmation that the mild overfitting we saw in the loss curve didn't actually hurt — the final model isn't regurgitating its training set, so there's no need to drop back to an earlier checkpoint. And the Tanimoto distribution is exactly the shape you want: 85.5% of the output sits in the 0.4–0.8 sweet spot (53% moderate + 32.5% close), with a median of 0.55. That's the central hump — close enough to your seeds to likely carry their properties, different enough to be worth making.

The only genuine judgment call in here is the **11.8% distant tail** (295 molecules under 0.4 similarity to their seed). Those are your model's most exploratory outputs — either interesting scaffold-hops or molecules that drifted far enough to lose what made the series matter. Nothing's wrong with 11.8%; it's a minority. You just get to decide: keep them if you want shots at novel scaffolds against your pocket, or drop them if you'd rather stay faithful to the known chemistry. For a first docking pass I'd lean toward keeping them — docking is the cheap arbiter of whether they're any good.

## What next?

sampled.csv came from plain sampling, which has none of that, so to dock you pick one of two paths:

**Path A — dock the static library (a one-shot virtual screen).** Take sampled.csv and screen it _outside_ REINVENT: prep your receptor once (PDB → add hydrogens/charges → define a box around the pocket → `receptor.pdbqt`), embed each generated SMILES to 3D and convert to a docking-ready format (RDKit/Corina/LigPrep), run a docking engine (Vina/smina/gnina, or DockStream in batch "console" mode), collect the best score per molecule, then rank and filter. A common hit cutoff is a docking score ≤ −8 kcal/mol with QED ≥ 0.7. Docking happens once, after sampling. This is the simple route.

**Path B — closed-loop optimization (the canonical REINVENT pipeline).** Here you _don't_ dock sampled.csv at all. You run `staged_learning` with your TL model as the `agent_file`, and a scoring function that contains a DockStream docking component. The model generates → DockStream docks against your pocket → the score steers the model → repeat for many steps, so it _learns_ to produce molecules that dock well. The docking component is configured inside the stage's scoring block, pointing at a DockStream JSON config and docker script. This is exactly the target-directed strategy you were thinking of, and it's where TL was always headed — TL is the warm-start for this RL run.
