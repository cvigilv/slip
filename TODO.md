# To-Do List

- [x] _Clean input_: Remove entries with similarity\_measure == -99 and with TP == 1.

- [x] _MySQL comparison_: Compare each ligand-protein pair with ChEMBL\_#\_db for temporal validation of predictions.

- [x] _Similarity Measure cutoff_: Filter out all entries with a similarity\_measure outside the given range. Range is inclusive for upper limit and non-inclusive for lower limit.

- [x] _Pfam comparison_: Count amount of common Pfam between query\_ligand and hit\_target, creating new column. Use this column to create a boolean for further filtering.

- [x] _Max Phase cutoff_: Filter out all entries with max\_clinical\_phase outside given range. Range is inclusive for upper limit and non-inclusive for lower limit.

- [x] _TopX selection_: Keep top _n_ selections of filtered output file. Use pandas _df.head(N)_ for this.

- [ ] _Add Molport ID for in-stock compounds_: Add identificator for easier search in Molport store.
