# Range-based MAF construction

```bash
panmanUtils -A all.maf -N tree.nwk -o all --backend auto
panmanUtils -I panman/all.panman --list-chromosomes
panmanUtils -I panman/all.panman --fasta --chromosome chrI -o chrI
panmanUtils -I panman/all.panman --maf -o exported
```

`auto` is the default. CPU-only builds, Macs, and machines without an available
CUDA device use CPU inference. `--backend cpu` forces CPU. `--backend cuda`
requires a CUDA-enabled build and an available device and fails clearly if
unavailable. In auto mode, a CUDA inference error reports the reason and retries
the batch on CPU. File reading, layout, refinement, and serialization use the
CPU. Independent inference columns run through TBB on CPU or a CUDA kernel.

## Progress logging

MAF construction prints phase boundaries and rate-limited progress to stderr.
Reading reports input bytes and percentage, completed blocks, rows, and genomes.
Later phases report chromosome paths, pseudo-root instances, refinement merges
and work used, inference blocks/columns, and elapsed seconds. Intermediate
updates appear at approximately five-second intervals when the current phase
reaches a reporting point. Completion messages are always flushed immediately.

```bash
panmanUtils -A all.maf -N tree.nwk -o all --backend auto 2>&1 | tee construction.log
```

A process already running uses its original executable; rerun the command after
rebuilding to see the new logging.

## Construction and refinement

Each MAF alignment record defines a homology group, independently of identical
DNA sequence or chromosome name. Source identifiers match a Newick tip exactly
(single chromosome `chr0`), or use `<tip>.<chromosome-or-scaffold>`. The longest
matching tip prefix is used, so tip names can contain dots. Differently named
chromosomes/scaffolds in one alignment record are supported without filtering
out genomes. The Newick input may span multiple lines.

The importer sorts each genome's chromosome occurrences by forward genomic
coordinate and progressively inserts those paths into a sequential pseudo-root.
Matching block identities are chained using sparse LCS. There is no dense
pseudo-root-length by path-length dynamic-programming matrix. For a bounded
matching-pair count the chain is exact; repetitive identities exceeding the
budget are omitted as anchors and represented by additional copies. No input
occurrence is dropped by that fallback.

Refinement builds ordering edges only from consecutive occurrences on each tip
chromosome. Homologous copies can merge only if no genome would use the merged
instance twice and neither node can reach the other in the current graph.
Reachability and candidate work have a budget; uncertainty rejects the merge.
Commits are sequential and checked against the current graph. A deterministic
topological ordering chooses between a stable ordering and a continuity-biased
ordering, accepting the latter only if it reduces tip chromosome range count.
This improves a progressive heuristic; it does not promise a globally minimal
pseudo-root. Block IDs and ranges are finalized before mutation inference.

Defaults and controls:

```bash
panmanUtils -A all.maf -N tree.nwk -o all \
  --backend auto --match-budget 1000000 --refine-passes 2 --refine-budget 2000000
```

`--refine-passes 0` disables copy merging. `--refine-budget 0` also disables
merging. Work budgets preserve sequence correctness but may increase duplication.
Input rows and final blocks still consume memory proportional to the data;
removing quadratic alignment does not make total memory constant. The legacy
sequence reconstruction representation can also be large for whole-genome
FASTA/MAF export. Full UCSC datasets should be measured on their intended host;
small fixtures are not a whole-genome memory benchmark.

## Representation and old files

The Cap'n Proto file ID and all existing ordinals and field types are unchanged.
Only the following fields are appended to `Tree`:

```capnp
sequenceChromosomeRanges @11: List(SequenceChromosomeRanges);
homologyGroups @12: List(HomologyGroup);
```

Ranges are half-open intervals of sequential **primary** block IDs. A range can
span blocks absent in that genome, but cannot overlap another chromosome range
in that genome. The existing block mutations encode presence and orientation;
new ranges do not duplicate nucleotide sequences. Copies have distinct block
IDs and explicit homology-group membership. Consensus storage may independently
deduplicate identical sequences; that does not assert biological homology.

Legacy `sequenceBlockOrders` remains populated with each genome's concatenated
chromosome paths. Legacy `chrLists` stores a disjoint pseudo-root catalog,
assigning a block to its first encountered owner; it is not authoritative for
per-tip chromosome assignment. Its names also supply the file-local numeric
`--chrID` catalog. Use `--chromosome NAME` for stable selection across MAF
export/re-import, which may change first-appearance order and numeric IDs.

Old files without the new fields retain the old reconstruction path. Tests
cover pre-range files, historical Text Newick encoding, and legacy protobuf
conversion. New files remain structurally readable by the old schema and keep
legacy whole-genome order. Unmodified old tools cannot interpret per-tip
chromosome assignments and may drop the new metadata on rewriting. Do not use
an old writer to edit a new range-based file.

FASTA, MAF, summaries, Newick, rerooting, annotation, and subnetwork extraction
preserve or use the new metadata. Global `--fasta-aligned` output remains an MSA
in pseudo-root order. VCF, GFA, Usher, tree-sequence export, chromosome-sensitive
indexing/translation, block deduplication, and network creation from these new
files currently fail explicitly rather than applying incompatible assumptions.
Legacy files retain their existing operations.

Input MAFs need not cover whole chromosomes. Reconstruction includes only the
supplied row occurrences, including repetitions/overlaps if present, and does
not invent uncovered bases. Exported coordinates describe those reconstructed
segments; original uncovered coordinate intervals are not preserved. A genome
with no input rows, or no rows on a requested chromosome, has an empty sequence.

## Build and validation

CUDA compilation is optional and detected by CMake. Use
`-DPANMAN_ENABLE_CUDA=OFF` to build CPU-only. On CUDA hosts, configure an
appropriate `CMAKE_CUDA_ARCHITECTURES` for the target GPUs. When no architecture
is set, CMake 3.23+ builds for its `all-major` architecture set; earlier CMake
versions require an explicit architecture for CUDA. The previous
standalone files in `gpu/` are retained; the integrated backend is
`gpu/constructionInference.cu`, using the same unit-cost Sankoff implementation
as CPU inference. Both support polytomies, missing block observations, DNA and
protein states, and deterministic lowest-state tie breaking. Batches bound the
inference workspace instead of constructing a genome-by-block matrix.

```bash
cmake --build build --target panman_layout_tests
ctest --test-dir build --output-on-failure
# Require actual GPU execution rather than allowing an unavailable-device skip:
python3 tests/test_construction.py ./build/panmanUtils \
  --capnp ./build/source-env/bin/capnp --require-cuda
```

Tests include randomized sparse-vs-dense LCS checks, 20,000-block sparse input,
repeat budgets, cycle rejection, copy-preserving refinement, exact DNA/protein
reconstruction, cross-chromosome/scaffold MAF round-trips, metadata corruption,
reroot/subnetwork preservation, and golden old formats. CPU/auto outputs are
compared including states and metadata. CUDA parity runs when a device exists;
`--require-cuda` turns an unavailable backend/device into a test failure. CI
builds CPU configurations and a CUDA development container; CUDA hardware tests
require an NVIDIA GPU host and cannot be established by a compiler-only job.

## FASTA output defaults

For files with per-tip chromosome ranges, `--fasta` writes one nonempty record
per genome/chromosome, using unique identifiers such as `>ce11.chrI` and
`>ce11.chrM`. Missing chromosomes are omitted. `--fasta --concatenate` retains
one whole-genome record per tip. An explicit `--chromosome NAME` or `--chrID N`
still selects one chromosome and retains the existing selected-output header.
Files without per-tip range metadata retain their previous default behavior.

```bash
panmanUtils -I panman/small.panman --fasta
panmanUtils -I panman/small.panman --fasta --concatenate
```

Progressive LCS alignment uses CPU TBB workers, including with `--backend auto`
and `--backend cuda`. `--threads` controls the shared TBB worker limit. One path
seeds the pseudo-root; subsequent batches align against an immutable root with a
shared homology index. Commits check sample ownership and rebuild the root once
per batch. Batches contain at most 64 paths and use a conservative 256 MiB
scratch estimate to limit simultaneous chaining jobs; one oversized path runs
alone. This estimate is not a process memory limit: the root, shared index, input,
and refinement graph require additional memory.

Batch boundaries and commit order are independent of worker scheduling and thread
count. Batching can create more initial copies than the previous sequential
heuristic because jobs cannot match insertions made by another job in the same
batch. Bounded refinement can merge compatible copies; it does not guarantee
removing every extra copy. Tip sequences, chromosome ownership, and ordering are
preserved. GPU inference and the on-disk schema are unchanged.

Adaptive layout batches estimate scratch from incoming path lengths and indexed
matching-pair counts, capped by the LCS match budget, rather than reserving the
full budget for every short path. Each batch has at most 64 paths and a 256 MiB
scratch estimate; one oversized job runs alone. Batches that insert no instances
reuse both the root order and its position index. Sample ownership is always
checked again, including when the index is reused. Progress and final layout
statistics report index builds, reuses, batch count and largest batch. Changed
batch boundaries can change block duplication relative to earlier builds.
