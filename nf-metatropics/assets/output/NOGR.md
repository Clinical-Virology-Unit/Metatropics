# NOGR (Non-Overlapping Genomic Regions)

Virasign reports **NOGR** as an extra per-virus evidence metric: it counts how many **distinct, non-overlapping regions of the reference genome** are supported by at least one mapped read, and the total number of reference bases spanned by those regions. This follows an mNGS reporting idea used before that can be very valuable for interpretation: requiring a minimum number of **non-overlapping viral reads/regions** to support a detection (e.g. Nature Communications [`s41467-024-51470-y`](https://www.nature.com/articles/s41467-024-51470-y), which uses a threshold of ≥3 non-overlapping viral reads/contigs as a positive criterion).

----

## Why this is useful

Coverage breadth can be low for true positives (degraded samples, low viral load), but it can also be low for technical artifacts such as **amplicon contamination**, where many reads stack on the same small viral genomic region. **NOGR helps distinguish real signals from false contamination**.

----

## What exactly is computed

Virasign looks at the virus BAM file (`{accession}.bam`) and does this:

- **Step 1**: For each mapped read, record where it maps on the virus genome (**start** and **end** coordinates).
- **Step 2**: Count how many reads map to **different parts** of the genome (reads that do not overlap each other on the genome).

This gives two values:

- **`nogr_regions`**: how many **separate (non-overlapping) genomic regions** are supported
- **`nogr_bases`**: how many **virus genome bases** are covered by those separate regions (so it can never be larger than the genome length)

----

