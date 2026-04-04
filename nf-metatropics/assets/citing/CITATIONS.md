# Citations

To cite **Metatropics** itself, use the metadata in [`CITATION.cff`](../../../CITATION.cff) at the repository root (GitHub shows this under **Cite this repository**). The rest of this file lists citations for **nf-core**, **Nextflow**, **container engines**, and **third-party software invoked by this pipeline** (see `nf-metatropics/modules/` and `workflows/`).

---

## nf-core / Nextflow / containers

### [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x). PubMed PMID: 32055031.

### [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820). PubMed PMID: 28398311.

### [Docker](https://linuxjournal.com/content/docker-lightweight-linux-containers-consistent-development-and-deployment)

> Merkel D. Docker: lightweight Linux containers for consistent development and deployment. *Linux Journal*. 2014;2014(239):Article 2, p. 2. Published 2014 Mar 1. Alternate listing: [ACM Digital Library](https://dl.acm.org/doi/10.5555/2600239.2600241).

### [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

> Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: [10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459). PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

---

## Pipeline software

Alphabetical. Match names to modules under `nf-metatropics/modules/`.

### bam-readcount

Used for per-position base counts in BAMs. Project page: [github.com/genome/bam-readcount](https://github.com/genome/bam-readcount). If you need a journal reference for the stack, cite the host toolkit context you rely on (e.g. your variant-calling workflow); the tool itself is distributed as open source from the McDonnell Genome Institute.

### BBMap / BBTools (`reformat.sh`, rarefaction)

BBTools suite (JGI). Primary distribution and documentation: [SourceForge: BBMap](https://sourceforge.net/projects/bbmap/). Typical attribution: Bushnell B (2014) *BBMap* / BBTools, Lawrence Berkeley National Lab/JGI. For a related peer-reviewed component of the suite, see e.g. Bushnell B. BBMerge: Accurate paired shotgun read merging via overlap. *PLOS ONE* 9(4):e93362. doi: [10.1371/journal.pone.0093362](https://doi.org/10.1371/journal.pone.0093362).

### Dorado (ONT basecalling / demultiplexing)

Oxford Nanopore Technologies. Source code and release notes: [github.com/nanoporetech/dorado](https://github.com/nanoporetech/dorado). Cite the **Dorado version** you ran and, where appropriate, the Oxford Nanopore Technologies basecalling documentation for your chemistry/flow-cell.

### fastp

> Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*. 2018 Sep 1;34(17):i884-i890. doi: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560). PubMed PMID: [30423086](https://pubmed.ncbi.nlm.nih.gov/30423086/); PubMed Central PMCID: [PMC6129281](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/).

### Homopolish

> Huang YT, Liu PY, Shih PW. Homopolish: a method for the removal of systematic errors in nanopore sequencing by homologous polishing. *Genome Biol*. 2021;22:95. doi: [10.1186/s13059-021-02282-6](https://doi.org/10.1186/s13059-021-02282-6). Code: [github.com/ythuang0522/homopolish](https://github.com/ythuang0522/homopolish).

### iVar

> Grubaugh ND, Gangavarapu K, Quick J, Matteson NL, De Jesus JG, Main BJ, Tan AL, Paul LM, Brackney DE, Grewal S, Gurfield N, Van Rompay KKA, Isern S, Michael SF, Coffey LL, Loman NJ, Andersen KG. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. *Genome Biol*. 2019 Jan 8;20(1):8. doi: [10.1186/s13059-018-1618-7](https://doi.org/10.1186/s13059-018-1618-7). PubMed PMID: [30621750](https://pubmed.ncbi.nlm.nih.gov/30621750/); PubMed Central PMCID: [PMC6325816](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6325816/). Software: [github.com/andersen-lab/ivar](https://github.com/andersen-lab/ivar).

### Medaka

Oxford Nanopore Technologies neural-network consensus/polishing. Repository: [github.com/nanoporetech/medaka](https://github.com/nanoporetech/medaka). The Medaka README asks users to cite the **repository** and the **Oxford Nanopore Technologies** product documentation; check the README of the Medaka version bundled in your run for the exact wording.

### MetaMaps

> Dilthey AT, Jain C, Koren S, Phillippy AM. Strain-level metagenomic assignment and compositional estimation for long reads with MetaMaps. *Nat Commun*. 2019 Jul 8;10:2926. doi: [10.1038/s41467-019-10934-2](https://doi.org/10.1038/s41467-019-10934-2). PubMed PMID: [31296857](https://pubmed.ncbi.nlm.nih.gov/31296857/).

### minimap2

> Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*. 2018 Sep 15;34(18):3094-3100. doi: [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191). PubMed PMID: [29750262](https://pubmed.ncbi.nlm.nih.gov/29750262/). GitHub: [github.com/lh3/minimap2](https://github.com/lh3/minimap2).

### NanoPlot (NanoPack family)

> De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics*. 2018 Aug 1;34(15):2666-2669. doi: [10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149). PubMed PMID: [29547981](https://pubmed.ncbi.nlm.nih.gov/29547981/). Repository: [github.com/wdecoster/nanoplot](https://github.com/wdecoster/NanoPlot).

### R (figures / summaries: e.g. MetaMaps plots, read counts, coverage PDFs)

> R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. [https://www.R-project.org/](https://www.R-project.org/) (year of the R version you used).

For **ggplot2**-based graphics (common in tidyverse stacks), many workflows cite:

> Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. ISBN 978-3-319-24277-4. [https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org).

### samtools

> Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *Gigascience*. 2021 Feb 16;10(2):giab008. doi: [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008). PubMed PMID: [33590861](https://pubmed.ncbi.nlm.nih.gov/33590861/). Project: [www.htslib.org](https://www.htslib.org/). For the original SAM/BAM format and early SAMtools: Li H, et al. *Bioinformatics*. 2009;25(16):2078-2079. doi: [10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352).

### seqtk

> Li H. seqtk: Toolkit for processing sequences in FASTA/FASTQ format. [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk) (cite the commit or release tag used in your container/conda environment).

---

## Notes

- **Versions**: For reproducibility, cite the **software versions** recorded in your run under `Summary/pipeline_info/` (`software_versions.yml`), not only the papers above.
- **Perl / shell / Python** one-liners (e.g. header cleanup) are standard utilities; no separate citation is usually required beyond this file and your methods text.
- If you add or swap tools in a fork, update this list to match `modules/` and `workflows/`.
