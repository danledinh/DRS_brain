# Direct RNA Sequencing reveals epitranscriptomic regulation of brain cells and Alzheimer’s Disease pathology

## Abstract
Alternative mRNA splicing and post-transcriptional RNA modification are key mechanisms that regulate transcript function; however, their role in neuronal activity and neurodegenerative disease remains poorly defined. In this study, we evaluated two nanopore-based long-read sequencing (LR-seq) formats: cDNA-PCR sequencing (CPS) and direct RNA sequencing (DRS). We then applied DRS to profile both full-length isoforms and RNA modifications in major brain cell types derived from induced pluripotent stem cells (iPSCs) and post-mortem Alzheimer's disease (AD) brains. Relative to CPS, DRS achieved higher accuracy and sensitivity for transcript quantification, de novo transcript model construction, and open reading frame (ORF) annotation across neuropathological gene sets. Focusing on iPSC-derived neurons, we built a multi-omic atlas to connect transcriptional output with translational engagement and protein abundance, by integrating DRS-based mRNA abundance, N6-methyladenosine (m6A) status and poly(A) tail length with ribosome profiling (Ribo-seq) and mass spectrometry (MS). The combination of DRS and Ribo-seq data demonstrated synergism in predicting protein abundance. This analysis also uncovered a significant inverse relationship between m6A modification and mRNA abundance, which was dependent on the engagement of the ribosomal A-site. Lastly, we applied DRS to the epitranscriptomic analysis of AD brain samples, demonstrating that m6A profiles can be used to distinguish early- versus late-stage disease. 

## Software versions and references
Pandas v1.5.3
IsoQuant v3.7.0
pycoQC v2.5.2
NumPy v1.26.4
SciPy v1.12.0
Modkit v0.3.1 [Modkit v0.5.0 for stats submodule]
Dorado v1.0.0 
Samtools v1.18
Gffread v0.12.7
GffCompare v0.12.6
TD2 v1.0.6
GSEApy v1.1.1
Logomaker v0.8.7
GeneStructureTools v1.28.0
bedtools v2.30.0
RiboCode v1.2.11
gppy v0.1.4
pysam 0.22.1
biopython v1.84
scikit-learn v1.4.0
scanpy v1.9.8
statsmodels v0.14.1
CoolBox v0.4.0

Reference genome GRCh38.p14 with GENCODE v46 annotations. The following synthetic unmodified control sequences were used: External RNA Controls Consortium (ERCC) Spike-in mix by Thermo Fisher Scientific and Spike-in RNA Variant (SIRV) by Lexogen. The reference sequences and expected titer levels were found at https://www.lexogen.com/sirvs/download/. 
