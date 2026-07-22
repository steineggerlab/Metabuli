[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/metabuli/README.html) 
![Platform](https://img.shields.io/badge/platform-Mac%20%7C%20Windows%20%7C%20Linux-brightgreen)
# Metabuli

### [Documentation/Manual](https://jaebeom-kim.github.io/metabuli-doc/).

### Metabuli v1.2.0 is released! Please check [major improvements](https://jaebeom-kim.github.io/metabuli-doc/introduction/improvements/).

### New preprint about v1.2.0 is available [here](https://www.biorxiv.org/content/10.1101/2024.06.17.545419v1).



---

***Metabuli*** classifies metagenomic reads by comparing them to reference genomes. You can use Metabuli to profile the taxonomic composition of your samples or to detect specific (pathogenic) species. 

***Sensitive and Specific.*** Metabuli uses a novel k-mer structure, called *metamer*, to analyze both amino acid (AA) and DNA sequences. It leverages AA conservation for sensitive homology detection and DNA mutations for specific differentiation between closely related taxa.

***A laptop is enough.*** Metabuli operates within user-specified RAM limits, allowing it to search any database that fits in storage. A PC with 8 GiB of RAM is sufficient for most analyses.

***A few clicks are enough.*** Metabuli App is now available [here](https://github.com/steineggerlab/Metabuli-App). With just a few clicks, you can run Metabuli and browse the results with Sankey and Krona plots on your PC.

***Short reads, long reads, and contigs.*** Metabuli can classify all types of sequences.


---

  
For more details, please see
[Nature Methods](https://www.nature.com/articles/s41592-024-02273-y), 
[PDF](https://www.nature.com/articles/s41592-024-02273-y.epdf?sharing_token=je_2D5Su0-xVOSjuKSAXF9RgN0jAjWel9jnR3ZoTv0M7gE7NDF_xi_3sW8QdRiwfSJNwqaXItSoeCvr7cvcoQxKLt0oROgWc6urmki9tP80cXEuHPN0D7b4y9y3i8Yv7sZw8MxxhAj7W6p9eZE2zaK3eozdOkXvwADVfso9cXIM%3D), 
[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.05.31.543018v2), or [ISMB 2023 talk](https://www.youtube.com/watch?v=vz2fuRcVwyk).

Please cite: [Kim J, Steinegger M. Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid and DNA. Nature Methods (2024).](https://doi.org/10.1038/s41592-024-02273-y)

---
### 🖥️  [Metabuli App](https://github.com/steineggerlab/Metabuli-App) for Windows, MacOS, and Linux are now available!
> Run taxonomic profiling in just a few clicks and explore results with Sankey and Krona plots.

> Download the app for your OS [here](https://github.com/steineggerlab/Metabuli-App/releases)—no separate Metabuli installation needed.
<p align="center"><img src="https://raw.githubusercontent.com/jaebeom-kim/Metabuli/master/.github/metabuli.jpg" style="max-height: 500px; width: auto;" /></p>




## Reference
- **Taxonomy dump**: [Shen W, Ren H. TaxonKit: a practical and efficient NCBI Taxonomy toolkit. Journal of Genetics and Genomics (2021).](https://doi.org/10.1016/j.jgg.2021.03.006)
- **FASTA format validation**: [Edwards R.A. fasta_validate: a fast and efficient fasta validator written in pure C. Zenodo.](https://doi.org/10.5281/zenodo.2532044) 
- **FASTQ format validation**: [Fonseca N, Manning J. nunofonseca/fastq_utils: 0.25.2. Zenodo.](https://doi.org/10.5281/zenodo.7755574)
