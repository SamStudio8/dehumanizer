# dehumanizer
Human DNA where it shouldn't be? Expunge it from your samples with the `dehumanizer`. Just point at a FASTQ or BAM that you suspect is contaminated with human DNA, and `dehumanizer` will rifle through your file, throwing your reads at as many aligning processes as you will allow, to yield a clean file, free of uninvited humans. I am only supporting use of this tool for CLIMB-COVID.

### How does it work?
Reads from your FASTQ or BAM are spewed as quickly as possible into a queue, where a series of hungry [`minimap2`](https://github.com/lh3/minimap2) aligners are waiting. Reads are tested against one (or more) pre-indexed genomic references of your choice, and only enjoy a second life in a new file if there are no hits to anything you don't want to see. Like most things in bioinformatics, the heavy lifting is performed by things Heng Li wrote: `minimap2`, via [`mappy`](https://pypi.org/project/mappy/), cheers Heng.

### How do I install it?

```
pip install git+https://github.com/SamStudio8/dehumanizer
```

This is a work in progress so proceed with caution. Thanks for your continued interest in our products.

### How do I run it?

#### Get yourself a manifest

```
wget https://sam.s3.climb.ac.uk/dehumanizer/20210906/dehuman.20210906.map-ont.mmi
wget https://sam.s3.climb.ac.uk/dehumanizer/20210906/dehuman.20210906.sr.mmi

echo "dh20210906 $(pwd)/dehuman.20210906.map-ont.mmi map-ont" >> manifest.txt
echo "dh20210906 $(pwd)/dehuman.20210906.sr.mmi sr" >> manifest.txt
```

#### Go go go
```
dehumanise <manifest> --fastx <fastq> --preset <minimap2_preset> -o <out.bam> --log <log>
```
```
dehumanise <manifest> --bam <bam> --preset <minimap2_preset> -o <out.bam> --log <log>
```

### Hang on, what's in the index?

The latest index (`20210609`) is built from:

* Humgen (`GCA_000001405.28_GRCh38.p13_genomic.fna`)
* Decoy (`GCA_000786075.2_hs38d1_genomic.fna`)
* HLA (`ipd-imgt-3_45_1.hla_gen.fasta`)

Each index must be built for a specific minimap2 preset, the above manifest will be suitable
for `--preset sr` (short read Illumina) and `--preset map-ont` (long read Nanopore).
If you need something else, you can easily build your own index by downloading
the raw data for the `20210609` index and building it yourself like so:

```
wget https://sam.s3.climb.ac.uk/dehumanizer/20210906/dehuman.20210906.fasta
minimap2 -x <preset> -d dehuman.20210906.<preset>.mmi dehuman.20210906.fasta
```

Then add a line to your `manifest`:

```
dh20210906 $(pwd)/dehuman.20210906.<preset>.mmi <preset>
```

#### md5sums

```
f8c8fcaf162c5f3462dac72be2df9912  dehuman.20210906.fasta
0aca1f615baf43402f47740610dca15c  dehuman.20210906.map-ont.mmi
3eb0a5f89d69e930a31fd6fb2b86cc51  dehuman.20210906.sr.mmi
8c3e0ada224d3e0c5b50f53d2e70b187  GCA_000001405.28_GRCh38.p13_genomic.fna
198bed2665ed0999af0afd8cbc2cb80e  GCA_000786075.2_hs38d1_genomic.fna
f4db8808276606e1365510c2d888e641  ipd-imgt-3_45_1.hla_gen.fasta
```
