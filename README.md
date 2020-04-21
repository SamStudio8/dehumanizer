# dehumanizer
Human DNA where it shouldn't be? Expunge it from your samples with the `dehumanizer`. Just point at a FASTQ or BAM that you suspect is contaminated with human DNA, and `dehumanizer` will rifle through your file, throwing your reads at as many aligning processes as you will allow, to yield a clean file, free of uninvited humans.

### How does it work?
Reads from your FASTQ or BAM are spewed as quickly as possible into a queue, where a series of hungry [`minimap2`](https://github.com/lh3/minimap2) aligners are waiting. Reads are tested against one (or more) pre-indexed genomic references of your choice, and only enjoy a second life in a new file if there are no hits to anything you don't want to see. Like most things in bioinformatics, the heavy lifting is performed by things Heng Li wrote: `minimap2`, via [`mappy`](https://pypi.org/project/mappy/), cheers Heng.

### How do I run it?

#### Get yourself a manifest
```
wget https://sam.s3.climb.ac.uk/dehumanizer/20200421/GCA_000786075.2_hs38d1_genomic.mmi
wget https://sam.s3.climb.ac.uk/dehumanizer/20200421/GCF_000001405.39_GRCh38.p13_genomic.mmi
wget https://sam.s3.climb.ac.uk/dehumanizer/20200421/ipd-imgt-3_39_0.hla_gen.mmi

echo "hs38d1 $(pwd)/GCA_000786075.2_hs38d1_genomic.mmi" >> manifest.txt
echo "GRCh38 $(pwd)/GCF_000001405.39_GRCh38.p13_genomic.mmi" >> manifest.txt
echo "HLA $(pwd)/ipd-imgt-3_39_0.hla_gen.mmi" >> manifest.txt
```

#### Go go go
```
python run.py manifest.txt <dirty_thing> (--bam|--fastx) -o <clean_thing> --preset <mappy_preset>
```

### How do I install it?

```
pip install git+git@github.com:SamStudio8/dehumanizer.git
```

This is a work in progress so proceed with caution. Thanks for your continued interest in our products.
