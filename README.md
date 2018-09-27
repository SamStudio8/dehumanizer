# dehumanizer
Human DNA where it shouldn't be? Expunge it from your samples with the `dehumanizer`. Just point at a directory of FASTQ that you suspect are contaminated with human DNA, and yield a directory of clean FASTQ files that are free of uninvited humans.

### How does it work?
Reads from your FASTQ files are spewed as quickly as possible into a queue, where a series of [`minimap2`](https://github.com/lh3/minimap2) aligners are waiting. Reads are tested against one (or more) pre-indexed genomic references of your choice, and output into a new FASTQ file if there are no hits. Like most things in bioinformatics, the heavy lifting is performed by things Heng Li wrote: `minimap2`, via [`mappy`](https://pypi.org/project/mappy/), cheers Heng.

### How do I run it?
```python run.py <dirty_dir> <clean_dir>```

### How do I install it?
You don't, it's not production-ready. Thanks for your continued interest in our products.
