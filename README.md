# dehumanizer
Human DNA where it shouldn't be? Expunge it from your samples with the `dehumanizer`. Just point at a FASTQ or BAM that you suspect is contaminated with human DNA, and `dehumanizer` will rifle through your file, throwing your reads at as many aligning processes as you will allow, to yield a clean file, free of uninvited humans.

Currently I am only supporting use of this tool on the CLIMB-COVID platform as I am in the process of making breaking changes to dehumanizer to improve performance.
