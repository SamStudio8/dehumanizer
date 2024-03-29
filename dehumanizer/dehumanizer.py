from multiprocessing import Process, Queue, Array
import ctypes
import sys
import os
import mappy as mp
import numpy as np
import pysam
from datetime import datetime
from datetime import date

from . import version

def load_manifest(path, preset):
    manifest = {
        "preset": preset,
        "references": [
        ],
    }
    manifest_fh = open(path)
    for line_i, line in enumerate(manifest_fh):
        fields = line.strip().split() # split on any whitespace if you have whitespace in your ref name you have bigger problems 

        if line[0] == '#':
            continue

        if len(fields) < 3:
            sys.stderr.write("[FAIL] Manifest did not contain a third column mapping a reference to a preset\n")
            sys.stderr.write("       Consult the README to ensure you are using a manifest suitable for dehumaniser >= 0.9.0\n")
            sys.exit(78) # EX_CONFIG

        if fields[2] != preset:
            continue

        manifest["references"].append({
            "name": fields[0],
            "path": fields[1],
        })

    if len(manifest["references"]) == 0:
        sys.stderr.write("[FAIL] Manifest did not contain any references for preset=%s\n" % preset)
        sys.stderr.write("       Consult the README to ensure your manifest is correctly configured and for\n")
        sys.stderr.write("       instructions on how to build your own indexes if needed\n")
        sys.exit(65) # EX_DATAERR
    else:
        sys.stderr.write("[NOTE] Detected %d references in manifest for preset=%s\n" % (len(manifest["references"]), preset))
    manifest_fh.close()
    return manifest


def dh_bam(log, manifest, bad_set, args):
    dirty_bam = pysam.AlignmentFile(args.dirty)

    dirty_header = dirty_bam.header.as_dict()

    pg_date = date.today().strftime("%Y%m%d")
    if args.pg_date:
        if len(args.pg_date) > 0:
            pg_date = args.pg_date

    if "PG" not in dirty_header:
        dirty_header["PG"] = []
    dirty_header["PG"].append({
        "ID": 'dehumanizer.%s' % pg_date,
        "PN": 'dehumanizer',
        "VN": version.__version__,
        "CL": " ".join(sys.argv),
    })
    clean_header = pysam.AlignmentHeader.from_dict(dirty_header)
    clean_bam = pysam.AlignmentFile(args.clean, "wb", header=clean_header)
    break_first = not args.nobreak # break on first hit, otherwise we can use this to 'survey' hits to different databases

    aligners = []
    each_dropped = []
    for ref_i, ref_manifest in enumerate(manifest["references"]):
        sys.stderr.write("[INFO] Init minimap2 aligner: %s (%s)\n" % (ref_manifest["path"], manifest["preset"]))
        aligners.append( mp.Aligner(ref_manifest["path"], preset=manifest["preset"]) )
        each_dropped.append(0)
    sys.stderr.write("[INFO] minimap2 aligners ready.\n")


    n_seqs = 0
    n_good = 0
    n_trash = 0
    n_known = 0
    n_collateral = 0
    n_baddies = 0

    bad_seen = set([])

    if dirty_bam.has_index():
        n_seqs = dirty_bam.mapped + dirty_bam.unmapped
    else:
        # First pass to get the number of sequences without an index
        for read in dirty_bam.fetch(until_eof=True):
            n_seqs += 1
    dirty_bam.close()

    bad_mask = np.zeros(n_seqs, dtype=np.bool)

    # Second pass to establish a bit mask of what to keep
    dirty_bam = pysam.AlignmentFile(args.dirty)
    for r_i, read in enumerate(dirty_bam.fetch(until_eof=True)):

        if not read.query_sequence:
            continue # supp alignment or something, its up to the user to trash these

        read_is_bad = False

        for ref_i, ref_manifest in enumerate(manifest["references"]):
            for hit in aligners[ref_i].map(read.query_sequence):

                if not args.minlen or not args.minid:
                    # a hit is a hit
                    read_is_bad = True
                else:
                    if args.minlen:
                        st = min(hit.q_st, hit.q_en)
                        en = max(hit.q_st, hit.q_en)
                        if ((en - st) / len(read.query_sequence)) * 100 >= args.minlen:
                            read_is_bad = True

                    if args.minid:
                        # http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
                        # "In the PAF format, column 10 divived by column 11 gives the BLAST identity."
                        bscore = hit.mlen / hit.blen
                        if bscore * 100 >= args.minid:
                            read_is_bad = True

                # Criteria satisifed
                if read_is_bad:
                    each_dropped[ref_i] += 1
                    if break_first:
                        break

            else:
                # Continue the outer loop to the next aligner, as no hit was found
                continue
            # Break the aligner loop as we've already break'ed a hit
            break

        if read_is_bad:
            n_baddies += 1


        # Check if the read is trash instead
        if not read_is_bad:
            if args.trash_minalen:
                try:
                    if (read.reference_length/read.query_length)*100.0 < args.trash_minalen:
                        read_is_bad = True
                        n_trash += 1
                except ZeroDivisionError: 
                    read_is_bad = True
                    n_trash += 1

        # Check if the read is on the shitlist
        if not read_is_bad:
            if read.query_name in bad_set:
                read_is_bad = True
                n_known += 1

        if read_is_bad:
            bad_mask[r_i] = 1
            bad_seen.add(read.query_name)

    dirty_bam.close()

    # Third and final pass to write
    dirty_bam = pysam.AlignmentFile(args.dirty)
    for r_i, read in enumerate(dirty_bam.fetch(until_eof=True)):

        # If the read really is good, write it out
        if not bad_mask[r_i]:
            # Finally, check if the QNAME has been tossed out already
            if read.query_name in bad_seen:
                n_collateral += 1
                continue

            n_good += 1
            clean_bam.write(read)

    sys.stderr.write("[INFO] %d sequences in, %d sequences out\n" % (n_seqs, n_good))
    log.write("\t".join([str(x) for x in [
        os.path.basename(args.clean),
        n_seqs,
        n_seqs - n_good,
        n_good,
        n_baddies,
        n_trash,
        n_known,
        n_collateral,
        "-"
        ]] + [str(x) for x in each_dropped]) + '\n')

    dirty_bam.close()
    clean_bam.close()


#TODO FUTURE Would be good to have another layer of multiproc that poured reads from multiple files to any available aligners
#               Need to think carefully about this however; as the mp.Aligner is primed to a particular reference and shared
def dh_fastx(log, manifest, args):

    fastx_path = args.dirty
    break_first = not args.nobreak # break on first hit, otherwise we can use this to 'survey' hits to different databases

    n_seqs = 0
    if args.n:
        n_seqs = args.n
    else:
        for name, seq, qual in mp.fastx_read(fastx_path):
            n_seqs += 1

    sys.stderr.write("[INFO] Preparing memory for flags.\n")
    super_flag_matrix = np.frombuffer(Array(ctypes.c_bool, n_seqs*len(manifest["references"]), lock=False), dtype=ctypes.c_bool)
    super_flag_matrix = super_flag_matrix.reshape(n_seqs, len(manifest["references"]))
    sys.stderr.write("[INFO] Raised %d x %d flags.\n" % (n_seqs, len(manifest["references"])))


    #aligners = []
    #for ref_i, ref_manifest in enumerate(manifest["references"]):
    #    aligners.append([])
    #    sys.stderr.write("[%d/%d] Booting minimap2 aligners.\n" % (ref_i+1, len(manifest["references"])))
    #
    #    for _ in range(args.threads):
    #        aligners[ref_i].append( mp.Aligner(ref_manifest["path"], preset=manifest["preset"]) )

    def map_seqs(work_q, manifest, break_first, block_i):
        aligners = []
        for ref_i, ref_manifest in enumerate(manifest["references"]):
            #sys.stderr.write("[%d:%d/%d] Booting minimap2 aligners.\n" % (block_i, ref_i+1, len(manifest["references"])))
            aligners.append( mp.Aligner(ref_manifest["path"], preset=manifest["preset"]) )
        sys.stderr.write("[%d:] minimap2 aligners ready.\n" % (block_i))

        while True:
            work = work_q.get()
            if work is None:
                return

            for ref_i, ref_manifest in enumerate(manifest["references"]):
                super_flag_matrix[ work["i"] ][ref_i] = 0

            for ref_i, ref_manifest in enumerate(manifest["references"]):
                for hit in aligners[ref_i].map(work["seq"]):

                    if args.minlen:
                        st = min(hit.q_st, hit.q_en)
                        en = max(hit.q_st, hit.q_en)
                        if ((en - st) / len(work["seq"])) * 100 < args.minlen:
                            continue

                    if args.minid:
                        # http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
                        # "In the PAF format, column 10 divived by column 11 gives the BLAST identity."
                        bscore = hit.mlen / hit.blen
                        if bscore * 100 < args.minid:
                            continue

                    # Criteria satisifed
                    super_flag_matrix[ work["i"] ][ref_i] = 1
                    if break_first:
                        break
                else:
                    # Continue the outer loop to the next aligner, as no hit was found
                    continue
                # Break the aligner loop as we've already seen a hit
                break

    sys.stderr.write("[INFO] Counted %d sequences\n" % (n_seqs))
    sys.stderr.write("[INFO] %s\n" % (fastx_path))

    work_queue = Queue(maxsize=args.threads*5000) # Queue N seqs per process
    processes = []

    for _ in range(args.threads):
        p = Process(target=map_seqs, args=(work_queue,manifest,break_first,_))
        processes.append(p)

    for p in processes:
        p.start()

    # Begin adding seqs
    sys.stderr.write("[INFO] Feeding sequences to queue\n")
    start_clock = datetime.now()
    for read_i, read_tuple in enumerate(mp.fastx_read(fastx_path)):
        if read_i % args.blockrep == 0:
            end_clock = datetime.now()
            sys.stderr.write("[NOTE] Queued Read#%d. Last block pushed in %s (%s pseq.)\n" % (read_i, str(end_clock - start_clock), str((end_clock-start_clock)/args.blockrep) ))
            start_clock = datetime.now()
        if args.n:
            if read_i+1 > args.n:
                break

        # Align
        # queue will block until there's room
        work_queue.put({"i": read_i, "seq": read_tuple[1]})

    sys.stderr.write("[INFO] Finished feeding sequences\n")

    # Add sentinels to kill off processes
    sys.stderr.write("[INFO] Wait for queues to empty... be patient\n")
    for _ in range(args.threads):
        work_queue.put(None)

    # Wait for processes to complete work
    for p in processes:
        p.join()



    flat_dropped = ( super_flag_matrix.sum(axis=1) > 0 )
    total_dropped = flat_dropped.sum()
    sys.stderr.write("[INFO] Dropped %d sequences\n" % (flat_dropped.sum()))

    # Now...
    clean_fq_p = args.clean
    if args.clean == "-":
        clean_fq = sys.stdout
    else:
        fp = os.path.basename(fastx_path).split(".")
        clean_fq = open(clean_fq_p, 'w')
        sys.stderr.write("[INFO] Writing FASTX %s\n" % (clean_fq_p))


    # Output FASTX
    n_good = 0
    for read_i, read_tuple in enumerate(mp.fastx_read(fastx_path)):
        if not flat_dropped[read_i]:
            n_good += 1
            if read_tuple[2] is None:
                out_read = ">%s\n%s\n" % (read_tuple[0],
                                          read_tuple[1])
            else:
                out_read = "@%s\n%s\n+\n%s\n" % (read_tuple[0],
                                                 read_tuple[1],
                                                 read_tuple[2])
            clean_fq.write(out_read)
    clean_fq.close()

    each_dropped = list( super_flag_matrix.sum(axis=0) )
    log.write("\t".join([str(x) for x in [
        os.path.basename(clean_fq_p),
        n_seqs,
        n_seqs - n_good,
        n_good,
        total_dropped,
        0,
        0,
        0,
        "-"
        ]] + [str(x) for x in each_dropped]) + '\n')


def cli():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("manifest", help="reference manifest")
    parser.add_argument("dirty", help="input dirty file")

    parser.add_argument("--known", help="new-line delimited list of reads known to be dirty")

    type_p = parser.add_mutually_exclusive_group(required=True)
    type_p.add_argument("--bam", action="store_true")
    type_p.add_argument("--fastx", action="store_true")

    parser.add_argument("--preset", help="mappy aligner preset", required=True)

    parser.add_argument("-o", "--clean", help="output clean file [default -]", default="-")
    parser.add_argument("--log", help="log path [default <dirty>.dehumanizer.log.txt]", default=None)

    parser.add_argument("-t", "--threads", help="number of minimap2 process queues to spawn PER REFERENCE [1]", default=1, type=int)
    parser.add_argument("-n", help="number of reads (prevents having to count)", type=int)
    parser.add_argument("--minid", help="min %%proportion of (L-NM)/L to determine a hit [use all hits]", type=float, default=None)
    parser.add_argument("--minlen", help="min %%proportion of read aligned to accept a hit [use all hits]", type=float, default=None)

    parser.add_argument("--nobreak", help="dont break on the first database hit [False]", action="store_true", default=False)
    parser.add_argument("--blockrep", help="report progress after a block of N sequences [100000]", default=100000, type=int)

    # Not really the place for it, but whatever
    parser.add_argument("--trash-minalen", help="trash reads whose alignment length is less than this %%proportion of their size [keep everything] ignored if not BAM", type=float, default=None)

    parser.add_argument("--pg-date", help="datestamp to insert into BAM PG header [default today in format YYYYMMDD]", default="")

    parser.add_argument("--version", action="version", version="%(prog)s " + version.__version__)

    args = parser.parse_args()

    #if not args.minid and not args.minlen:
    #    sys.stderr.write("You must set a minimum identity (--minid) and/or minimum length (--minlen).\n")
    #    sys.exit(1)

    if not args.log:
        log = open(args.dirty + ".dehumanizer.log.txt", 'w')
    else:
        log = open(args.log, 'w')

    manifest = load_manifest(args.manifest, args.preset)

    log.write("\t".join([
        "name",
        "seqs_in",
        "seqs_total_dropped",
        "seqs_out",
        "n_hits",
        "n_clipped",
        "n_known",
        "n_collateral",
        "-"
    ] + [x["name"] for x in manifest["references"]]) + '\n')

    if args.fastx:
        dh_fastx(log, manifest, args)
    elif args.bam:
        bad_set = set([])
        if args.known:
            bad_set = set([x.strip() for x in open(args.known)])
        dh_bam(log, manifest, bad_set, args)

    log.close()

if __name__ == "__main__":
    cli()
