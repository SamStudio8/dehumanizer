from multiprocessing import Process, Queue, Array
import argparse
import ctypes
import sys
import os
import mappy as mp
import numpy as np
import pysam
from datetime import datetime


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

        manifest["references"].append({
            "name": fields[0],
            "path": fields[1],
        })
    manifest_fh.close()
    return manifest


def dh_bam(log, manifest, args):
    dirty_bam = pysam.AlignmentFile(args.dirty)
    clean_bam = pysam.AlignmentFile(args.clean, "wb", template=dirty_bam)
    break_first = not args.nobreak # break on first hit, otherwise we can use this to 'survey' hits to different databases

    aligners = []
    each_dropped = []
    for ref_i, ref_manifest in enumerate(manifest["references"]):
        aligners.append( mp.Aligner(ref_manifest["path"], preset=manifest["preset"]) )
        each_dropped.append(0)
    sys.stderr.write("[INFO] minimap2 aligners ready.\n")


    n_seqs = 0
    n_trash = 0
    n_baddies = 0
    for read in dirty_bam.fetch(until_eof=True):
        n_seqs += 1
        read_is_bad = False

        # Check if the read is trash
        if args.trash_minalen:
            if (read.reference_length/read.query_length)*100.0 < args.trash_minalen:
                read_is_bad = True
                n_trash += 1

        # If the read is not already trash, see if it maps against the references
        if not read_is_bad:
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

        if not read_is_bad:
            clean_bam.write(read)

    sys.stderr.write("[INFO] Took out %d trash sequences, and disposed of %d contaminating sequences\n" % (n_trash, n_baddies))
    log.write("%s\t%d\t%s\t%d\t%d\t-\t%s\n" % (os.path.basename(args.clean), n_seqs, n_baddies, n_trash, n_seqs-(n_baddies+n_trash), "\t".join([str(x) for x in each_dropped])))

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
    for read_i, read_tuple in enumerate(mp.fastx_read(fastx_path)):
        if not flat_dropped[read_i]:
            clean_fq.write("@%s\n%s\n+\n%s\n" % (read_tuple[0], read_tuple[1], read_tuple[2]))
    clean_fq.close()

    each_dropped = list( super_flag_matrix.sum(axis=0) )
    log.write("%s\t%d\t%d\t%d\t%d\t-\t%s\n" % (os.path.basename(clean_fq_p), n_seqs, total_dropped, 0, n_seqs-total_dropped, "\t".join([str(x) for x in each_dropped])))


def cli():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", help="reference manifest")
    parser.add_argument("dirty", help="input dirty file")

    type_p = parser.add_mutually_exclusive_group(required=True)
    type_p.add_argument("--bam", action="store_true")
    type_p.add_argument("--fastx", action="store_true")

    parser.add_argument("--preset", help="mappy aligner preset", required=True)

    parser.add_argument("-o", "--clean", help="output clean file [default -]", default="-")
    parser.add_argument("--log", help="log path [default <dirty>.dehumanizer.log.txt]", default=None)

    parser.add_argument("-t", "--threads", help="number of minimap2 process queues to spawn PER REFERENCE [1]", default=1, type=int)
    parser.add_argument("-n", help="number of reads (prevents having to count)", type=int)
    parser.add_argument("--minid", help="min %proportion of (L-NM)/L to determine a hit [use all hits]", type=float, default=None)
    parser.add_argument("--minlen", help="min %proportion of read aligned to accept a hit [use all hits]", type=float, default=None)

    parser.add_argument("--nobreak", help="dont break on the first database hit [False]", action="store_true", default=False)
    parser.add_argument("--blockrep", help="report progress after a block of N sequences [100000]", default=100000, type=int)

    # Not really the place for it, but whatever
    parser.add_argument("--trash-minalen", help="trash reads whose alignment length is less than this %proportion of their size [keep everything] ignored if not BAM", type=float, default=None)

    args = parser.parse_args()

    #if not args.minid and not args.minlen:
    #    sys.stderr.write("You must set a minimum identity (--minid) and/or minimum length (--minlen).\n")
    #    sys.exit(1)

    if not args.log:
        log = open(args.dirty + ".dehumanizer.log.txt", 'w')
    else:
        log = open(args.log, 'w')

    manifest = load_manifest(args.manifest, args.preset)
    log.write("name\tn_sequences\tn_dropped\tn_trash\tn_saved\t-\t%s\n" % "\t".join([x["name"] for x in manifest["references"]]))

    if args.fastx:
        dh_fastx(log, manifest, args)
    elif args.bam:
        dh_bam(log, manifest, args)

    log.close()

if __name__ == "__main__":
    cli()
