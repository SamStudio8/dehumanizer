from multiprocessing import Process, Queue, Array
import argparse
import ctypes
import sys
import os
import pysam
import mappy as mp
import numpy as np
from datetime import datetime

def load_manifest(path):
    manifest = {
        "preset": "",
        "references": [
        ],
    }
    manifest_fh = open(path)
    for line_i, line in enumerate(manifest_fh):
        fields = line.strip().split("\t")

        if line[0] == '#':
            continue

        if line_i == 0:
            manifest["preset"] = fields[0]
        else:
            manifest["references"].append({
                "name": fields[0],
                "path": fields[1],
            })
    manifest_fh.close()
    return manifest



#TODO FUTURE Would be good to have another layer of multiproc that poured reads from multiple files to any available aligners
#               Need to think carefully about this however; as the mp.Aligner is primed to a particular reference and shared
def main(args):

    log = open(os.path.join(args.clean, "dehumanizer.log.txt"), 'a')
    if args.extension:
        seqs = sorted([
            os.path.join(args.dirty, x) for x in os.listdir(args.dirty) if x.endswith(args.extension)
        ])
    else:
        seqs = [args.dirty]

    manifest = load_manifest(args.manifest)
    log.write("fastx\tn_sequences\tn_dropped\tn_saved\t-\t%s\n" % "\t".join([x["name"] for x in manifest["references"]]))

    break_first = not args.nobreak # break on first hit, otherwise we can use this to 'survey' hits to different databases
    sys.stderr.write("Booting minimap2 aligners.\n")
    aligners = []
    for ref_i, ref_manifest in enumerate(manifest["references"]):
        aligners.append( mp.Aligner(ref_manifest["path"], preset=manifest["preset"]) )
    sys.stderr.write("Aligners at the ready.\n")

    for fastx_i, fastx_path in enumerate(seqs):
        sys.stderr.write("Counting sequences...\n")
        sys.stderr.write("[%d/%d] %s\n" % (fastx_i+1, len(seqs), str(fastx_path)))
        n_seqs = 0
        for name, seq, qual in mp.fastx_read(fastx_path):
            n_seqs += 1

        sys.stderr.write("Preparing memory for flags.\n")
        super_flag_matrix = np.frombuffer(Array(ctypes.c_bool, n_seqs*len(aligners), lock=False), dtype=ctypes.c_bool)
        super_flag_matrix = super_flag_matrix.reshape(n_seqs, len(aligners))

        def map_seqs(work_q, ref_i, aligners, break_first):
            #TODO Do we need to sleep?
            while True:
                work = work_q.get()
                if work is None:
                    return

                for a_i, a in enumerate(aligners):
                    super_flag_matrix[ work["i"] ][a_i] = 0
                    for hit in a.map(work["seq"]):
                        super_flag_matrix[ work["i"] ][a_i] = 1

                        if break_first:
                            break
                    else:
                        # Continue the outer loop to the next aligner, as no hit was found
                        continue
                    # Break the aligner loop as we've already seen a hit
                    break

        sys.stderr.write("[%d/%d] Counted %d sequences\n" % (fastx_i+1, len(seqs), n_seqs))
        sys.stderr.write("[%d/%d] %s\n" % (fastx_i+1, len(seqs), fastx_path))

        work_queue = Queue(maxsize=args.threads*1000) # Queue up to 1000 seqs per process
        processes = []

        for _ in range(args.threads):
            p = Process(target=map_seqs, args=(work_queue,ref_i,aligners,break_first))
            processes.append(p)

        for p in processes:
            p.start()

        # Begin adding seqs
        sys.stderr.write("Feeding sequences to queue\n")

        start_clock = datetime.now()
        for read_i, read_tuple in enumerate(mp.fastx_read(fastx_path)):
            if read_i % args.blockrep == 0:
                end_clock = datetime.now()
                sys.stderr.write("Queued Read#%d. Last block pushed in %s (%s pseq.)\n" % (read_i, str(end_clock - start_clock), str((end_clock-start_clock)/args.blockrep) ))
                start_clock = datetime.now()

            # Align
            # queue will block until there's room
            work_queue.put({"i": read_i, "seq": read_tuple[1]})

        sys.stderr.write("Finished feeding sequences\n")

        # Add sentinels to kill off processes
        for _ in range(args.threads):
            work_queue.put(None)

        sys.stderr.write("Wait for queues to empty... be patient\n")
        # Wait for processes to complete work
        for p in processes:
            p.join()


        flat_dropped = ( super_flag_matrix.sum(axis=1) > 0 )
        total_dropped = flat_dropped.sum()
        sys.stderr.write("[%d/%d] Dropped %d sequences\n" % (fastx_i+1, len(seqs), flat_dropped.sum()))

        # Now...
        fp = os.path.basename(fastx_path).split(".")
        clean_fq_p = os.path.join(args.clean, "%s.%s.%s" % (".".join(fp[:-1]), "dehumanizer.clean", fp[-1]))
        clean_fq = open(clean_fq_p, 'w')

        if args.keepdirty:
            dirty_fq_p = os.path.join(args.clean, "%s.%s.%s" % (".".join(fp[:-1]), "dehumanizer.dirty", fp[-1]))
            dirty_fq = open(dirty_fq_p, 'w') #TODO gzip?


        # Output FASTX
        sys.stderr.write("[%d/%d] Writing FASTX %s\n" % (fastx_i+1, len(seqs), clean_fq_p))
        if args.keepdirty:
            sys.stderr.write("[%d/%d] Writing FASTX %s\n" % (fastx_i+1, len(seqs), dirty_fq_p))
        for read_i, read in enumerate(pysam.FastxFile(fastx_path)): #TODO this is garbage for gz
            if not flat_dropped[read_i]:
                clean_fq.write(str(read)+'\n')
            else:
                if args.keepdirty:
                    dirty_fq.write(str(read)+'\n')
        clean_fq.close()
        if args.keepdirty:
            dirty_fq.close()

        each_dropped = list( super_flag_matrix.sum(axis=0) )
        log.write("%s\t%d\t%d\t%d\t-\t%s\n" % (os.path.basename(clean_fq_p), n_seqs, total_dropped, n_seqs-total_dropped, "\t".join([str(x) for x in each_dropped])))
    log.close()


def cli():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", help="reference manifest")
    parser.add_argument("dirty", help="dirty human file or folder, if folder you must specify --extension")
    parser.add_argument("clean", help="clean non-human folder")

    parser.add_argument("-t", "--threads", help="number of minimap2 process queues to spawn [1]", default=1, type=int)
    #parser.add_argument("--minid", help="minimum sequence id to determine a hit [minimap2 default]")
    #parser.add_argument("--minlen", help="minimum alignment length to determine a hit [minimap2 default]")

    parser.add_argument("--nobreak", help="dont break on the first database hit [False]", action="store_true", default=False)
    parser.add_argument("--extension", help="dirty is a folder, containing files to clean with a particular extension", default=None)
    parser.add_argument("--blockrep", help="report progress after a block of N sequences [100000]", default=100000, type=int)
    parser.add_argument("--keepdirty", help="make a file of the dirty sequences in the same location as clean [False]", action="store_true", default=False)

    main(parser.parse_args())


if __name__ == "__main__":
    cli()
