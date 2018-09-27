from multiprocessing import Process, Queue, Array
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
        "reads": sorted([os.path.join(DIRTY_DIR, x) for x in os.listdir( os.path.dirname(DIRTY_DIR) ) if "fastq" in x]),
    }
    for line_i, line in enumerate(open(path)):
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
    return manifest



#TODO FUTURE Would be good to have another layer of multiproc that poured reads from multiple files to any available aligners
#               Need to think carefully about this however; as the mp.Aligner is primed to a particular reference and shared
def main(DIRTY_DIR, CLEAN_DIR, log):
    manifest = load_manifest(os.path.join(DIRTY_DIR, "manifest.txt"))
    log.write("fastq\tn_reads\tn_dropped\tn_saved\t-\t%s\n" % "\t".join([x["name"] for x in manifest["references"]]))

    break_first = True # break on first hit, otherwise we can use this to 'survey' hits to different databases
    sys.stderr.write("Booting minimap2 aligners.\n")
    aligners = []
    for ref_i, ref_manifest in enumerate(manifest["references"]):
        aligners.append( mp.Aligner(ref_manifest["path"], preset=manifest["preset"]) )
    sys.stderr.write("Aligners at the ready.\n")

    for fastq_i, fastq_path in enumerate(manifest["reads"]):
        sys.stderr.write("Counting reads...\n")
        n_reads = 0
        for name, seq, qual in mp.fastx_read(fastq_path):
            n_reads += 1

        sys.stderr.write("Preparing memory for flags.\n")
        super_flag_matrix = np.frombuffer(Array(ctypes.c_bool, n_reads*len(aligners), lock=False), dtype=ctypes.c_bool)
        super_flag_matrix = super_flag_matrix.reshape(n_reads, len(aligners))

        def map_reads(work_q, ref_i, aligners, break_first):
            #TODO Do we need to sleep?
            while True:
                work = work_q.get()
                if work is None:
                    return
                
                for a_i, a in enumerate(aligners):
                    super_flag_matrix[ work["i"] ][a_i] = 0
                    for hit in a.map(work["seq"]):
                        #print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
                        super_flag_matrix[ work["i"] ][a_i] = 1

                        if break_first:
                            break
                    else:
                        # Continue the outer loop to the next aligner, as no hit was found
                        continue
                    # Break the aligner loop as we've already seen a hit
                    break

        sys.stderr.write("[%d/%d] Preparing %d reads from FASTQ\n" % (fastq_i+1, len(manifest["reads"]), n_reads))
        sys.stderr.write("[%d/%d] %s\n" % (fastq_i+1, len(manifest["reads"]), fastq_path))

        work_queue = Queue(maxsize=PROCESSES*100) # Queue up to 100 reads per process
        processes = []

        for _ in range(PROCESSES):
            p = Process(target=map_reads, args=(work_queue,ref_i,aligners,break_first))
            processes.append(p)

        for p in processes:
            p.start()

        # Begin adding reads
        sys.stderr.write("Feeding reads to queue\n")

        start_clock = datetime.now()
        for read_i, read_tuple in enumerate(mp.fastx_read(fastq_path)):
            if read_i % BLOCK_SIZE == 0:
                end_clock = datetime.now()
                sys.stderr.write("Queued Read#%d. Last block completed in %s (%s pseq.)\n" % (read_i, str(end_clock - start_clock), str((end_clock-start_clock)/BLOCK_SIZE) ))
                start_clock = datetime.now()

            # Align
            # queue will block until there's room
            work_queue.put({"i": read_i, "seq": read_tuple[1]})

        sys.stderr.write("Finished feeding sequences\n")

        # Add sentinels to kill off processes
        for _ in range(PROCESSES):
            work_queue.put(None)

        sys.stderr.write("wait for queues to empty... be patient\n")
        # Wait for processes to complete work
        for p in processes:
            p.join()


        flat_dropped = ( super_flag_matrix.sum(axis=1) > 0 )
        total_dropped = flat_dropped.sum()
        sys.stderr.write("[%d/%d] Dropped %d reads from FASTQ\n" % (fastq_i+1, len(manifest["reads"]), flat_dropped.sum()))

        # Now...
        clean_fq_p = os.path.join(CLEAN_DIR, os.path.basename(fastq_path).replace(".gz", ""))
        clean_fq = open(clean_fq_p, 'w') #TODO gzip?
        sys.stderr.write("[%d/%d] Writing FASTQ %s\n" % (fastq_i+1, len(manifest["reads"]), clean_fq_p))
        for read_i, read in enumerate(pysam.FastxFile(fastq_path)):
            if not flat_dropped[read_i]:
                clean_fq.write(str(read)+'\n')
        clean_fq.close()

        each_dropped = list( super_flag_matrix.sum(axis=0) )
        log.write("%s\t%d\t%d\t%d\t-\t%s\n" % (os.path.basename(clean_fq_p), n_reads, total_dropped, n_reads-total_dropped, "\t".join([str(x) for x in each_dropped])))
        #TODO Report Time?


if __name__ == "__main__":
    PROCESSES = 12
    BLOCK_SIZE = 100000
    DIRTY_DIR = sys.argv[1]
    CLEAN_DIR = sys.argv[2]
    log = open(os.path.join(CLEAN_DIR, "manifest.txt"), 'w')

    main(DIRTY_DIR, CLEAN_DIR, log)

