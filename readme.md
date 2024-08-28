# tamipami:  a tool for finding the TAM and PAM sites of novel endonucleases

This tools takes a control fastq seqeuncing library and an experimental libray treeted with 
the endonuclease you are trying ti finda TAM /PAM site for.  The two librarues are generated 
from a pool of sequences containing a Gguide site and a mixture of random bases flankning the 
guide region. When a endonuclease recignizes the PAM/TAM site it cuts and those sequenced are 
depleted in the sample library relative to the control.

# Analysis

Tamipami processes the data in several steps.

1. Paied end I'llumina sequences are merges using BBtools employing error correction
2. 
    mergedfile = os.path.join(tempdir, "merged.fastq")
    logging.info("Merging reads.")
    stdout = merge_reads(fastq=fastq, fastq2=fastq2, outfile=mergedfile)
    logging.info(stdout)
    logging.info("Counting PAM/TAM sites.")
    refcount = count_pam(pamlen=pamlen, spacer=spacer, fastq=mergedfile, orientation=orientation)
    logging.info("Runing multiplicitive replacment in case there are zeros")
    refmc = composition.multiplicative_replacement(list(refcount.values()))
    logging.info("Calculating a center log ratio to change from Aitchison geometry to the real number space")
    refclr = composition.clr(refmc)
    ref_n = check_N(list(refcount.values()))
    logging.info("Poisson noise is expected to be {:5.1f} % of total".format(ref_n * 100))
    return refcount, refclr