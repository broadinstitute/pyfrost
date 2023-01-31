from pyfrost import all_minimizers
#from pyfrostcpp import minhash_iter, Minimizer

seq = "TACGTACTGCTGACGTCTGCAGTCGTACGTCGTACCGTA"
for minimizer, pos in all_minimizers(seq, 11, 5):
    print(minimizer, pos, hash(minimizer))

"""
minhashes = list(minhash_iter(seq, 11, 5))

for pos, minhash_results in enumerate(minhashes):
    kmer = seq[pos:pos+11]
    print(kmer)

    for result in minhash_results:
        minimizer = Minimizer(seq[result.pos:result.pos+5])
        print("-", minimizer, hash(minimizer), result)
"""
