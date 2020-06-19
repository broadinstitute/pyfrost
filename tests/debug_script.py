import pyfrost

g = pyfrost.build_from_refs(['tests/data/mccortex.fasta'], k=5, g=3)
for e in g.edges:
    print(str(e))
