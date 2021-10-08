import pyfrost
from pyfrost import links


def test_add_links(mccortex):
    db = links.MemLinkDB()
    db.color = 0

    links.add_links_from_single_sequence(mccortex, db, "TTTCGATGCGATGCGATGCCACG")
    db.save("test.links")

    db2 = links.MemLinkDB.from_file("test.links")
    print(db2.color)


if __name__ == '__main__':
    g = pyfrost.build_from_refs(['tests/data/mccortex.fasta'], k=5, g=3)
    test_add_links(g)
