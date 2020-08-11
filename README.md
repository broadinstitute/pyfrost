Pyfrost
=======

A Python library for creating and analyzing compacted colored de Bruijn Graphs powered by 
[Bifrost](https://github.com/pmelsted/Bifrost) with a [NetworkX](https://github.com/networkx/networkx) compatible API. 

This library is still in an early development stage, and the API is still subject to change. Furthermore, not all
functions from NetworkX are implemented yet.

Requirements
------------

* Python >= 3.6
* C++14 compatible compiler (GCC >5, Clang >3.4)

Installation
------------

Clone the repository, including third-party dependencies as git-submodules:

```bash
git clone --recurse-submodules https://github.com/broadinstitute/pyfrost
```

Install the library:

```bash
cd pyfrost
pip install .
```

Usage 
-----

### Building and loading graphs

Build a graph from a list of references:

```python
import pyfrost
g = pyfrost.build_from_refs(['data/ref1.fa', 'data/ref2.fa'])
```

Build a graph from reads:

```python
g = pyfrost.build_from_samples([
    'data/sample1.1.fq.gz', 'data/sample1.2.fq.gz',
    'data/sample2.1.fq.gz', 'data/sample2.2.fq.gz'
])
```

Build a graph from both references and reads:

```python
g = pyfrost.build(
    ['data/ref1.fa', 'data/ref2.fa'],
    ['data/sample1.1.fq.gz', 'data/sample1.2.fq.gz']
)
```

All calls above also accept optional parameters `k` and `g` for the k-mer and
minimizer size respectively.

Loading a *colored* graph created earlier with Bifrost:

```python
g = pyfrost.load('graph.gfa', 'graph.bfg_colors')
```

Access some graph metadata:

```pycon
>>> g.graph
{'k': 31, 'color_names': ['data/ref1.fa', 'data.ref2.fa']}
```

You can set custom metadata too:

```pycon
>>> g.graph['custom_attr'] = "test"
```

### Access unitigs and k-mers

Iterate over all unitigs (nodes) in the graph. Nodes are keyed by the first
k-mer of the unitig.

```pycon
>>> list(g.nodes)
[<Kmer 'TCGAA'>,
 <Kmer 'CCACG'>,
 <Kmer 'CGATG'>,
 <Kmer 'ATGCG'>,
 <Kmer 'GTGGC'>,
 <Kmer 'ATCGA'>]

>>> len(g.nodes)
6
```

Access unitig sequence and node metadata

```pycon
>>> str(g.nodes["TCGAA"])
TCGAAATCAGT

>>> g.nodes["TCGAA"]['unitig_sequence']
TCGAAATCAGT

>>> g.nodes["TCGAA"]['head']
<Kmer 'TCGAA'>

>>> g.nodes["TCGAA"]['tail']
<Kmer 'TCAGT'>

>>> g.nodes["TCGAA"]['strand']
Strand.FORWARD

>>> for k, v in g.nodes['TCGAA'].items():
   ...:     print(k, "-", v)
   ...:
colors - <bifrost_python.UnitigColors object at 0x10d79d5b0>
unitig_sequence - TCGAAATCAGT
tail - TCAGT
length - 7
strand - Strand.FORWARD
mapped_sequence - TCGAAATCAGT
pos - 0
is_full_mapping - True
unitig_length - 7
head - TCGAA
```

`length` and `unitig_length` are in terms of number of k-mers, not nucleotides.
Mapped sequence and length will be discussed later.

Iterate over nodes and its neighbors:

```pycon
>>> for n, neighbors in g.adj.items():
   ...:     print("Current node:", str(g.nodes[n]))
   ...:     for nbr in neighbors:
   ...:         print("-", str(nbr))
   ...:
Current node: TCGAAATCAGT
Current node: CCACGGTGG
- GTGGC
Current node: CGATGC
- ATGCC
- ATGCG
Current node: ATGCGAT
- CGATG
Current node: GTGGCAT
- GCATC
Current node: ATCGA
- TCGAA
- TCGAT
```

Find a specific k-mer (not necessarily the head of a unitig):

```pycon
>>> mapping = g.find('AAATC')
>>> for k, v in mapping.items():
    ...:     print(k, "-", v)
    ...:
colors - <bifrost_python.UnitigColors object at 0x10d79d1b0>
unitig_sequence - TCGAAATCAGT
tail - TCAGT
length - 1
strand - Strand.FORWARD
mapped_sequence - AAATC
pos - 3
is_full_mapping - False
unitig_length - 7
head - TCGAA
```

As you can see, this k-mer is located on the same unitig used in an earlier
example. The metadata, however, now shows different values for
`mapped_sequence` and `length`, and `is_full_mapping` is now False, because
this k-mer doesn't represent a whole unitig. `head` and `tail` still refer to
the head and tail k-mer of the whole unitig.

### Access unitig colors

Bifrost builds *colored* compacted De Bruijn graphs, and keeps track which k-mers are observed in which references
/samples. To obtain which colors are associated with a unitig or a k-mer, access 'colors' in the unitig metadata dict.

```python
for n, data in g.nodes(data=True):
    for c in data['colors']:
        print("Node", n, "has color", c)

```

**Note:** By default `data['colors']` iterates over *all* colors associated with any of the k-mers of that unitig. 
However, it is possible that not all k-mers share the same colors. To access the colors of a specific k-mer of the
unitig use `g.find()` or index:
 
```python
for n, data in g.nodes(data=True):
    print("Colors of first k-mer in the unitig:", set(data['colors'][0]))
```

### Build the DNA sequence for a path through the graph:

```pycon
>>> from pyfrost import Kmer, path_sequence, path_nucleotide_length

# Three nodes: ACTGATTTCGA, TCGAT, CGATGC
>>> path_sequence(g, [Kmer('ACTGA'), Kmer('TCGAT'), Kmer('CGATG')])
ACTGATTTCGATGC
>>> path_nucleotide_length(g, [Kmer('ACTGA'), Kmer('TCGAT'), Kmer('CGATG')])
14
```

### K-mer counter

Pyfrost includes a separate k-mer counter. It's still pretty unoptimzed and slow though, but it works.


```pycon
>>> from pyfrost import KmerCounter
>>> counter = KmerCounter(31).count_kmers("ACTGCTAGCTAGCTACGTACGTACGATCGTACATGCATGC")
>>> counter["ACTGCTAGCTAGCTACGTACGTACGATCGT"]
1

# You can k-merize (gzipped) FASTA/FASTQ files directly
>>> counter = KmerCounter(31).count_kmers_files(["data/sample1.fq.gz", "data/sample.2.fq.gz"])

# Save counts for later use
>>> counter.save("sampe.counts")

# Load k-mer counts
>>> counter = KmerCounter.from_file("sample.counts")
```

