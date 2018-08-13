dentist
=========

Close assembly gaps using long-reads with focus on correctness.


Building
--------

```sh
git clone https://github.com/a-ludi/dentist.git
cd dentist
dub build
```


Developer Notes
---------------

```sh
# run the program for development
dub run -- ARGS...

# run tests
dub test

# clean up
dub clean
```


Functional Description
-------------------

TODO

### Difficulties

#### Low Complexity Regions/Repeats

1. Low Complexity Regions
    ```fasta
    ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
    tcctgtgacatgtcactgttgcggtcgaaccggtcgtgcaatccgac
    gtcccaatgcccgccgcattaacggtagccatAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAcgc
    atcaccgatcggggtcggtaataaaaggacaaagttagtgttggcca
    cgaacttctcacgaataagttccctggttttgcgagggaatgcatct
    gctaggcgtcactggacacagtgggaaagctgccgggggcga
    ```
2. Low Complexity Regions
   ```fasta
   ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
   tcctgtgacatgtcAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGT
   AGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAG
   TAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTccacgagctggag
   cctaaaacaattccatgagactggtctaggttacgcagtgtagccgc
   atcaccgatcggggtcggtaataaaaggacaaagttagtgttggcca
   cgaacttctcacgaataagttccctggttttgcgagggaatgcatct
   gctaggcgtcactggacacagtgggaaagctgccgggggcga
   ```
3. Tandem Repeats
   ```fasta
   ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
   tcctgtgacatgtcactgttgcggtcgTTGTGTATCGAGGGGACTTA
   TAGTGCTCCTGTTTGTGTATCGAGGGGACTTATAGTGCTCCTGTTTG
   TGTATCGAGGGGACTTATAGTGCTCCTGTTTGTGTATCGAGGGGACT
   TATAGTGCTCCTGTTTGTGTATCGAGGGGACTTATAGTGCTCCTGTT
   TGTGTATCGAGGGGACTTATAGTGCTCCTGTTTGTGTATCGAGGGGA
   CTTATAGTGCTCCTGTaagttccctggttttgcgagggaatgcatct
   gctaggcgtcactggacacagtgggaaagctgccgggggcga
   ```
4. Transposable elements
   ```fasta
   ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
   tcctgtGACATGTCACTGTTGCGGTCGAACCGGTCGTGcaatccgac
   gtcccaatgcccgccgcattaacggtagccataGACATGTCACTGTT
   GCGGTCGAACCGGTCGTGtagcgcgacaaaaaccccacgagctggag
   cctaaaacaattccatgagactggtctaggGACATGTCACTGTTGCG
   GTCGAACCGGTCGTGtcgtaaaggtctgtcatagtttgtgtgtgtga
   gcggaagtataaacgaaaagaggaccagaaaaGACATGTCACTGTTG
   CGGTCGAACCGGTCGTGacagtgggaaagctgccgggggcga
   ```


License
-------

This project is licensed under [MIT License](./LICENSE).
