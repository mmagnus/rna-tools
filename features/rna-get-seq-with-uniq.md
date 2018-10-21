181021

Imagine that you have a folder with hundreds of models from a few simulations. Now, you want to get preview sequences for one example model from each simulation. So instead of running:

    $ rna_pdb_toolsx.py --get_seq *

and getting a few hundreds of results, you can select uniqueness of your hits, like in the following examples. Show me only sequences of models whose names differ for each other in the first 5 characters ;-)

    [mm] structures$ git:(master) ✗ rna_pdb_toolsx.py --get_seq --uniq '[:5]' *
    # rp13_20569fa1_ALL-000001_AA
    > A:1-71
    GGGUCGUGACUGGCGAACAGGUGGGAAACCACCGGGGAGCGACCCCGGCAUCGAUAGCCGCCCGCCUGGGC
    # rp13cp0016_min.out.1
    > A:1-123
    ACCUUGCGCGACUGGCGAAUCCUGAAGCUGCUUUGAGCGGCUUCGACAUGUGGGGCAGCGACCACGAGGAAGCGCAAGGUUUCUGGCGUCAUGCACAACGGCGCCUGCCGCUCGCCUGGGCAG
    # rp13nc3295_min.out.1
    > A:1-121
    ACCUUGCGCAACUGGCGAAUCCUGGGGCUGCCGCCGGCAGUACCCGGCAGUGGGCGUUGACCACGAGGAAGCGCAAGGUCUCUGGCGUCAUGCACAACGACGCCUGCCGCUCGCCUGGGCA
    # solution
    > A:1-45 57-71
    GGGUCGUGACUGGCGAACAGGUGGGAAACCACCGGGGAGCGACCCGCCGCCCGCCUGGGC
    # zaa_6a68812b_ALL-000001_AA
    > A:1-80
    GCCCGUUCGCGUGACUGGCGCUAGUGAUGGGGAACCAUCGGGGAGCGCGAACCACAUCGCCGCGCGCCUGGGCUCCUCGA
    # zbaa_684ef8ce_ALL-000001_AA
    > A:1-86
    UGAGUUUUCUGCGACUGACGGAUUAUUGCAGAGCACUGCAAGGGAACAGAAAAACUCUUUUUCAGCCGACCGUCUGGGCACACCUG
    # zcp_6537608a_ALL-000001_AA
    > A:1-123
    ACCUUGCGCGACUGGCGAAUCCUGAAGCUGCUUUGAGCGGCUUCGACAUGUGGGGCAGCGACCACGAGGAAGCGCAAGGUUUCUGGCGUCAUGCACAACGGCGCCUGCCGCUCGCCUGGGCAG
    # znc_a1ea6711_ALL-000001_AA
    > A:1-72
    GCUCUCGCGCGACUGGCGACUUUGGAUGGAGCACCAUCGGGGAGCGCGGGAUCGACCGCCGUGCGCCUGGGC

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">Imagine that you have a folder with hundreds of models from a few simulations. Show me only sequences of models whose names differ for each other in the first 5 characters ;-)<br>rna_pdb_toolsx.py --get_seq --uniq &#39;[:5]’ *<a href="https://t.co/nfqrG4ZRmh">https://t.co/nfqrG4ZRmh</a></p>&mdash; rna-pdb-tools (@rna_pdb_tools) <a href="https://twitter.com/rna_pdb_tools/status/1053978809362534400?ref_src=twsrc%5Etfw">October 21, 2018</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 
