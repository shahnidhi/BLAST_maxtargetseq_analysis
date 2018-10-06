# Analyzing max_target_seqs parameter behavior in BLAST command line tool 
## Testing on BLAST 2.6.0

### Difference in BLAST output when using max_target_seqs = 1 vs. not using this parameter
For testing this, please download database file and the query sequences.
```
makeblastdb -in db.fasta  -out db.fasta  -dbtype nucl
blastn -query example.fasta  -db db.fasta -outfmt 6 -max_target_seqs 1 -out example_db_top.blast.out
blastn -query example.fasta  -db db.fasta -outfmt 6  -out example_db_full.blast.out
```
We observe that the first hit for each query sequence is different in both the cases
```
S50_7242:pb_10cov_000M-001M     s_1013:COG0088  84.598  870     33      89      423     1258    1       803     0.0     771
S50_7242:pb_10cov_000M-001M     s_15833:COG0185 87.583  451     8       42      2352    2785    1       420     2.40e-135       479

S208_335:pb_10cov_001M-002M     s_9883:COG0088  86.420  648     35      49      1406    2038    1       610     0.0     660
S208_335:pb_10cov_001M-002M     s_6114:COG0087  85.065  616     32      55      809     1406    3       576     5.90e-164       573

S216_4030:pb_10cov_001M-002M    s_9886:COG0090  82.671  906     46      101     2118    2995    3       825     0.0     701
S216_4030:pb_10cov_001M-002M    s_401:COG0088   82.232  681     40      73      1154    1818    1       616     2.33e-145       512

S2153_228:pb_10cov_022M-023M    s_1068:COG0092  82.841  880     43      98      12227   13084   1       794     0.0     689
S2153_228:pb_10cov_022M-023M    s_7936:COG0087  82.336  702     38      74      8837    9516    1       638     2.27e-150       531

NC_006448.1-15016:454_5cov_045M-050M    s_16779:COG0087 98.540  137     1       1       1       136     492     628     2.46e-65        243
NC_006448.1-15016:454_5cov_045M-050M    s_7827:COG0088  100.000 130     0       0       161     290     1       130     8.87e-65        241

S52_804:pb_10cov_000M-001M      s_14626:COG0090 87.126  769     19      74      911     1656    7       718     0.0     798
S52_804:pb_10cov_000M-001M      s_7931:COG0088  86.025  644     23      62      1       618     195     797     1.42e-180       628

S190_1420:pb_10cov_001M-002M    s_16963:COG0090 84.294  885     49      83      3789    4646    2       823     0.0     782
S190_1420:pb_10cov_001M-002M    s_8639:COG0088  83.639  709     39      72      2762    3452    3       652     3.42e-170       595

S232_2558:pb_10cov_001M-002M    s_10253:COG0201 85.682  1348    63      120     7042    8359    14      1261    0.0     1301
S232_2558:pb_10cov_001M-002M    s_3834:COG0090  84.407  885     51      84      28      883     3       829     0.0     789

S188_1416:pb_10cov_001M-002M    s_17130:COG0201 84.922  1479    47      149     5687    7120    1       1348    0.0     1334
S188_1416:pb_10cov_001M-002M    s_15600:COG0094 86.311  599     15      63      2224    2801    2       554     2.46e-168       590

S1170_2543:pb_10cov_011M-012M   s_1966:COG0201  82.612  1409    94      133     4088    5459    3       1297    0.0     1105
S1170_2543:pb_10cov_011M-012M   s_8667:COG0094  86.038  573     24      53      606     1162    3       535     1.06e-160       564
```
