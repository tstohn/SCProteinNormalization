\ENCODING=ISO-8859-1
\Problem name: 

Minimize
 obj1: - 3.31889707543919 cell1_beta - 3.48152620336399 cell2_beta
       - 3.60384756819145 cell3_beta - 3.89972888356815 cell4_beta
       + 0 neg_covThres + 0 pos_covThres + [ 0.999 cell1_beta ^2
       + 0.999 cell2_beta ^2 + 0.999 cell3_beta ^2 + 0.999 cell4_beta ^2
       + 0.001 neg_covThres ^2 + 0.001 pos_covThres ^2 ] / 2
Subject To
 c1: cell1_beta + cell2_beta + cell3_beta + cell4_beta  = 0
 cov_inverted#0: 3.30102999566398 cell1_beta + 3.47712125471966 cell2_beta
                 + 3.60205999132796 cell3_beta + 3.90308998699194 cell4_beta
                 - neg_covThres + [ - cell1_beta ^2 - cell2_beta ^2
                 - cell3_beta ^2 - cell4_beta ^2 ] >= 0
 cov_inverted#1: 3.30102999566398 cell1_beta + 3.47712125471966 cell2_beta
                 + 3.60205999132796 cell3_beta + 3.90308998699194 cell4_beta
                 + pos_covThres + [ - cell1_beta ^2 - cell2_beta ^2
                 - cell3_beta ^2 - cell4_beta ^2 ] >= 0
 cov_inverted#2: 2 cell1_beta + 1.74036268949424 cell2_beta
                 + 1.69897000433602 cell3_beta + cell4_beta - neg_covThres
                 + [ - cell1_beta ^2 - cell2_beta ^2 - cell3_beta ^2
                 - cell4_beta ^2 ] >= 0
 cov_inverted#3: 2 cell1_beta + 1.74036268949424 cell2_beta
                 + 1.69897000433602 cell3_beta + cell4_beta + pos_covThres
                 + [ - cell1_beta ^2 - cell2_beta ^2 - cell3_beta ^2
                 - cell4_beta ^2 ] >= 0
Bounds
-0.581413221350318 <= cell1_beta <= 0.581413221350318
-0.581413221350318 <= cell2_beta <= 0.581413221350318
-0.581413221350318 <= cell3_beta <= 0.581413221350318
-0.581413221350318 <= cell4_beta <= 0.581413221350318
 -Inf <= neg_covThres <= 0
End
