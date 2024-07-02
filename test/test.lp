\ENCODING=ISO-8859-1
\Problem name: 

Minimize
 obj1: 0 Cov_0 + 0 Cov_1 + 0 Cov_2 + 0 Cov_3 + 0.99999 I_neg_0
       + 0.99999 I_neg_1 + 0.99999 I_neg_2 + 0.99999 I_neg_3 + 0.99999 I_pos_0
       + 0.99999 I_pos_1 + 0.99999 I_pos_2 + 0.99999 I_pos_3
       - 0.0000430608174563838 cell1_beta - 0.0000548049560871447 cell2_beta
       - 0.0000660459321155117 cell3_beta - 0.0000790435844410001 cell4_beta
       + [ 0.00000999999999995449 cell1_beta ^2
       + 0.00000999999999995449 cell2_beta ^2
       + 0.00000999999999995449 cell3_beta ^2
       + 0.00000999999999995449 cell4_beta ^2 ] / 2
Subject To
 c1: cell1_beta + cell2_beta + cell3_beta + cell4_beta  = 0
 CovarianceEquation_0: 3.22283861080625 cell1_beta
                       + 4.10410480130889 cell2_beta
                       + 4.94704149293573 cell3_beta
                       + 5.92192651072684 cell4_beta
                       + [ - 0.74931786433838 cell1_beta ^2
                       - 0.74931786433838 cell2_beta ^2
                       - 0.74931786433838 cell3_beta ^2
                       - 0.74931786433838 cell4_beta ^2 ]  = 0
 CovarianceEquation_1: 1.78885438199983 cell1_beta
                       + 2.68328157299975 cell2_beta
                       + 3.57770876399966 cell3_beta
                       + 4.47213595499958 cell4_beta
                       + [ - 0.894427190999916 cell1_beta ^2
                       - 0.894427190999916 cell2_beta ^2
                       - 0.894427190999916 cell3_beta ^2
                       - 0.894427190999916 cell4_beta ^2 ]  = 0
 CovarianceEquation_2: 1.78885438199983 cell1_beta
                       + 2.68328157299975 cell2_beta
                       + 3.57770876399966 cell3_beta
                       + 4.47213595499958 cell4_beta
                       + [ - 0.894427190999916 cell1_beta ^2
                       - 0.894427190999916 cell2_beta ^2
                       - 0.894427190999916 cell3_beta ^2
                       - 0.894427190999916 cell4_beta ^2 ]  = 0
 CovarianceEquation_3: 1.3697963858457 cell1_beta + 2.26422357684561 cell2_beta
                       + 3.15865076784553 cell3_beta
                       + 4.05307795884544 cell4_beta
                       + [ - 0.894427190999916 cell1_beta ^2
                       - 0.894427190999916 cell2_beta ^2
                       - 0.894427190999916 cell3_beta ^2
                       - 0.894427190999916 cell4_beta ^2 ]  = 0
 I_pos_0: I_pos_0 = 0 -> Cov_0 =< 1.33454712291287
 I_neg_0: I_neg_0 = 0 -> Cov_0 => -1.33454712291287
 I_pos_1: I_pos_1 = 0 -> Cov_1 =< 1.11803398874989
 I_neg_1: I_neg_1 = 0 -> Cov_1 => -1.11803398874989
 I_pos_2: I_pos_2 = 0 -> Cov_2 =< 1.11803398874989
 I_neg_2: I_neg_2 = 0 -> Cov_2 => -1.11803398874989
 I_pos_3: I_pos_3 = 0 -> Cov_3 =< 1.11803398874989
 I_neg_3: I_neg_3 = 0 -> Cov_3 => -1.11803398874989
Bounds
      Cov_0 Free
      Cov_1 Free
      Cov_2 Free
      Cov_3 Free
 0 <= I_neg_0 <= 1
 0 <= I_neg_1 <= 1
 0 <= I_neg_2 <= 1
 0 <= I_neg_3 <= 1
 0 <= I_pos_0 <= 1
 0 <= I_pos_1 <= 1
 0 <= I_pos_2 <= 1
 0 <= I_pos_3 <= 1
-3.59827669847801 <= cell1_beta <= 3.59827669847801
-3.59827669847801 <= cell2_beta <= 3.59827669847801
-3.59827669847801 <= cell3_beta <= 3.59827669847801
-3.59827669847801 <= cell4_beta <= 3.59827669847801
Binaries
 I_neg_0  I_neg_1  I_neg_2  I_neg_3  I_pos_0  I_pos_1  I_pos_2  I_pos_3 
End
