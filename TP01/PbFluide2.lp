\* Problem: ModeleFluide2 *\

Minimize
 CoutsTotaux: + 10 F(F1,M1,D1) + 2 F(F1,M1,D2) + F(F1,M2,D1)
 + 20 F(F1,M2,D2) + 5 F(F1,M3,D1) + 10 F(F1,M3,D2) + 10 F(F2,M1,D1)
 + 2 F(F2,M1,D2) + F(F2,M2,D1) + 20 F(F2,M2,D2) + 5 F(F2,M3,D1)
 + 10 F(F2,M3,D2) + 110 B(M1,D1) + 110 B(M1,D2) + 90 B(M2,D1)
 + 90 B(M2,D2) + 100 B(M3,D1) + 100 B(M3,D2)

Subject To
 RespectStockMagasins(F1,M1): + F(F1,M1,D1) + F(F1,M1,D2) <= 1
 RespectStockMagasins(F1,M2): + F(F1,M2,D1) + F(F1,M2,D2) <= 2
 RespectStockMagasins(F1,M3): + F(F1,M3,D1) + F(F1,M3,D2) <= 3
 RespectStockMagasins(F2,M1): + F(F2,M1,D1) + F(F2,M1,D2) <= 1
 RespectStockMagasins(F2,M2): + F(F2,M2,D1) + F(F2,M2,D2) <= 3
 RespectStockMagasins(F2,M3): + F(F2,M3,D1) + F(F2,M3,D2) <= 2
 RespectDemandes(F1,D1): + F(F1,M1,D1) + F(F1,M2,D1) + F(F1,M3,D1) = 2
 RespectDemandes(F1,D2): + F(F1,M1,D2) + F(F1,M2,D2) + F(F1,M3,D2) = 1
 RespectDemandes(F2,D1): + F(F2,M1,D1) + F(F2,M2,D1) + F(F2,M3,D1) = -0
 RespectDemandes(F2,D2): + F(F2,M1,D2) + F(F2,M2,D2) + F(F2,M3,D2) = 3
 Contrainte(F1,M1,D1): - F(F1,M1,D1) + 1000 B(M1,D1) >= -0
 Contrainte(F1,M1,D2): - F(F1,M1,D2) + 1000 B(M1,D2) >= -0
 Contrainte(F1,M2,D1): - F(F1,M2,D1) + 1000 B(M2,D1) >= -0
 Contrainte(F1,M2,D2): - F(F1,M2,D2) + 1000 B(M2,D2) >= -0
 Contrainte(F1,M3,D1): - F(F1,M3,D1) + 1000 B(M3,D1) >= -0
 Contrainte(F1,M3,D2): - F(F1,M3,D2) + 1000 B(M3,D2) >= -0
 Contrainte(F2,M1,D1): - F(F2,M1,D1) + 1000 B(M1,D1) >= -0
 Contrainte(F2,M1,D2): - F(F2,M1,D2) + 1000 B(M1,D2) >= -0
 Contrainte(F2,M2,D1): - F(F2,M2,D1) + 1000 B(M2,D1) >= -0
 Contrainte(F2,M2,D2): - F(F2,M2,D2) + 1000 B(M2,D2) >= -0
 Contrainte(F2,M3,D1): - F(F2,M3,D1) + 1000 B(M3,D1) >= -0
 Contrainte(F2,M3,D2): - F(F2,M3,D2) + 1000 B(M3,D2) >= -0

Bounds
 0 <= B(M1,D1) <= 1
 0 <= B(M1,D2) <= 1
 0 <= B(M2,D1) <= 1
 0 <= B(M2,D2) <= 1
 0 <= B(M3,D1) <= 1
 0 <= B(M3,D2) <= 1

Generals
 F(F1,M1,D1)
 F(F1,M1,D2)
 F(F1,M2,D1)
 F(F1,M2,D2)
 F(F1,M3,D1)
 F(F1,M3,D2)
 F(F2,M1,D1)
 F(F2,M1,D2)
 F(F2,M2,D1)
 F(F2,M2,D2)
 F(F2,M3,D1)
 F(F2,M3,D2)
 B(M1,D1)
 B(M1,D2)
 B(M2,D1)
 B(M2,D2)
 B(M3,D1)
 B(M3,D2)

End
