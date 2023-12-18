\* Problem: ModeleFluide1 *\

Minimize
 CoutsTotaux: + F(F1,M1,D1) + F(F1,M1,D2) + 2 F(F1,M2,D1)
 + 2 F(F1,M2,D2) + 3 F(F1,M3,D1) + 3 F(F1,M3,D2) + F(F2,M1,D1)
 + F(F2,M1,D2) + 3 F(F2,M2,D1) + 3 F(F2,M2,D2) + 2 F(F2,M3,D1)
 + 2 F(F2,M3,D2)

Subject To
 RespectStockMagasins(F1,M1): + F(F1,M1,D1) + F(F1,M1,D2) <= 2.5
 RespectStockMagasins(F1,M2): + F(F1,M2,D1) + F(F1,M2,D2) <= 1
 RespectStockMagasins(F1,M3): + F(F1,M3,D1) + F(F1,M3,D2) <= 2
 RespectStockMagasins(F2,M1): + F(F2,M1,D1) + F(F2,M1,D2) <= 1
 RespectStockMagasins(F2,M2): + F(F2,M2,D1) + F(F2,M2,D2) <= 2
 RespectStockMagasins(F2,M3): + F(F2,M3,D1) + F(F2,M3,D2) <= 1
 RespectDemandes(F1,D1): + F(F1,M1,D1) + F(F1,M2,D1) + F(F1,M3,D1) = 2
 RespectDemandes(F1,D2): + F(F1,M1,D2) + F(F1,M2,D2) + F(F1,M3,D2) = 1
 RespectDemandes(F2,D1): + F(F2,M1,D1) + F(F2,M2,D1) + F(F2,M3,D1) = -0
 RespectDemandes(F2,D2): + F(F2,M1,D2) + F(F2,M2,D2) + F(F2,M3,D2) = 3

End
