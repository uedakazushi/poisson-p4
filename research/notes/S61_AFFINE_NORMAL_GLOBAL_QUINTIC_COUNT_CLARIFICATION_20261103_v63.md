# Cofactor項数の表記補足

2026-11-03。凍結済み `S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63` のproducerは `len(Poly(f).terms())` を用いて31/31/52を記録した。SymPyでは零多項式にも一個の零termを返し、三つの恒等式で第五cofactorが零である。従って**実際の非零項数は30/30/51**である。旧noteの「非零項数31/31/52」という語はこの補足で訂正する。

全QQ多項式identity、p-unit minor、全134式のdivision、global quinticおよび5×670 moduleの結論には影響しない。再現可能性のため、凍結producer/data/reportのhashは変更しない。独立reviewは両方の数を区別して保存する。
