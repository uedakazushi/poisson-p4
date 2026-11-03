# 元の被約交叉の導体商：深さ19と重複度1

2026-11-03。generic B18 の形式的 germ での系。
全 normalization conductor の独立受理 report を実行 gate とする。
追加 jet、source 再計算、CAS は不要である。

## 入力と短完全列

受理済み N の実像を A_N、正規化を B_N とし、その導体を IB_N とする。
元の被約交叉 F の正規化導体は (J,T⁴) であり、J⊂IB_N。
全導体の結果は、IB_N/J が実際の共通 G19 を戻した後に正則19変数環上
rank4自由であることを証明している。一方 N の全完成型から

```
C=A_N/IB_N=R/I=C[[g1,...,g19,w,v]]/(wv,v²)
```

が得られる。従って実際の有限 module の短完全列

```
0 → M=IB_N/J → E=A_N/J → C → 0
```

がある。これは環の直和分裂を主張するものではない。

## 深さと次元

M は次元19、深さ19である。C の二変数 core では v≠0 かつ (w,v)v=0 なので
深さ0、次元1。19正則変数を戻すと dim C=20、depth C=19。
よって dim E=20。深さの補題から depth E≥19。

同じ極大イデアルでの局所コホモロジー完全列を使うと

```
H^18(C)=0,    0 → H^19(M) → H^19(E)
```

となる。次元・深さ19の非零 module M について H^19(M)≠0 なので H^19(E)≠0。
従って **depth E=19** である。等しい深さ19の二 module の extension だからと
いうだけでは上界を得られず、この単射を使う必要がある。

## 重複度と導体商

C の被約商は正則 C[[g19,w]]、余分な nil class v の台は次元19。
従って C の Hilbert–Samuel 重複度は1。同じ極大イデアルで短完全列の
Hilbert–Samuel 重複度を加法的に取ると、M は E より低次元なので
e(E)=e(C)=1。

F=A_N×_(C4)C[[T]] で導体は (J,T⁴) だから

```
F/conductor ≅ (A_N/J)×_(C4)C4 ≅ A_N/J.
```

従って全導体商の (次元,深さ,重複度) は cone で **(20,19,1)**。
既証の共通非零 scalar 因子を除けば projective 版は **(19,18,1)** となる。

これは generic B18 の元の被約交叉の導体商だけの結論である。
embedding dimension、環の分裂、元の非被約 D の深さや全環、全域の導体商へは
拡張しない。

## 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_Dred_conductor_quotient_depth_review_20261103_v63.py
```

導体独立受理 report が存在しない、または入力 hash が変わった場合は PASS を生成しない。
二変数 core の小さな monomial Hilbert 値だけを独立に確認し、上記の純粋な
完全列証明を receipt に束縛する。
