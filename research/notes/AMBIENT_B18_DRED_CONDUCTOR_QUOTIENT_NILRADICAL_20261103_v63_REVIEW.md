# 導体商の被約環と五次元の横断 nilradical

2026-11-03。受理済みの全導体と導体商深さの系。generic B18 の形式的 germ に限る。
元の非被約 D の nilradical の主張ではない。

G=C[[g1,…,g19]] とし、E=D_red/conductor=A_N/J と書く。
既証の短完全列は

```
0 → M=IB_N/J → E → C=G[[w,v]]/(wv,v²) → 0,
```

で、M は G 上 rank4自由。実際の G19 積分解を取った一次元横断 core では
M は長さ4の ideal である。

一次元局所環の有限長 ideal M は nilpotent である。実際 M は極大イデアル内に
あり、ある n で mⁿM=0 となるので Mⁿ⁺¹=0。従って M は nilradical に入る。
C の nilradical は G*v（v²=wv=0）であり、その被約商は G[[w]]。
ゆえに (v) の E への逆像はちょうど Nil(E) で、

```
0 → M → Nil(E) → G*v → 0,
E_red=G[[w]].
```

この有限 G-module の商 G*v は自由 rank1なので、module として列は分裂する。
従って **Nil(E) は G 上自由 rank5**。ここで環の分裂を主張していない。
被約導体商は正則 cone20次元／projective19次元であり、nilradical の被約台は
共通 G19（projectiveでは18次元）である。

横断 core では Nil(E) は長さ5。商が C[[w]] で零次局所コホモロジーを持たない
ため、H⁰_m(E)=Nil(E)、長さは5である。これも導体商だけの invariant であり、
D_red 自身や元の D の nilpotency index へ読み替えない。

nilpotency の正確な指数、embedding dimension、環の直和分裂、全域の導体商は
求めていない。追加 jet、点評価、CAS は用いない。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_Dred_conductor_quotient_nilradical_review_20261103_v63.py
```

受理済み三 report と全入力 hash を gate とする。
