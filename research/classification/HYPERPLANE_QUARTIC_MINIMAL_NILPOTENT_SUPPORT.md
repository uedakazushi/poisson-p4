# 最小冪零・定数項零の四次因子：全支持を MINIMAL_M へ被覆

2026-11-07、課題 `E_reduced_nonnormal`。有限 cofactor 表の純冪零型8、
Jordan型(2,1,1)、translation零を扱う。実退化式は単純Hと既約四次Fの積である。
結論はこの型の**全Poisson源**の既存MINIMAL_M成分への帰属であり、因子の形だけの認識ではない。
成分の唯一性は既存 `legacy/v07/minimal_thickness.tex` 定理minimal-fullを入力とする。
その非被約指数を新たに計算する必要はなく、唯一の被約成分という弱い帰結だけを使う。

## 1. 全二次係数と既約性を保つ標準形

cofactorを `Z=y∂x` とする。`[Z,P]=0` の全係数比較により

    Pxy=a, Pxz=xl+b, Pxw=xm+d,
    Pyz=yl, Pyw=ym, Pzw=f

となる。a,b,d,fはxに依存しない二次以下の多項式、l,mは一次以下である。
定数項や一次項を落としていない。div P=Zから

    l_z+m_w=0, f_w=2l+y l_y, f_z=−2m−y m_y,
    −a_y−b_z−d_w=y.                                  (1)

実Pfaffianは `F=af+y(dl−bm)` で、xに依存しない。
規約 `Z^j=Σ_i ∂i Pij` に対するmodular恒等式は、行列表示で `P grad F=−F Z`。
そのy,z,w成分はFが三変数部分のCasimirであることを表す。

Fは既約四次なのでF|y=0は零でない。y=0では `F0=a0 f0` でありf0も零でない。
modular恒等式のz,w成分から `f0 ∂z F0=f0 ∂w F0=0`。
従ってF0は非零定数で、a0,f0=cも非零定数となる。
(1)をy=0へ制限するとl,mはyで割れ、`l=αy,m=βy`、α,βは定数である。
積分して

    f=c+uy+vy²+3y(αw−βz),   c≠0.                    (2)

α=β=0ならCasimir条件からFはyだけの多項式であり、四次斉次化は二変数式として
複素数体上で分解する。既約四次にはならない。
それ以外では `aα+bβ=1` を選び、`z'=az+bw,w'=αw−βz` とする。
この変更の行列式は1で、Zを保ちl=y,m=0へ移す。α非零・β非零の二主開で全てを覆う。
以下ではprimeを外す。

## 2. 四次以下のCasimirを全て解く

Casimir条件はF_z=0と

    (c+uy+vy²+3yw) F_w−y² F_y=0.                    (3)

である。Fを `Σ_j A_j(y)w^j` と書き、最大のjを取る。
最高w係数から `y² A'_j=3jy A_j` なのでA_jは定数倍のy^(3j)。
全次数四以下という条件は4j≤4を与え、j≤1となる。
従って(3)の全解は

    F=F0+k I,   I=y³w+c y²/2+u y³/3+v y⁴/4.         (4)

独立に再検査できる有限証明書では、四次以下の全15単項式の係数行列を保持した。
完全核は1,I、rank13のminorは非零定数14929920である。
従ってu,vの特別値やcの分母に隠れた例外はない。

Pfaffian恒等式F=af+dy²を使うと、modular恒等式のx成分は

    a F_y+d F_w=k y(af+dy²)=k yF=−yF.

整域性とF≠0からk=−1である。単に同じ因子を持つだけの余分なtensorは含めない。

## 3. 全Poisson fiberの既約源

F= a0 c−I、a0,c非零とする。F_z=0と `F=af+dy²` から
`f a_z+y²d_z=0`。gcd(f,y)=1、次数a_z≤1なのでa_z=d_z=0。
yの零次・一次係数を比較すると、a,dは正確に

    a=a0−a0u y/c−3a0yw/c+A y²,
    d=(F−af)/y².                                    (5)

最後の商は二次以下の多項式である。(1)の残る式から

    Pxz=xy+z(−y−a_y−d_w)+B₂(y,w),                    (6)

B₂は任意の二次以下六係数となる。他の五括弧は(2),(5)とPyz=y²,Pyw=0である。
これで全係数を可逆に復元した。除外したのは既約Fと矛盾する枝だけである。

逆に、任意の `a0,c∈C*` とu,v,A、全六係数B₂に対して、(2),(5),(6)は
div P=y∂x、全四Jacobi零、Pf(P)=a0c−Iを満たす。
計十一変数の源は `(C*)²×A⁹` で既約である。
全式の分母はcのみ。個々の固定源の稠密点や標本のrankから全fiberを推測していない。

## 4. 実際のMINIMAL_M点と全閉包への帰属

旧MINIMAL_Mを座標(X,Y,Z,W)で書き、実際の線形変更

    (x,y,z,w)=(Y/6, X, −W/2, Z)

を施す。このとき(2),(5),(6)の点

    a0=−1/6, c=1/2, u=v=0, A=−1/2, B₂=−y²/2

と全六括弧が一致する。係数・符号・定数項まで確認した。
旧定理minimal-fullは、この点における全Poisson空間の被約既約成分が唯一で、
23次元MINIMAL_M軌道閉包であることを与える。

上記十一変数源の像閉包は既約なので、ある全空間成分Cに含まれる。
その像には旧MINIMAL_M点が含まれるため、唯一性からCはMINIMAL_M成分である。
従って全源がこの既知閉成分へ入る。複素GL変更・平行移動・非零scalarとPGL飽和を戻しても、
各成分は連結群で保たれるので同じ帰属となる。
これはFをMINIMAL_Mの因子へ変えただけでπも同値と推測する議論ではない。

sourceの全Fは、係数wのy³と定数項a0cに共通因子がない原始一次式なので既約である。
四次斉次化にもhの共通因子はない。全源で一般rank四・単純Hの仮定を保持する。

## 5. 証明書と射程

`code/classification_hyperplane_quartic.py --nil-rank1-geometry` は全PDE行列と核を、
`--nil-rank1-support` は全十一変数源・Jacobi・divergence・実Pfと旧点の全六括弧一致を検証する。
報告は `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_RANK1_GEOMETRY_20261107.json` と
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_RANK1_SUPPORT_20261107.json`。
各有限jobは実行前にgateを通過した。gate自体を数学的証明とはしない。
本稿の完全性の証明と既存の唯一成分定理は、producerの恒等式計算とは別の入力である。
独立再構成は `code/classification_minimal_quartic_review.py` と
`CLASSIFICATION_HYPERPLANE_QUARTIC_MINIMAL_NILPOTENT_REVIEW_20261107.json` に保存した。
producer/CASをimportせず全PDE行列、疎Laplaceによるminor14929920、全十一変数源、
全Jacobi・modular恒等式・実Pf、旧Mの全六括弧を再構成した。完全性と唯一成分の適用も
独立読解確認済みである。他のtranslation・Jordan型の帰属はここから従わない。
