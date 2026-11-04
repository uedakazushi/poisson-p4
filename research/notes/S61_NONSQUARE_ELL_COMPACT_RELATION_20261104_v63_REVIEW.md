# compact cubic relation / 実 source kernel38 の独立監査

2026-11-04。対象は凍結済み `S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63`。
generic canonical11 rank10 は先の BORDERS REVIEW で受理済み。本監査は、未知 pivot
分母を含まない新しい polynomial source vector と、既知37方向からの独立性を別々に検証する。

## 1. 元 source に直接戻した literal identities

元の六 quad を source sparse coefficients から再構成し、保存された六行174係数と照合する。
元 normalized cubic43,49 は、actual135 cubic20,61、凍結した remaining57 source の行42,47。
既存 row-basis metadata と actual labels を確認する。319 行定義は accepted 354→319 certificate
の retained rowsから独立に再構成し、29/91/239 weighted monomial 列を全て照合する。
quad/cofactor の source frame と genuine cubic の source frame は既存 frame REVIEW のものを使う。

新 L-frame の G6−G1 は全40係数中、最初の五つだけ非零で、残り35は literal zero。
その五係数は L²*(L,x0,y1,y2,y3) に当たる。実際に L=ΣLj*旧xj と戻して E8 を再構成する。
両 genuine 行の clearing は同じ den³/C=−64000C²=delta であり、元八変数では

    delta*(G61−G20)−8100*q1*E = Z*A+W*B,
    Z=3v0−v1+v6, W=v3

の **全91 coefficient** を独立な polynomial arithmetic で照合する。
A,B の保存値を入力として恒等式を直接検証し、author の quotient 抽出算法を再実行しない。
ell は別の transformed-coordinate source を仮定せず、実六 quad から
−(de+8)q1+8q2+dq3 の W係数として独立に回収する。その ell が保存された ell8 と一致し、

    a*ell*Z=(10a+3e²−3e)q1−3(e−1)q3−10a*q4,
    ell*W=−(de+8)q1+8q2+dq3

を満たすことを全58係数で検証する。従って cubic identity に a*ell を掛けた式は、
四つの quad 行の quadratic cofactors と二つの genuine cubic の linear cofactorsだけで
literal な weight4 source-kernel vector を与える。

保存された六 quad cofactor 行全174係数、cubic43/49 の cofactors ±a*delta*ell、
全33 quartic cofactors zero を検証する。さらにこの recipe を使った和だけで済ませず、
保存 **319×1 vector の各行定義** から元 source polynomial を引き当てて直接 dot product を作り、
weight4 の **全239係数が零**であることを照合する。quartic source は cofactor が零なので展開不要。
vector は123非零 entry・19694 core terms。全て QQ[a,e,d] 多項式で、
分母は非零有理数のみであり、generic border pivot の未知因子は反転していない。

## 2. 既知37方向との独立性を証明する別 gate

FULLU=D(adFBCJ) 上の既知 kernel は quadratic22 と direct Bianchi15 の直和で rank37。
そのうち quadratic22 は cubic/quartic source coordinate を持たない。
accepted I11 quartic pivotsにより、quartic-zero の D15 組合せは残りの四方向 R0,…,R3 へ移る。
この statement と actual319 への source transport は既存二つの独立 REVIEW の exact hashes に束縛する。

正規化された residual4 は、固定 cubic source coordinate

    p0=176 (G18*v2), p1=177 (G18*v3),
    p2=237 (G86*v0), p3=181 (G20*v0)

に I4 を持つ。本レビューでもこの16 scalar identityを確認する。
新 vector は p0,p1,p2 で零。S(v)=v[182]+v[186]、すなわち G20*v1 と G20*v6 の和とすると、
最後の residual R3 は S(R3)=−1+1=0。一方、新 vector N は

    S(N)=−1024000*a*C³*J/243

を持つ。この式を full319 vector と元六 quadからの ell により直接確認する。
これは FULLU の既知 unit で、未知 minor の非零性から単元と推論したものではない。
I4 と S の計24 selected rational entries の分母も、非零有理数と既知 a,F,B,C の積だけであると確認する。

従って c*N が旧37の組合せなら、quartic projectionを消して residual4へ移り、
p0,p1,p2を見ると最初の三係数が零になる。残りに S を適用すると c*S(N)=0、よって c=0。
これは各 FULLU fiber と base module 上で成立する。
より明示的には、まず I11 quartic projector で quartic component を消し、
その後 S から Σ(i=0..2) S(Ri)*coordinate[pi] を引くと、旧37を消して N を既知 unitへ送る functional となる。
既存の37方向の独立性と合わせ、**全 FULLU 上で38方向が独立**である。

## 3. scope

この結論は generic bordered rank10 だけからは出ない。新 polynomial vector の literal239 identity と
既知 unit separator が追加の根拠である。旧22の scope に従って combined38 は FULLU 上に限る。
d=0/J=0やaFBC全体へ拡張しない。
残る kernel target80 まで42方向。全 kernel exhaustion、full quartic239、全点 cubic57は未証明。
新しい free281 finite circuit の構成は別 packet であり、本 REVIEW は実行も受理もしていない。

再実行（新しい空 output directory）：

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_ell_compact_relation_review_20261104_v63.py --output work/marked-square-v63/ell-compact-review-REPLAY
```
