# actual reduced normalizer 比較 operator：独立レビュー

2026-11-04。一次証明書
`AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63` を、独立な sparse QQ
多項式演算と typed semantic fingerprints で監査する。primary report は
`79b94d1cf15d3ad05555f78c6cc06afc2c90833cd04134930f3ce72d4a6a65bf`、
evidence は
`376d5fd2b8f2ab8a45651b0374d5b983ebe2df8f8a507108675942e853adbba5`。
source/jet solve、old r/g の再除算、scalar series の展開をしない。

## radical の実九関係と完全な native resolution

直前に独立監査した正規化 E4 の全 radical から

\[
E_{\mathrm{red}}=C/[dL],\qquad
L=(u_1,\ldots,u_5,h,cr,cl,rl),\qquad d=c+r+l
\]

を使う。同じ regular C9 と同じ physical R0/free4 interface を保持する。
実 Q0,Q1,f の保存された全 cofactor から、\(G_8=D_{8\times9}R_9\) の
全八本を literal QQ polynomial identity として再検査する。
proper/core の変数番号、GL4、actual u5 の exact core roots は、受理済み
complete-normalizer interface に直接束縛されている。

Koszul6 on `(u5,h)` と、三軸 ideal `(cr,cl,rl)` の Hilbert–Burch resolution
を tensor する。後者は二列 `(l,-r,0)`, `(l,0,-c)`。
`(u5,h)` は三軸 quotient 上で regular なので tensor complex は exact。
保存された全 basis labels、全 differential を独立に tensor rule から
確認し、全374 column の chain product を QQ[9vars] 上で再検査する。

Betti列は

\[
(1,9,35,77,105,91,49,15,2).
\]

ambient C は regular domain、d は非零なので `ker(dR9)=ker(R9)`。
最初の differential だけ d 倍したものが Ered の完全解像度となる。
全 entries は maximal ideal 内だから C 上では minimal。
C/R0 が同じ free4 algebra であるため restriction は exact、rank は四倍。
R0-minimal とはしない。特に Qrad4×36 の**全 syzygy**は36×140である。

## 同じ actual source と transpose

accepted source unit matrix は `N12=T(kappa)*R12`、逆は同じ
`R12^-1*T(kappa)^-1`。新しい逆行列や未知 regular function を導入しない。

\[
A=N_{12}^{-1}E_8D,\qquad F_{\rm old}=A(dR_9),\qquad
Q_{\rm old}=Q_{\rm rad}\operatorname{Res}_4(A^t).
\]

block(j,i) は \(\mu(A_{ij})\) であり、mu の内部を transpose しない。
独立 checker は全108 A entries、Res4(A^t) の全48列、Qrad36列、complete
Srad140列を、元の actual coefficient interface と基底
`(1,tauA,tauB,tauA*tauB)` から構造的に再構成して比較する。

全53610 typed nodes について C/R0 の型と参照を検査する。内部 child は
既存 local node を参照する一方、old B operator ID と honest C0 scalar ID
は**別 artifact の namespace**である。foreign ID を新53610 nodesの欠番と
誤認したり、同番号の local scalar に置き換えたりしない。
元 r/g ancestry は受理済み foreign regular-function root のまま保持する。

全640 oldB列を十二の actual C-value に group し、A^t を掛けた各九 C-valueと
全36 R0-componentsを独立 semantic fingerprints で照合する。これにより
\(B_{\rm rad}=\operatorname{Res}_4(A^t)B_{\rm old}\) の全640列を確認する。
scalar coefficient の未計算値を評価したと装うのではなく、同じ源に束縛された
有限和積・component extraction の operator identity を監査している。

## exact homology と単位消去

P9×640 は受理済み honest D0 presentation、U4×9 は同じ比較写像。

\[
V=(Q_{\rm rad},-U),\qquad
X=\begin{pmatrix}B_{\rm rad}&S_{\rm rad}\\P&0\end{pmatrix}.
\]

kernel pair `(b,p)` を D0 の class `[p]` に送る。
`p=Pz` なら `(Brad*z,Pz)` を引き、残る `(b-Brad*z,0)` は completeSrad140
で消せる。逆に D0 から Ered へ零となる class には Qrad-lift がある。
従って **全 P640 と全 Srad140** により H1 が比較写像の kernel と等しい。

`U0=e0` のため source列 `-U0=-e0` と target row0 を単位消去する。
他列の row0 を消す source change の逆を X に掛けると、chain恒等式が
除くべき X の pivot row を零にし、他行は変わらない。残りは正確に

\[
R_0^{780}\xrightarrow{X_{\rm red}}R_0^{44}
\xrightarrow{V_{\rm red}}R_0^3.
\]

全 P8×640、V3×44、X の二つの完全 block recipe も照合する。

## Nil(D0) との同定に使う追加 gate

有限 integral・Spec上全射の根拠は
`research/notes/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63.md` §4
と同名 `_REVIEW/report.json` にある。actual二つの有限 normalizer cover を
同じ base 上で tensor したものが E4 である。free4の rank や支持次元だけ
から Spec 全射を推測していない。

D0→Ered も integral・Spec上全射であるため、その kernel は D0 の全 prime
に入り nil。逆向きには Ered が被約なので D0 の全 nilpotent が零へ行く。

\[
\boxed{\operatorname{Nil}(D_0)=\ker(D_0\to E_{\rm red})
       =H_1(X_{\rm red},V_{\rm red}).}
\]

これは actual nilideal の**完全な homology model**である。H1 の生成元を
抽出した cokernel presentation、principality、cube、socle、深さ、finiteT長
を計算したものではない。Nil(E4) からの全生成元の lift も仮定しない。

独立 code:
`code/ambient_B18_radical_comparison_operator_review_20261104_v63.py`。
算術・型 audit、証明・固定 report は同名 `_REVIEW` directory に保存する。
著者のファイルは変更せず、producerや旧640 chain builderは再実行しない。
