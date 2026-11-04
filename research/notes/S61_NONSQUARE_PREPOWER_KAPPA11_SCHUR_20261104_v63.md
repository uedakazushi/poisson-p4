# κ・tail2による未展開actual11pivotと4×11 Schur

2026-11-04。新しいscalar coefficient/minor/rank計算をしない。
受理対象はprojected prePOWER15×22であり、raw H7へのeta/HZ/HW/Cdelta liftは別gate。

## 実8minorの多項式恒等式

source8順はcurrent15 `[0,1,3,6,2,4,5,8]`。pure六列は元32
`[22,11,13,25,15,2]`、最後にheavy `[27,31]` を置く。列2はB倍し、二heavy列は
FULLU既知単元L=aFC²、L*gammaで割る。この実行列のdetをD8tildeとする。
四pair-only rowsの正規化heavy entriesはc_i*N_i,jである。

κはこの四行について `[H0,H1,pure22]` のsigned cofactor、κ3はrow current15=8の係数。
独立監査済みの六pure列rank一から、κは最初五pure列と二heavy列を全点恒等式で消す。
最後のsource rowをκ結合へ変えるrow operatorのdetはκ3。最後の行の唯一の残存成分は
列index5の T=BκM2 なので、符号(7+5)は正で、

    κ3*D8tilde = T*D7.

D7はsource順 `[P4,pair1,pair3,u4]`、pure五列と正規化heavy二列の7minorである。
既存rho_i=c_i*LとBhat定義を逐語的に使うと

    D7=c1*c3*c4*Saux,
    κ3*D8tilde=(c1*c3*c4)*Saux*T.                 (★)

7minorの符号はheavy pair二行を四pure行の前へ動かす8回、heavy二列を五pure列の前へ
動かす10回の交換がともに偶数である。初めD(delta2)でheavy pairを消して導き、分母を
既知単元で払ってintegral domainへ戻せるため、(★)自体にpの逆は不要。
Saux非零とT非零は既存二つの全QQ source計算の結果である。

## κ3を追加で反転しない

D(Saux*T)では(★)右辺が単元だからκ3,D8tildeは両方とも単元となる。明示的に

    D8tilde^-1 = κ3 / ((c1*c3*c4)*Saux*T).

従ってκ3の零集合を別chartとして捨てない。Adj8と上のinverseを有限operatorにし、
D8の巨大な係数展開を計算済みと呼ばない。

## 三つのB sourceを加える

既存B3 sources current15 `[7,13,10]`、heavy targets元32 `[28,30,29]`を後へ足す。
前八sourceのこの三列はliteral zero。11×11 blockは

    K11=[[A8,0],[C38,B3actual]],
    det K11=D8tilde*(rho6*rho12*rho9)*Delta407.

B3actual=diag(rho6,rho12,rho9)*N3である。N3の順は既存selected B3と完全に同じ。
そのinverseはAdj3(N3)/Delta407とknown rho inverseであり、p/R/h22/Delta5を反転しない。
新しい領域は厳密に FULLU∩D(Delta407*Saux*T)。全補集合を保持する。

全source/target permutation、known column units、K11 inverse、Schur左右演算と各逆を
保存し、元15×22をdiag(I11,residual4×11)にする。残るsourceはcurrent15
`[9,11,12,14]`、targetは元32 `[1,3,4,6,7,8,9,14,17,18,19]`。
追加tail2は単位pivotとして正しく消費され、他七tailは残る。元八tailの効果を無視した
切捨てではなく、full target inverseで全て復元できる。

4×11のscalar値とrank4は未計算。全FULLU57やfinite bad locusを主張しない。
source covectorsは既存current15→old25へ戻すが、raw H7 Cdelta liftを暗黙に合成したとは
せず、元sourceへの十分条件の適用はその別の依存scopeに従う。
