# 線形42×42行列による全未知残余のdegree42方程式

2026-11-03。元134式の普遍三平面に沿う疎な構造から、
**一つの線形42×42行列の行列式R42**で全未知残余をスキームとして正確に定義する。
I_res=I+(R42)=I+(e)で、eは凍結済み大域polarization idempotentである。
未知残余の空性や一般長さ5を結論しない。
未知qを反転せず、一般係数はbetaだけの有限逆行列recipeで指定する。

## 1. 129個の二次式を線形normal行列へ因数化

一般betaで正則な普遍三平面の枠を使い、affine座標を(z1,z2,z3,n1,…,n42)とする。
既知五点はn=0にあり、その断面の実idealは既存五本の二次式で生成される。
元134式の平面制限の行空間は正確に5次元なので、五行を選び、
残る129行からそのspanを引くと、純zの項を全て消した式G_iを得る。
五行のrestriction minorを反転するだけで一般beta上でもこの操作を定義できる。

各G_iは次数2以下で、全単項式に少なくとも一つのnormal n_jを含む。
単項式中で最初のnormalを因子に選ぶ決定的な規則によって

    G_i(q)=sum_j A_ij(q)*n_j

と書く。**Aは129×42、全成分はqに関して次数1以下**である。
producerは中心の全129行をQQ Fractionで実際に因数化・再乗算し、
元134式と保存された五行span係数に対して正確に一致することを確認した。
一般betaでは普遍三平面枠と五行restriction逆行列で同じ有限recipeを使う。
普遍三平面の存在と5次元restrictionは先行の全16beta証明を入力とする。

## 2. 一つの共通42行とadjugate

129行のうち、中心の既知五点全てで同じ42行が可逆となった。
保存selectorで選んだ正方行列をA_square、R42=det A_squareとする。
五点での行列式はmod32003で

    1680,23338,12913,13967,26226

と全て非零。従って一般beta上でR42の既知五点代数への像は単元である。
そのnormを反転する条件はbetaだけに依存する。
matrix entriesが一次なのでdeg_q R42≤42。
保存された最高次方向の行列式が31743≠0 mod32003であり、
一般betaでこのrecipeは**正確に次数42**を持つ。

A_square*nは元134式の組合せなので、adjugateにより

    R42*n_j ∈ I  (j=1,…,42)

という実ideal membershipを得る。
これは局所的なnormal Jacobianの階数だけの主張ではなく、全q多項式の正確な因数化である。

## 3. 全未知残余とのideal等号

全アフィン代数A=K[q]/Iは凍結済みpolarization splitterからE_known5×A_resへ分かれる。
normal ideal N=(n1,…,n42)Aの商は、三平面との実交差として既知E_known5そのものなので

    N=(1−e)A。

従ってR42*N=0からR42=eR42を得る。
一方R42は既知因子E_known5上で単元だから、A内で(R42)=(e)である。
つまり正確に

    **I_res=I+(R42)=I+(e)**。

これはradicalの等号ではなくスキームとしてのideal等号である。
未知の孤立点・正次元・非被約構造は全て保持する。
cyclic determinant零やquintic H零の未知補集合を別に切り落としたものではない。
全追加局所化は普遍枠・五行restriction・既知normに関するbetaだけの非空開である。

## 4. 保存と次の有限課題

producer `code/s61_affine_residual42_20261103_v63.py`、
selector `code/data/s61_affine_residual42_20261103_v63.json`、
receipt `reports/S61_AFFINE_RESIDUAL42_20261103_v63/`。
中心QQの線形129×42係数全体を圧縮JSONで保存した。
一般betaの42次多項式を巨大展開したとの主張ではない。

尽くしの具体的目標は全134式とR42からunitを証明することであり、未完了である。
さらに強い十分条件として、線形長方形A(q)の多項式left inverseを構成できれば、
n=0が元I内で従い、既知五点が全fiberを尽くす。
そのleft inverseの存在・不存在もまだ証明していない。
この有限線形module問題は、以前の一般関数体GBを同条件で繰り返すものではない。
