# 一般full PGL stabilizerを尽くす独立監査

2026-11-04。凍結 `Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63` を監査する。
元の全affine tensorと保存sigma変換に束縛し、producerを呼ばない別の疎多項式環演算で
全旗付きcaseを確認する。対象はQQ(u,v)、C(u,v)とその代数閉包上の一般stabilizer。
特殊parameterの群の跳びは扱わない。

## 全PGL変換を覆う理由

一般tensorのquinticは3H+Q。唯一の三重Hを保つPGL変換はHの補集合ではaffineである。
Q|Hの唯一のvertexであるz方向を保ち、既受理の唯一のcompatible constant bivector
line C=∂y∧∂zも保つ。affine変換はconstant tensorとSchouten compatibilityを保つので、
後者の一意性を使うことができる。この三条件から全候補は

    x'=ax+bw+tx, w'=cx+dw+tw,
    y'=ex+fy+gw+T, z'=Hx+Iy+Jz+Kw+Z

となる。Poisson scalar κも自由に残し、最初からκ=1やb=64cを仮定しない。

独立計算した元tensorのdivergenceは `(0,4u/3,ux+uvy+1/48,0)`。
affine pushforwardに対するdivergence covarianceからκf=1。
全Pxy,Pywの一次係数は、κf=1を代入するとb=64c,d=aを強制する。
したがってunorderedな二本のinvariant lineの一方のcaseを黙って捨てていない。
Pxwのconstant係数からκD=1、D=a²−64c²。従ってf=D。
divergenceのz成分のy,x,w係数を順に比較すれば
`J=D²,e=(D−a)/v,g=−64c/v`。この別の導出はa≠0またはc≠0のchartを選ばない。
Pxy,Pywのconstant係数とdivergenceのz constantからtx,tw,Iも一意に決まる。

## 二次caseと残る全係数

QQ(u,v,a,c)係数の疎多項式環で、元の全六二次bracketを独立に変換した。
Pyzのx²,xw,w²だけが残り、その三分子Q1,Q2,Q3はprimaryと一致する。
二つの多項式恒等式から、一般のu²−256と32u+2vが非零ならac=0,c=0。
D≠0とQ2=a²(1−a²)からa=±1だけになる。この消去は代数閉包上でも有効である。

各aについてH,K,Zを未知数、Tを自由なy平行移動として、元の全raw homogeneous tensorを
五変数・QQ(u,v,T,H,K,Z)係数環で変換する。全150係数のうち五つだけが未消去の線形式。
保存式と全て照合し、独自に非零3×3 minorを計算する。このrank3により解は一意。
実際の解は

    a=1:  exp(3T N/(4u)),
    a=−1: exp(3(T−4/(3v²))N/(4u)) E

であり、双方の全150係数が零となる。これは二つの正しい行列を見つけただけの
主張ではなく、intrinsicな全affine旗から始めた候補の尽くしである。

## 群の結論

両sigmaの差Rについて `Mminus*R=Mplus`、`R²=Ga(−1/(128uv))` を独立検査した。
`E=Ga(1/(256uv))*R` はE²=I、EN=NEを満たす。Gaのx,w linear blockはIでEは−Iなので
EはGaの外。Gaのparameterはy平行移動から `t=3T/(4u)` と回収できる。
二つのcomponentについて正則な逆parameterがあり、標数零のstabilizer subgroupはsmooth。
従ってactual geometric generic projective stabilizerは

    Stab_PGL5([P(u,1,v)]) ≅ Ga × μ2

であり、QQ(u,v)上にもsplitしている。必要条件は最終的にD=κ=1を強制する。
このμ2は同じtensorのisotropyであり、parameter間の同定を与えるKlein fourとは別物。
特殊parameter、boundary、全family normalizerの群構造は含まない。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_generic_stabilizer_review_20261104_v63.py`
（一回60秒/2GB上限、新しい点・prime・GB・solver探索なし）。
