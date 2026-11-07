# nil31 image内 n7=0：r=−1/2 全源は新しい25次元成分

2026-11-07。主課題 `E_multiple_hyperplane`。
`MODULAR_NIL31_COMPONENT.md` §4 に残した二つの全源を扱う。
**r=−1/2 の全十変数源は、新しい25次元の全 Poisson 成分を与える。**
r=−1/3 の全源の初回計算は25次元下界・接上界27だった。続く真正二法線の根基原点証明により、別の25次元成分と判明した。後続証明は `NIL31_IMAGE_THIRD_COMPONENT.md` を参照。
以下の成分性は元の全源を直接扱い、非線形座標同値や CA/CB への単なる接触を使わない。

## 1. 全源を保つ射影座標と六括弧

旧 affine modular 場は `y∂x+z∂y`、定数成分は image 内である。
旧全15変数源の残余二枝は

| 枝 | 固定条件 | 非零 rank4 条件 |
|---|---|---|
| half | n10=−1/2,n6=n7=n9=n11=0 | p=n5≠0 |
| third | n10=−1/3,n5=n7=n9=n11=0 | q=n6≠0 |

両枝で n0,n1,n2,n3,n4,n8,n12,n13,n14 および p 又は q の全十変数を残す。
もう一つの重複超平面 z を無限遠に取る射影置換

    (h,a,b,c,d) = (old_z,old_h,old_x,old_y,old_w)

を行う。これは行列式−1の PGL 変更である。置換後の lift から

    E ∧ (0,0,−rc,−rh,0)

を引くと h-Casimir 代表が得られる。これは Euler 鉛直 gauge である。
全旧源との全係数恒等式を producer が保存したので、z=0 上の点を捨てる操作ではない。

全磁気二次式と二つの一次係数を

    A = A0 h²+A1 ah+A2 dh+A3 a²+A4 ad+A5 d²,
    K = 3r+t+1,          S = ua+vh+(t−r)d

と置く。元パラメータとの対応は

    (n0,n1,n2,n3,n4,n8,n12,n13,n14)
      = (A3,A0,A2,A5,A1,A4,u,v,t).

六括弧は正確に

    {a,b}=r ac,              {a,c}=r ha,       {a,d}=0,
    {b,c}=A+K(hb−c²/2),      {b,d}=g+cS,       {c,d}=hS,

ただし half では r=−1/2,g=p a²、third では r=−1/3,g=q ha である。
全十パラメータの strict Jacobi が恒等的に零であり、実 affine homogeneous Pfaffian は

    half:  (p/2) h a³,       third: (q/3) h²a².

従って実射影退化式はそれぞれ `(p/2)h²a³` と `(q/3)h³a²`。
指定された非零 p,q 以外の条件を付けず、全十変数で一般階数四となる。

## 2. half 全源の全空間成分性

固定点を

    (A0,A1,A2,A3,A4,A5)=(1,2,3,5,7,11),
    u=2, v=3, t=5, p=13

とし、tensorを12倍して整数係数を使う。
Dは全350射影 Jacobi 式の全150係数方向への微分である。
Vは全 Euler 鉛直方向、Gはそれに全 PGL 作用と scalar を加えたもの、
Fは G と上記**全十パラメータ**の実際の微分を合わせた行列である。
QQ 上の非零小行列式と全 `DF=0` により

    rank D=100, rank V=24, rank G=47, rank F=50.

κを追加した拡大家族の方向による下界ではない。元の全十源そのものの F50 である。
よって射影化・Euler 商の後で

    dim(actual source saturation closure) ≥ 50−24−1 =25,
    dim(full ambient tangent) =150−100−24−1=25.

源の PGL 飽和像閉包は既約である。従って上記点は全 Jacobi 空間の滑らかな25次元点で、
そこを通る全成分は一意である。既約な像閉包はその成分と同次元なので一致する。
これで half の**全源閉包**を一つの全 Poisson 成分として同定した。
全一般横断代数、有限 jet の積分又は形式積分解を追加仮定していない。
p=0 も同じ既約多項式源の閉境界としてこの成分へ入る。

## 3. 既存全成分との相違

既存92成分のうち、射影次元25は T01,T39,T40、positive-weight NEW25、
および nil31 の独立三超平面(2,2,1)成分である。他の次元の成分とは成分性と次元から異なる。
一般 rank2 の成分とは閉 rank 条件でも異なる。

本新成分の全 rank4 部分は、proper な射影掛け算写像による閉条件

    D = H³ L²

を満たす。この条件は非被約因子を持つため SNC 点を持たず、前三つの MPS25成分とは異なる。
positive-weight NEW25 の一般因子は3Hと滑らかな rank5 quadricであり、H³L² 条件を満たさない。
既知 nil31 25成分の一般因子は三つの独立超平面の2,2,1型で、三乗因子を持たない。
従って本新成分はこれら同次元の五成分とも異なる。
既に成分性を得た後の、全閉包に保たれる条件による比較である。

## 4. third 全源の初回次元診断（後続支持証明への入口）

同じ全係数表示を r=−1/3,g=qha とした全十源で、同じ数値点 q=13 を使うと

    rank D=98, rank V=24, rank G=48, rank F=50.

従って全実源閉包の射影次元は少なくとも25であり、全接空間の上界は27である。
全 kernel は52次元、実 F50 に補な真正 normal は二本ある。
この初回rank計算だけではこの枝の成分性又は包含を主張しない。後続の `NIL31_IMAGE_THIRD_COMPONENT.md` は真正二法線の全二次余核と根基恒等式で、別の25次元全成分であることを証明した。
次の課題は二法線の一回の最弱支持判定、又は全源を保つ実幾何接続である。
同じ退化因子の次数・重複度だけで half に含めてはならない。

## 5. exact 証明書と review 範囲

- Producer: `code/classification_nil31_image_zero_components.py`
- Report: `reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json`
- 有限契約: `jobs/nil31_image_zero_components.json`

全旧15源の parameter substitution、全 PGL/Euler 恒等式、普遍 strict Jacobi/Pf、
固定点、全十微分、全D/Fのhashと非零 QQ minors の行列添字・値を保存した。
実測4.535秒。これらは producer の exact 確認であり、独立作者による再実装と区別する。
統合担当は独立 full-ambient/source review 後にのみ新成分を計数する。
