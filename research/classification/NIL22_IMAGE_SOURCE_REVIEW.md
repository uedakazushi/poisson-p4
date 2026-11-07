# nil22 image 全十一源有限分類の独立査読

2026-11-07。`E_multiple_hyperplane`。
`NIL22_IMAGE_SOURCES_AND_COMPONENTS.md` 第1–5節の全支持還元を、別作者が元90係数から
独立に再構成した。**固定した nil22 image cofactor の rank4 被約支持が、実 GL2 正規化を除いて
全十一個の明示的既約源に尽くされることを受理する。**
十一源が相異なる十一成分であるとは主張しない。成分性や既知閉和への帰属は別の証明に依存する。
今回までに別途受理した二新23成分と三つの regular limit 被覆から、未被覆は六源となる。

独立実装は `code/classification_nil22_image_source_review.py`、
報告書は `reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json`、
有限契約は `research/classification/jobs/nil22_image_source_review.json`。
producer・producer helper を import せず、独立 Fraction 全線形消去と SymPy の厳密係数演算を用いた。
最終実行は606恒等式で PASS、1.486秒、最大RSS54944 KiB。
既存 producer 報告・コード及び既済成分の独立証明書は変更していない。初回604恒等式の証明書も元のpathで保持し、V2はweighted cubicのstabilizer保存を明示する二恒等式を加えた606検査の別実行記録である。初回・最終の数学的分岐結果は同じで、各実行source bytesをそのまま保存した。

## 1. 元90係数からの全19源

cofactor を `N=x∂q1+y∂q2` に固定する。
定数部が image 内であるため、最初の affine 平行移動で定数部を零にできる。
この定数部の正規化を終えた全源が今回の対象である。

全90係数の div/Lie 方程式110行に、全最高二次 Q の `Q∧N=0` を80行加えた。
零行も含めた全190×90行列を独立 Fraction 消去し、階数71、核19を得た。
元 source 報告第4 chart の全19係数が全190行を満たし、係数写像の階数も19である。
元の並び `(q1,x,q2,y)` を `(q1,q2,x,y)` へ実際に並べ替えてから比較した。
これにより全源は、追加の切片仮定なしで

    P=F∂q1∧∂q2+ΣMij(p)∂qi∧∂pj, p=(x,y),
    F=F0(p)+(y+div M2)q1+(−x−div M1)q2

となる。Ppp=0、M は p のみ、次数≤2、F0 は任意の二次以下の二変数式である。
F0 の六係数を最後まで保持する。`D=det M` の次数は3以下である。

## 2. Hessian 表示の必要十分性

M の二行を planar vector field U,V とする。直接展開により

    [U,V]−(div V)U+(div U)V=(−Dy,Dx).

全19源について、二つの非自明な実 Jacobi が正確に

    J(q1,q2,x)=−[∇D+adj(M)p]2,
    J(q1,q2,y)= [∇D+adj(M)p]1

となることを独立に検証した。残る二 Jacobi は恒等零。
従って全 Jacobi は `adj(M)p=−∇D` と同値であり、ここでは D による割算すら不要である。
原点 p=0 で ∇D=0 なので、D の一次部分は零となる。

    D=D0+D2+D3,
    S=Hess(D2), H=S+Hess(D3)/2, Hp=∇D.

すると `(adj M+H)p=0`。
C[x,y] の二変数 Koszul syzygy により各行は `(y,−x)` の次数≤1の倍数である。
2×2 adjugate は行列の成分に線形であり、この syzygy の adjugate は
p と affine-linear covector の外積になる。従って全 M は

    M=−adj H+p⊗V(p), V=V0+V1p

に尽くされる。これが全必要表示である。
行列式補題と Hp=∇D により

    det M=det H−V·∇D.

従って必要条件は `D=det H−V·∇D`。定数項比較から D0=detS。
逆に、この式から detM=D が成立し、表示そのものから `adj(M)p=−Hp=−∇D`。
よって全 Jacobi も成立する。div から定めた F に任意の F0 を加えてよい。
この必要十分性は点標本、Hessian 可逆性、形式変形又は未証明の平坦性を使わない。

## 3. 実 GL2 作用と binary cubic の全四型

同時変換 `(q',p')=(Aq,Ap)` は cofactor N を保つ実 GL4 変換である。
その全 tensor と determinant の変換則は

    M'=A M(A⁻¹p') Aᵀ,
    D'(p')=det(A)² D(A⁻¹p'),
    S'=det(A)² A⁻ᵀ S A⁻¹,
    V'(p')=A V(A⁻¹p').

F0 は `det(A)F0(A⁻¹p')` へ移る。A 可逆なので、これは二次以下の六次元空間の可逆線形写像。
磁気係数を消したり、非零条件を課すことはない。

複素 binary cubic は、三単根、二重根＋単根、三重根、零の四型で尽くされる。
前の三つはそれぞれ `xy(x−y),x²y,x³` に正規化できる。
PGL2 の根への作用に加えて、scalar 行列 tI は重み付き cubic を t 倍するので、
非零係数も必ず正規化できる。零 cubic は別枝として扱う。

S=[[u,v],[v,w]] と置く。全三次係数を V1 の四未知数の線形系として解き、
階数は三単根4、二重根3、三重根2であることを独立に検証した。
表示した解はそれぞれ残自由度0,1,2を全て持つため、特殊な V1 は漏れない。

## 4. 三単根：四 Hessian を二源へ

D3=xy(x−y) では V1=−I/3。
二次係数から V0=`(1+w/6,−1−u/6)`、次いで `u+v+w+9=0`。
その代入後、全残式は正確に

    −f1 x/6+f2 y/6,
    f1=u²+2uw+15u−6w−54,
    f2=w²+2uw−6u+15w−54.

`f1−f2=(u−w)(u+w+21)` なので、全被約支持は u=w 又は u+w=−21。
前者で `f1=3(u+6)(u−3)`、後者で `f1=−(u−3)(u+24)`。
従って Hessian は元 note の四個だけである。Hessian の逆元を一度も使わず、
singular Hessian の補集合もこの解法で空になる。

第一から第四への `A=[[1,0],[1,−1]]`、第四から第三への `A=−swap` について、
実 GL4 tensor 変換と全六磁気 F0 の可逆移送を全括弧で再現した。
よって `stable_line_conic`,`stable_three_lines` の二つの全六源を残す。

## 5. 二重根：四分岐と全零補集合

D3=x²y では V1=diag(a,−1−2a)。二次係数が

    V0=(av/2,−1−(a+1/2)u)

を一意に与え、残る全係数は次の三式と同値である。

    (4a+1)w=0,
    v[2−(a+1)u]=0,
    av²−[4+(2a+1)u]w=0.

全 residual をこれら三式の明示線形結合として独立に一致させた。
従って、零係数の補集合を含む完全な枝分けは

| 条件 | 必要な値 | 残自由変数 |
|---|---|---|
| w≠0,v≠0 | a=−1/4,u=8/3,w=−3v²/64 | v≠0 |
| w≠0,v=0 | a=−1/4,u=−8 | w≠0 |
| w=0,v≠0 | a=0,u=2 | v≠0 |
| w=v=0 | 条件なし | a,u |

となる。
実 stabilizer A=diag(r,1) は D3 を保ち、`v'=rv,w'=r²w`。
V0,V1 の全移送も独立に検証した。
v非零なら r=1/v、v零/w非零なら r²w=1 を選び、
`double_A`,`double_B`,`double_C`,`double_D` の四つへ正規化する。
複素平方根は可逆な normalizer の選択で、w=0 は別枝に残す。
double_D の a,u は全域自由である。

## 6. 三重根：三分岐と判別式零の監査

D3=x³ では

    V1=[[-1/3,0],[c,d]], V0=(m,n),
    m=−u/18−cv/3,
    (d+2/3)v+cw=0, (d+1/2)w=0,
    um+vn=3w, vm+wn=0.

これら五式が全 residual の係数と正確に一致することを検査した。

w非零では d=−1/2,c=−v/(6w)。Δ=uw−v² とすると

    m=−Δ/(18w), n=−vm/w,
    Δ²+54w³=0.

最後の式は低次二条件から正確に導かれ、w非零なので Δ非零。
従って Δ の反転が新しい未処理補集合を作ることはない。
実 stabilizer

    A=[[r²,0],[t,r]], t=rv/w, r²=−18w/Δ

により S'=diag(−18,−6)、V0'=(1,0)、V1'=diag(−1/3,−1/2)。
この全 S,V の恒等式を独立に再現した。
具体的には r非零、`w=−6/r⁴,Δ=108/r⁶` の表示で全条件を検査し、
任意の元の w非零解から上の r² を取れるので、この表示は全支持を覆う。
v=0 の点も同じ式に含まれ、v を反転しない。

w=0,v非零では m=n=0,d=−2/3,c=−u/(6v)。
r³v=1、t=ru/(2v) により S'=[[0,1],[1,0]]、V0'=0、V1'=diag(−1/3,−2/3)。
全移送を再現した。
w=v=0 では m=−u/18 と um=0 より、被約支持上 u=m=0。
n,c,d は全域自由である。
以上は `triple_A`,`triple_B`,`triple_C` の完全三源であり、
平方・立方根の零境界を落としていない。

## 7. cubic 零：Hessian の全三階数

D3=0 では D=detS+pᵀSp/2。
必要十分な全条件は

    S V0=0,  (V1p)·Sp=−pᵀSp/2.

S の変換則は `S'=adj(A)ᵀ S adj(A)` でもある。
2×2行列では adj が GL2 の全単射なので、これは複素二次形式の通常の congruence と同じ階数分類を持つ。
従って rankS=2,1,0 を全て扱える。

rank2 では S=[[0,1],[1,0]] へ正規化でき、V0=0。
V1 の全線形解は `diag(λ−1/2,−λ−1/2)`。
rank1 では S=diag(1,0)、V0=(0,n)、V1=[[-1/2,0],[c,d]]。
各全係数線形系の階数5／3を検査して解の完全性を確認した。
rank0 では D=0 であり、元の rank4 仮定に反する。この排除理由以外の一般性は使わない。
よって `quadratic_A`,`quadratic_B` の二源が残る。

## 8. 十一源と帰属の射程

canonical S,V0,V1,D3 を全十一個とも別に再入力し、
`M=−adjH+p⊗V` と全六磁気 F0 から P を構成した。
全 Jacobi、divP=N、実 Pf=−D、及び元報告の全150斉次係数を独立に一致させた。
全追加パラメータに隠れた非零条件はなく、各源は既約 affine 空間の像である。

| 源 | 全自由係数数 | 今回と別の既済帰属 |
|---|---:|---|
| stable_line_conic | 6 | 未解決 |
| stable_three_lines | 6 | 未解決 |
| double_A | 6 | 未解決 |
| double_B | 6 | 未解決 |
| double_C | 6 | 未解決 |
| double_D | 8 | 全源 CN31_25 |
| triple_A | 6 | 全源の飽和閉包が新23成分 |
| triple_B | 6 | 全源の飽和閉包が別の新23成分 |
| triple_C | 9 | 全源 third25 |
| quadratic_A | 7 | 未解決 |
| quadratic_B | 9 | 全源 half25 |

成分性の根拠は `NIL22_IMAGE_COMPONENTS_REVIEW.md`、
三つの全源被覆は `NIL22_IMAGE_LIMITS_REVIEW.md` にある。
本有限分類だけから残る六源を被覆したり、新成分と数えることはしない。
これで E の残る nil22 image 義務を、全境界を保持したこの六源に限定できる。
