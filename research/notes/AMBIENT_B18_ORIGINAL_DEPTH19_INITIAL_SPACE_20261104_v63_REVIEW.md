# 元 B18 深さ19：全二次初期関係の尽くしと深さ推論の独立監査

2026-11-04。これは pure proof と固定証明書の独立監査であり、実 source の
係数輸送を再実行しない。その actual source premise は別の独立レビューを
必須 gate とする。対象は元の実 scheme 交叉 D=CA∩CB の B18 局所環であり、
正規化比較環 E4 の深さや全 Poisson Jacobi scheme の深さではない。

## 1. 監査する完全な actual presentation

S=k[[zA0,...,zA6,rA]]、最大 ideal m、F=S^5 とする。受理済み全環 Schur
operator は実有限長加群 T を

    S^825 --R--> F -> T ->0

と表す。111個の単位対を除いた後の全関係であり、部分的な initial matrix
ではない。825列は cycle657 と literal d3列168からなる。全列は mF に入る。
d3列は m^2F に入る。cycle 列の一次行列 L は rank37、従って

    dim_k(T/mT)=5,
    dim_k(mT/m^2T)=5*8-37=3。

有限長性は同じ actual T と conductor積二乗の有限底に関する既受理の証明を
使う。新しい radical 比較 operator の別の homology と T を取り違えない。

固定620個の cycle 定数結合は L の核の**基底**である。各結合には固有の
非選択657列座標の係数1があり、他の座標は37本の選択列だけ。各結合を
保存全column factorizationへ代入すると零である。非選択座標は全620個を
重複なく尽くす。従って独立性と生成性は、rank37とこの triangular形だけで
従い、未検査の kernel dimension を推測する必要はない。

## 2. 任意の形式係数を許した二次初期関係の完全性

R_j=L_j+Q_j+O(m^3)、a_j=a_j0+a_j1+O(m^2) と展開する。ここで L_j は
一次、Q_j は二次、a_j0 は定数、a_j1 は任意の八変数一次式。d3列では L_j=0。
任意の実関係 q=R*a が m^2F に入る条件は

    sum_j a_j0*L_j=0。

このとき二次項は厳密に

    q_2=sum_j a_j0*Q_j + sum_j a_j1*L_j。

二次以上の係数 a_j は関係列自体の次数が1以上なので三次からしか寄与しない。
有限表示の image は complete Noetherian module の閉 submodule であるため、
completionから隠れた追加関係も発生しない。

定数 a0 の cycle 部分は前節の620 kernel basisの結合、d3部分は任意定数。
一方、一次 a1 の寄与は一次関係と八変数の積で尽くされる。したがって
**全ての** degree2 initial relations は、次の三群の線形包に厳密に等しい。

1. rank37 の一次関係と八変数の積。
2. 全168 d3列の二次部。
3. 全620 zero-linear cycle結合の実二次部。

これは個々の列の initialを取るだけの不十分な議論ではない。定数係数による
一次 cancellation を全て含み、その後の任意の一次係数も含めている。
unit Schur消去、Euler lift、zero-origin補正で生じる二次変更は、source gate
の証明に従い既知 d3 二次関係と一次関係の変数倍を法として扱う。
未射影 X2 の等式を仮定しない。

## 3. 実第三双対は全二次初期関係を消す

二次係数空間は k^5 tensor Sym^2(k^8)、次元180。保存 monomial順の179番は
第5 target の rA^2 である。固定第三双対 lambda は文字どおり

    lambda=coefficient of e_4*rA^2,
    lambda(e_4*rA^2)=1。

固定証明書により、上の第1・第2群の全464列を lambda は消す。
actual source/cache輸送の別の独立監査が、第3群の全620列についても
lambda=0を証明する。その輸送では未知の正則商二次係数を保持したまま感度が
消え、第三 native functional も空であり、既知部分の全620値が零である。
未知二次係数を零と置いて得た値ではない。

従って lambda は前節の**全**二次 initial relationを消し、自身は非零。
もし e_4*rA^2 が T で m^3T に入れば、Fでの差は m^2F に入る実関係となり、
その二次項に lambda を適用して1=0を得る。よって

    gr^2_m T !=0,
    length_k(T/m^3T) >= 5+3+1=9,
    L:=length_k T >=9。

全二次商の次元が1、全長が9、m^3T=0とは主張しない。

## 4. 元 actual D の深さへ戻す

既受理の同じ actual 比較 K=ker(D0->E4)、C=coker(D0->E4) は DVR
Delta=k[[ell]] 上有限であり、**rank_Delta K=4**。これは rank C=4と異なる
源 theorem で、criticalCA1+Pa1+CBcritical0+Qb2を合計した証明を使う。
dim_k C[ell]=4 と正確な完全列

    0 -> K/ellK -> T -> C[ell] ->0

により

    mu_Delta K=L-4>=5 > 4=rank_Delta K。

DVRの有限加群構造定理から K の torsion は非零。これは torsion の長さ・指数
を求めた結論ではなく、少なくとも一つの torsion generator が必要という結論。

actual E4 は H0=0なので、image(D0->E4)も H0を持たず、H0(D0)=H0(K)。
conductor支持上で ell は Noether parameter だから、これは上の非零DVRtorsion。
全正次数成分は受理済み T1=T2=T3=0 により深さ1の安定 actual tail に一致する。
従って全相対 graded D の H0 は零次数 H0(D0)に等しい。元の共通 regular19
群方向を戻すと、既受理の上下界 {19,20} のうち深さ19が成立し、scalar chartで
一方向を除いて projective深さ18となる。

結論は元 B18 の **actual CA/CB intersection** の cone depth19/projective18。
Tの全長、H0長、torsion指数、D0の具体的socle代表、nil三乗の零・非零、全分類、
一般beta尽くしは決めていない。数値source監査とこの初期空間・深さ証明を
分離して記録し、どちらもPASSして初めて結果を採用する。
