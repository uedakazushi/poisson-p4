# generic W枝に沿うCA・CB・Nの全完成環

2026-11-03。全critical idealと実odd二次jetの受理済み証拠だけを再利用する。
新しい高次source計算やGBは行わない。本corollary自体の独立監査は別に記録する。

## 全W曲線と九次minorの最初の非零係数

CA・CBの全monic critical表示の十source変数を

    (s0,...,s9)=(0,0,W,0,0,0,0,0,0,0)

と置く。七monic式、四core式、全150 Schur critical式の全てがQQ[W]で零となる。
両者のphysical tensor曲線は一致する。
実quadratic-base R29の座標でも、この曲線の微分は受理済みT_Wである。
以下で近似直線だけをactual critical曲線と取り違えない。

原点で独立な七actual odd functionをf_jとし、その非零firstjet minorをBとする。
この七function自体をregular R-coordinate s7に取る。
残る二actual odd generatorからf_jの定数係数結合を引いたh_46,h_95は、Rの全firstjetが零である。
共通kernel T=g19+WUV上の受理済み二次項は

    CA: [[4/25,−4/25],[9/25,1/25]] * [WU,WV]^T,
    CB: [[7/45,−1/5],[−1/15,−1/5]] * [WU,WV]^T.

二つの2×2 determinantはそれぞれ8/125、−2/45である。
七functionはactual critical W曲線上で恒等的に零となる。
この曲線の接線はT_Wなので、h_iのT_U,T_V微分はW*C+O(W²)となる。
h_iの全firstjetが零であるため、曲線の高次accelerationはこの先頭係数に寄与しない。
七regular constraintを保つためのdirection補正もO(W)であり、h_i側のO(W)微分との積はO(W²)にしか寄与しない。

従って七coordinateを正規化した九actual odd minorのstaircaseは

    [[I7+O(W), O(W)],
     [O(W),    W*C+O(W²)]]

であり、determinantの最初の項は **det(C)*W²** で非零である。
未計算の高次係数を零に置いてはいない。
raw七方向を使えばdet(B)はCAが−8、CBが1/128なので、先頭項はそれぞれ−64/125*W²、−1/2880*W²である。
保存certificateは七actual function、正規化九方向、全二次residual、staircaseとこの四係数を含む。

## 全成分の完成像へ適用する

全proper P0 fiber長2から得るfinite・fiber長2以下の開集合はW曲線のW=0を含む。
従ってgeneric Wでは全proper fiberが有限で長2以下となる。
既知の全critical idealはWをinvertするとU=V=0の滑らかな20次元枝となる。
そのsource-to-target微分はcorank1だから、generic Wの全fiberは一台長2である。

R29上のfinite quadratic source B=R[t]/(t²−nu)とactual image A=R⊕I*tを使う。
critical quotientに対応するJ=(I,nu)は正則九生成idealである。
上の非零minorはactual Iがgeneric WでJ/mJを全て張ることを示す。
NakayamaによりI=Jとなり、非零c∈Iからt=(ct)/c∈Frac(A)を得た後にBをnormalizationと同定する。
よってCA・CBはgeneric Wでともに

    C[[g20,s8,t²,s_i*t,t³]] ⊂ C[[g20,s8,t]]

となる。全81 kernelは九生成正則列のKoszul関係とpair productsから得る。
coneではdim29/edim38/depth21/e2、projectiveではg19となりdim28/edim37/depth20/e2である。

ここでformal minorの非零からalgebraic generic Wへ移すのは、未計算係数の数値外挿ではない。
finite source、quadratic trace、odd moduleとそのFitting条件は、必要ならétale局所化した有限代数上のものとして定義できる。
actual odd rankとregular critical quotientという開条件がformal W枝で非空であるため、元の代数W曲線のgeneric pointでも成立する。
以後の全完成型はこの二つの開条件とNakayamaから得る。
任意の特定のW≠0がその開集合に入るとは主張しない。

## Nは受理済み全idealを局所化する

Nについては、既に証明された全odd ideal

    I_N=(nu,s1,...,s6,WV,V²)

を使う。WがunitならI_N=(nu,s1,...,s6,V)は正則八生成idealである。
Nの元の九source変数でも(t0,...,t8)=(0,W,0,...,0)を全critical式へ代入して零を確認し、
そのphysical tensor曲線がCA・CBの全W曲線と一致することも保存する。
従ってgeneric WのN completionは

    C[[g20,s7,t²,s_i*t,t³]] ⊂ C[[g20,s7,t]].

coneではdim28/edim36/depth21/e2、projectiveではdim27/edim35/depth20/e2となる。
八生成正則列に対する28 Koszul関係と36 pair products、計64関係が全kernelである。
三つのnormalizationはいずれもt²,t³∈Aによるelementary subintegral extensionである。

これはgeneric W枝だけの比較である。
B18原点でのN completionは別の既知型であり、CA・CB原点の全odd idealは依然未解決である。
全W、global bad locus全体、mixed cubic、actual intersection Dの完成型はここで分類しない。
