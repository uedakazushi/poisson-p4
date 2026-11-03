# B18原点でのCA・CB conductor：実七制約を解いた後のdim20と重複度上界4

2026-11-03。受理済みactual odd二次jetとU³ classを再利用する。
本稿は七actual odd regular constraintの消去を明示し、monic critical source sectionをその代わりに置くことを避ける。
全odd ideal、binary cubicの混合係数、conductorのdepthやCM性はまだ決定しない。

## 1. 最小表示でのG19方向と二次式

各actual imageの有限normalizationを B=R[t]/(t²−nu)、A=R⊕I*t と書く。
Rは正則29次元である。conductor商はintrinsicに

    C=A/conductor=R/I

である。normalizationとconductorはcanonicalなのでCはGL作用で保たれる。
R自体がG-stableという未証明の仮定を置かない。

独立な七actual odd function f7をregular R-coordinateに取り、これを消去する。
S=R/(f7)は正則22次元で、C=S/L、L⊂m_S²が最小表示である。
保存した全physical GL列から、Cのorbit方向g19とWUV方向を合わせた22個が全minimal tangent kernelを張ることをQQで検査する。
これらの七方向との29×29 unitも保存する。

**二次式のG依存を除く補題。** C上の十九orbit derivationは、Sの各座標の像を任意にliftすることでSのderivationにliftでき、Lを保つ。
そのconstant tangentはそれぞれ∂/∂g_iである。
h∈Lの二次初期形をqとすると、D_i h∈L⊂m²の一次項は∂q/∂g_iだから零である。
従って全minimal quadratic relationはg19に依存しない。
受理済みWUV restrictionがその全二次式を決め、spanは正確にWU,WVである。

## 2. actual七制約を解いてもU³が残る

全150 actual odd functionの定数係数結合h63を保存する。
これは七linear classと二quadratic classを引いた、以前の受理済みcubic witnessそのものである。
そのRの全firstjetは零、WUV上の二次restrictionも零である。
前節の補題により、七actual constraintを消去した後のh63は**全22変数で三次以上**となる。

受理済みのexact monic U-sectionをgamma(U)とする。
ここでは全150 actual odd functionのgamma上の値がO(U³)であり、h63の値は

    CA: (12/125)*U³+O(U⁴),
    CB: (−4/243)*U³+O(U⁴)

である。branch derivativeも三次まで零なので、Schur値とactual odd値の区別はこの次数で既に監査済みである。

七actual f7を本当に零へ解くformal IFT補正はO(U³)から始まる。
保存certificateはその先頭29成分も示し、七functionのcubicを打ち消すことを確認する。
h63のfirstjetが零なので、このO(U³)補正がh63へ与える影響はO(U⁴)である。
従って上の非零cubic係数はactual七制約の消去後にも残る。
しかもh63は全変数で三次以上だから、gammaの二次acceleration、G-gauge、tangentを保つformal coordinate changeもそのU³係数を変えない。
これでJ/mJのclassだけでなく、actual conductor最小表示の本当のcubic initialが得られる。

## 3. genuine initial idealと上界

最小表示S22の初期idealはWU,WVを含み、さらにU³係数a≠0のcubic h3を含む。
h63のorderは三次以上なので、前節と同じlifted derivationにより∂h3/∂g_iは二次初期ideal span(WU,WV)に入る。
g次数ごとのEuler公式を使うと、h3のg依存部分は(WU,WV)の倍数として除ける。
従ってformal productを新たに仮定せず、g19に依存しない本当のcubic initialを得られる。
actual critical W曲線とV曲線は全source式で証明済みであり、そのconductor内の接線はそれぞれW、Vである。
三次以上のh63をこれらの曲線へ制限すると、W³とV³の係数は零になる。
WU,WVの倍数を引いて、cubicを

    f3(U,V)=U*(aU²+bUV+cV²),       a≠0

と取れる。b,cを計算したとは主張しない。
従ってgenuine initial idealは

    I0=(WU,WV,f3(U,V))

のg19延長を含む。これはheight2のHilbert–Burch idealである。
行列 [[V,aU²+bUV+cV²],[-U,0],[0,−W]] の三つのmaximal minorがこの三生成元となる。
I0の商はCM一次元、Hilbert seriesは(1+T)²/(1−T)、重複度4である。
そのg19延長からactual gr(C)への全射があり、dim C≤20となる。

一方、真のG19×W critical familyは、保存したactual tangent20minorが非零なので、conductor内に滑らかな20次元formal imageを持つ。
従って **dim C=20** であり、actual multiplicityは4以下である。
G19×Vも同様で、二つの20次元接空間のunion rankは21である。
従ってW、Vの二つの閉じた滑らかな20次元reduced support germは異なるminimal support branchである。
他のbranchやembedded構造をここで排除しない。

## 範囲

cone conductorはdim20、embedding dimension22で、projectiveではdim19となる。
個々のnormalizationのoffdiagonal selfpairは、conductor上のfinite normalization fiber productに含まれるので、原点近傍でその次元も20以下である。
offdiagonalな開部分自体がfiniteだとは仮定しない。
CA同士・CB同士のselfpairであり、CAとCBのintersection Dの上界ではない。

actual multiplicity4、actual conductor CM性、depth、全Iの最小生成数上界は、ここでは証明していない。
特にCAの追加selfpair証拠を、まだsecond flagのないCBへ移してはいけない。
lower multiplicityとinitial上界の一致を使う後続theoremは別に監査する。
