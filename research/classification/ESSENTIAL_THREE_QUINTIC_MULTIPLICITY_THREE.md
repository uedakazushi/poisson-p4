# free平面五次・最大尖点重複度三の全曲線を三候補へ還元する

2026-11-07、PRIMARY E。全本質三変数既約五次Poisson定理から、
対象curve Fは有理cuspidal quintic、全μ=τ=12、free exponents(2,2)である。
ここではこの全curveのうち最大尖点重複度三のものを、三つの具体parametrizationへ還元する。
Nambaのnormal form一覧や未確認の一次文献の一意性を仮定しない。
候補の全十八係数Saito核・成分性・既知帰属は、それぞれ別に確認する。
first104用focused runnerと四reportは凍結したまま、新計算は行わない。

## 1. 重複度三の尖点はE6又はE8

重複度三の尖点をpとする。接線の交点重複度は三より大きく、Bézoutから五以下なので四又は五である。
三と四、三と五の最大公約数は1であり、Puiseux pairはそれぞれ(3,4)、(3,5)で尽きる。
従ってpの型はE6（δ=3）又はE8（δ=4）である。
算術種数六を使い切らないので他の尖点qがある。
別の重複度三尖点があれば、その二点を結ぶ直線との交点重複度は少なくとも六となり不可能である。
従って他の尖点は全て重複度二である。

## 2. 二つの尖点を用いる全degree5 linear series

正規化をP1_(s:t)、pの逆像をs=0、qの逆像をt=0とする。
target座標はp=[0:0:1]、q=[1:0:0]に取る。
pとqを結ぶ直線は接線にはなれず、交点重複度3+2=5でBézoutを使い切る。
従ってそのsectionをY=s³t²とできる。

Xはpで少なくともorder3、qで非零なので

    X=s³(a s²+b st+c t²),   a≠0.

XからYの倍を引き、target scaleを取るとX=s⁵+A s⁴tとなる。
Zはqで少なくともorder2、pで非零なので、同様にYの倍を引いて
Z=t⁵+B st⁴+C s²t³とできる。
従って全曲線は

    φ_(A,B,C)(s:t)=(s⁵+A s⁴t : s³t² : t⁵+B st⁴+C s²t³).       (1)

に入る。ここまでparameterのgeneric仮定は置いていない。
三sectionは同時零を持たず線形独立である。
line bundleのdegreeは素数五なので、非直線のimageはdegree5で、mapはbirationalである。
従ってこのsource上では算術種数六のδ budgetを使える。

pではYはorder3、XはA≠0ならorder4、A=0ならorder5である。
従ってA≠0がE6、A=0がE8に対応する。

## 3. p,q以外の尖点は高々一つ

r=s/t、Y≠0のchartで(1)は

    X/Y=r²+A r,
    Z/Y=r⁻³+B r⁻²+C r⁻¹.

と書ける。
p,qと異なる点が尖点になる必要条件は両座標の微分が零であることなので

    2r+A=0,
    3+2Br+C r²=0.                                  (2)

A=0では非零の候補rが無い。
A≠0では候補はr=−A/2だけで、高々一つである。
したがって曲線全体の尖点は高々三つである。

## 4. E8の場合：二候補のうちmonomial型はfreeでない

A=0ならpのδは四で、他にあるのはqだけである。
全δ=6、全特異点cuspidalという既済条件からqはδ=2のA4でなければならない。
qでの局所座標は

    Y=t²/(1+A t),
    Z=(C t³+B t⁴+t⁵)/(1+A t)

（s=1）なので、A4に必要なのはC=0である。
A=0,C=0のときはB=0でもB≠0でもqはA4となる。

B=0のmapは(s⁵:s³t²:t⁵)で、imageはY⁵−X³Z²=0である。
これは非trivialなdiagonal PGL3対称性を持ち、一次logarithmic derivationがEuler以外にもある。
free exponents(2,2)ではそのような一次場は無いため、このmonomial型は対象から除かれる。

B≠0ならsourceのdiagonal変更とtargetの三scaleでBを1へ移せる。
従ってE8の場合の全admissible curveは一候補

    (s⁵ : s³t² : t⁵+st⁴)                            (E8+A4)

へ入る。これはMoeのC3Bに一致し、DS Theorem4.6のk=2もこの型のfree curveを与える。
新しい全18核の再現は別義務として保持する。

## 5. E6の場合：qを選び直せばC=0で尽くされる

A≠0ならpはE6、δ_p=3である。
残りのδ budgetは三で、p以外の尖点は高々二つである。
尖点がqだけならqはA6で、qが二つならA4とA2である。
いずれの場合も、qとしてA6又はA4の点を選べる。
その選択を(1)へ用いるとC=0となる。
初めにC≠0であった場合を捨てたのでなく、もう一つのより高いδの尖点へqを付け替えた。

C=0、s=1で

    Y=t²/(1+A t),   Z=(B t⁴+t⁵)/(1+A t),
    Z−B Y²=t⁵((1+AB)+A t)/(1+A t)².                (3)

1+AB≠0ならqはA4である。
1+AB=0ならδ_q≥3で、δ_p=3と合わせてbudgetを使い切るのでqは正確にA6となり、他の特異点は無い。
従って二尖点E6+A6型はAB=−1であり、A,Bとも非零である。

AB≠−1ならδ_p+δ_q=5なので、残る尖点が一つ、δ=1でなければならない。
(2)をC=0で用いると、その存在条件はAB=3である。
AB=3の場合、二微分が零になる点はp,qと異なるのでδ少なくとも一の尖点であり、
全budgetから正確にA2となり他の特異点も無い。
従って三尖点E6+A4+A2型はAB=3である。

source rのdiagonal scaleでAを1へ移し、target scaleを併用するとABは保たれる。
ゆえに全admissible E6 curveは次の二候補に尽くされる。

    (s⁵+s⁴t : s³t² : t⁵−st⁴)                       (E6+A6),
    (s⁵+s⁴t : s³t² : t⁵+3st⁴)                      (E6+A4+A2).

## 結論と残義務

最大尖点重複度三の全free quinticは、上記三つのPGL3型のいずれかである。
monomial E8+A4の補集合も一次対称性で排除済みであり、連続parameterを残していない。
全fixed-F Poisson源は既済の全磁気inverseで対応する一つの23軌道へ入る。
各候補がfreeであることはDSの該当family/examplesで報告されるが、
今回の全18核・全cotangent/周囲・成分帰属の独立算術再現はまだ最初のmult4代表だけである。
本稿からさらに三成分を数えない。

最大尖点重複度二のA12、A8+A4、3A4、A6+3A2の候補と、
そのconfiguration・normal formの尽くしは別に残す。
FZの三尖点以上の一覧とそのdual次数誤記の補修は
`ESSENTIAL_THREE_QUINTIC_CURVE_FRONTIER.md` に記録してあり、本稿のlinear-series証明には使わない。
