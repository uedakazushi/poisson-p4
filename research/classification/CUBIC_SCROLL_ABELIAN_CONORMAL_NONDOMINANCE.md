# smooth cubic scrollの全rank零・既約quintic枝：abelian conormal源による次元二十三上界

2026-11-07、主義務G。
**結論。** smooth cubic scroll S⊂P4上で零となり、実退化因子が既約被約quinticである
全rank四Poisson tensorのprojective PGL5飽和像はdimension二十三以下である。
従って未知25〜27次元G優越成分の一般元にならない。
clean条件、generic ordinary double、normalization又は導手のsmoothnessを仮定しない。
この上界を固定sourceの全E帰属又は新成分性と取り違えない。

## 1. 新しい線形必要条件の幾何学的理由

`RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md` により、nonplane S⊂SingDではgeneric multiplicityは二。
rank零ならそのfirst-normal Lie algebraはFrobeniusで、conormal idealはabelianでなければならない。
従ってSのsmooth一般点で `{I_S,I_S}⊂I_S²`。
これは以前の三十係数sourceに無かった、新しい**線形**必要条件である。
old全sourceJacobi rowspaceの負結果、native timeout又は平方membershipの負結果を再解釈してはいない。

全sourceは既証明のEuler/Hilbert–Burchから `Π=Σ_(g=0)^2 qg Ag`、三十係数。
homogeneous q_iをlocal projective ideal generator q_i/ell²へ直したbracketの補正は、
Πの全係数がI_Sに入ることからI_S²に入る。
従って `H_ij=dqi*Π*dqj` のsymbolic-square条件が正確な必要条件である。

H_ijがS一般点でI_S²に入れば、その全五ambient first partialはS一般点で零となる。
これらはpolynomialなので、integral S全体で零となる。
逆向きはcharacteristic零のsymbolic squareとdifferential squareの一致でも確認出来るが、
今回の上界には上の**必要方向だけ**で十分である。
generator bracketから全ideal bracketへの拡張も、rank零によるLeibniz補正I²で保証される。

## 2. 全線形systemとその完全kernel

q0=x0x3−x1x2、q1=x0x4−x1x3、q2=x2x4−x3²とし、
dominant affine-cone parametrization

`x=(ss*u,tt*u,ss²*v,ss*tt*v,tt²*v)`

を使う。全三generator pair、全五first partialを調べる。
Πの一係数 qg*(∂i∧∂j) に対し、Sへのpartialのpullbackは正確に

`(∂k qg)|S * ((∂i qa)(∂j qb)−(∂j qa)(∂i qb))|S`

である。qg|S=0なので、他の微分項は消える。
`code/classification_cubic_scroll_abelian_conormal.py` は全三十係数と全pullback monomialを保持し、
159×30のQQ行列を構成した。
rank十八、kernel十二を得て、全matrix、全kernel、free座標十二、非零rank minorを
`CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json` のcases[0]へ保存した。
これはchosen parameter pointではなく、必要線形条件の全kernelである。

全kernelへ元の全350 Jacobi及び全126実D係数を制限した。
Jacobiのparameter二次row rankは二十四。
実Dの単なるrowspace包含はなおfalseであり、rank四不存在は主張しない。

## 3. 十二係数の全Poisson源のexact dimension上界

同じrunnerの追加dimension modeは、上記全350のexact independent二十四行を使う。
old線形reportはその当時のproducer bytesを埋込み、追加modeに合わせて書き換えていない。
`CLASSIFICATION_CUBIC_SCROLL_ABELIAN_DIMENSION_20261107.json` のcases[0]には

* 元の350行における二十四generator index、
* 六十二のlifted polynomialと、元のgeneratorからの三百四十二のQQ ideal-membership lift、
* 各lifted polynomialのgrevlex leading monomial、
* leading supportのminimum vertex-cover六

を保存した。producer自身は全ideal liftをQQで展開して確認した。
native数学計算と確認は約0.424秒で完了した。

initial idealはこの六十二leading monomialを含むので、
全Jacobi源のdimensionは `12−6=6` 以下である。
この十分な上界には、native output全体のstandard-basis正しさ又はradical性を仮定する必要が無い。
別作者が有限QQ arithmeticだけでmembershipとleading-coverを再検証出来る。

## 4. 全PGL像への橋

smooth cubic scrollのpreserving groupはdimension六。
base PGL2三次元、relative fiber scaling一次元、u→u+v(αs+βt)二次元が忠実に作用する。
全abelian-conormal条件はintrinsicなのでこの群に保たれる。
全必要Poisson源Xのdimension≤6より、
`PGL5×X→P H⁰(Λ²T)` のfiberにはpreserving group六次元とscalar一次元が入る。
従って全projective像は

`24+6−6−1=23`

以下となり、主張が従う。
実Dの非零一点patch又はirreducible-D open条件を追加しなくても、この上界で十分である。

## 5. twisted cubic coneは別の補集合

同runnerはconeについても全78×30線形system rank十二/kernel十八を保存した。
こちらは `Π=∂c∧Σqi Vi`（Vi任意constant base vector）という十二次元rank二linear subspaceを含む。
従って全源dimension≤9という目標は不可能であり、whole-source native出力の未認証parseを停止した。
その停止143と理由は `CLASSIFICATION_CUBIC_SCROLL_ABELIAN_DIMENSION_CONE_RUNTIME_20261107.json` に保存した。
smooth scrollの上記exact証明書は変更していない。
coneのG排除は、既証明のdense point orbitに対する単一actual D(p0)≠0局所化又は
base六/vertex十二の別構造議論を必要とし、このnoteの結論に含めない。
