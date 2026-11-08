# unimodular零最高部とfactor最高部：全E入口の復帰

2026-11-08、PRIMARY `E_multiple_hyperplane`。producer純証明、独立/root査読前。新算術無し。
全α=J+β+γ、係数次数0/1/2、dα=0、Pfα非零constantを使う。
volume dualでP=C+L+Qへ戻すと全C6/L24/Q60を保持した元unimodular sourceである。

## 1. γ零の全sourceはCA

γ=0では定数Pfの正次数条件がJ∧β=β²=0を与える。
β(V)はJ-primitiveな五次元二形式空間の完全等方線形空間。
wedge quadratic formは非退化で、そのWitt indexは二なのでdimβ(V)≤2。

零又は二次元以下のdecomposable二形式の線形空間で、任意二つのwedgeが零なら、
非零constant一形式aを共通因子として持つ。
実際、二つの独立decomposable形式の対応する二planeは一次元で交わり、
その共通一次元がa。dim一の場合も任意の一つのfactorを選べる。
従ってβ=a∧b(x)、bはlinear一形式と書ける。

i_vJ=aとなる非零constant vector vを選ぶ。a(v)=0、
J∧a∧b=0はb(v)=0と同値なのでi_vβ=0。
従ってi_vα=aがconstantであり、逆Poisson tensorのlinear Hamiltonian場はconstant。
`UNIMODULAR_LINEAR_HAMILTONIAN_WHOLE_MEMBERSHIP.md` §2の全CA族を適用して、
**γ零の全C/L sourceは閉CAに入る**。
β零も同じ議論で含まれる。一次評価rank不足だけをmembershipの根拠にしていない。

## 2. すべてのfactor最高部はlinear-Casimir又はpullback

γ非零で係数が非constant gcdを持つとする。
もしgcdにlinear factorが無ければ、γ=qΩ、qは既約quadratic、Ωはconstant非零二形式。
γ²=0からΩ=a∧bはdecomposableで、dγ=0はdq∧a∧b=0を与える。
qは二つのlinear coordinatesだけのbinary quadraticとなり、C上でlinear factorsに分裂して矛盾する。
従ってγ=ℓη、ℓはnonzero linear、ηは係数linearの二形式である。

γ²=0を整域上でℓ²で割るとη(x)²=0。偏極により
η(V)⊂Λ²V*は任意二つのwedgeが零となる線形空間である。
decomposable二形式の線形空間の二rulingの初等分類により、次のどちらかになる：

1. あるnonzero constant aに対しη(V)⊂a∧V*。
2. ある三次元constant一形式空間Wに対しη(V)⊂Λ²W。

この分類は、独立な二元をa∧b,a∧cへ取り直し、第三元とのwedge零を比較して、
第三元もaで割れるか、全てがΛ²span(a,b,c)に入るかを判定すれば得られる。
低いimage-rank、重複factor、ηの零係数も同じ二つのsupersourceに入る。

第一の場合γ∧a=0で、volume dualなQはaをlinear Casimirに持つ。
従って新しい全linear-Hamiltonian membership定理を**全低次延長**に適用出来る。

第二の場合、Wを消すnonzero constant vector vに対してi_vγ=0。
γは閉だからCartan公式でLie_vγ=0も成立する。
従ってγは三変数からの実pullbackである。
係数が第四変数に依存する未知fiberを仮定で除いたのではなく、閉性から全依存が消える。

これでfactor最高部の全sourceは、既知閉和へ入るlinear-Casimir枝と、
**full pullback最高部の全低次延長**という正確な未処理枝へ帰着した。
後者をGのβ非優越だけでE完了とは呼ばない。

## 3. 六foliation閉包と今後のwhole-source義務

既受理の全六最高部閉包の入口を変更しない。
R13全閉包はlinear-Casimir、primitive regular R22は全T16又は空、
全退化R22のcommon-kernel部分はpullback、common-kernel無しは既空。
S23全閉包とE3のrank一nil-incidenceはpullbackへ入る。
E3 regular境界のlinear-Casimir枝にも新定理を適用出来るが、
rank二nil-incidenceから閉L112へ入る枝は、全log/pole源のE membershipを別に要する。

したがって、零最高部とfactor最高部を落とさずに、残るunimodular Eの主な全sourceは
「full三変数pullback最高部」「閉log integrating-factorの非empty残源」に整理される。
generic最高部だけの稠密性から未知低次fiberを被覆したとは主張しない。
各残sourceの全C/L、元Jacobi/Pf、全translation、non-split及び全unit零は引き続き保持する。
