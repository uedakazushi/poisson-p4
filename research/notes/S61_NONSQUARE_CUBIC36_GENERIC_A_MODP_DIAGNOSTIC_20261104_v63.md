# Safe36 coker の generic-a finite-field diagnostic

2026-11-04。source は frozen V6 `after21`、SHA `c65b3d80b2207e75ba794821bb41ebc3ffc864959dab27307af77cdeab9d8fd1`。read-only preflight と同じ literal36×70 coefficient matrix を使う。許可された一回だけの diagnostic であり、準備・sanity・std/profile を合わせた wall budget は60秒、子 process は2GB、artifact 合計32MiBを上限とする。途中終了なら input/hash/log を保存し、自動再試行しない。

最初に既存 Singular runtime の F32003(a) 算術、literal2×3 right-inverse identity、std/reduce と vdim=0 の小さい sanity を実行する。成功 marker がない場合は main GB を始めない。新しい install、環境の再設定はしない。

main ring は F32003(a)[e,d,z]。a は coefficient-field の transcendental parameter として残す。全421,667 rational coefficients を literal にmod32003へ還元する。source の denominator は全て単元であり、零になる14 numerator termsも件数を保存する。同じ(e,d) monomialの a-powers を一つの coefficient polynomial として括るだけで、行・列・係数を省略しない。

localizer は dFBC であり、F=500a+90e²+3e−24、B=(50a+18e²−9e−2)d+60e−40、C=100a+18e²−9e−2。`(z*d*F*B*C−1)*I36` を36列として付加するため column map は R^106→R^36。Jを単元にはしない。d factor は alternative ad57 row basis の比較条件から必要である。

std が完了すれば quotient unit、dim、vdim、basis column/term count を記録する。actual std と leading module は小さい場合だけ literal coefficient付きで保存する。有限標準 module monomial basis も上限内なら保存する。multiplication operators、annihilator、QQ certificate の追加計算はこの予算に含めない。

GF の unit/finite result を QQ に昇格しない。QQ の exact identity、または monic finite Z_(p)(a) model、全 localizer の no-escape と exact special-fiber gate 等が別途必要である。timeout/OOM の場合も数学的な rank/drop 結論はない。この packet は新しい成分や actual fullquartic surjectivity を主張しない。
