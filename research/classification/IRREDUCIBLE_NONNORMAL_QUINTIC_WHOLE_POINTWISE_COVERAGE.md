# 被約既約・非正規actual quintic：全点の既知成分被覆

2026-11-08、PRIMARY E。純合成定理。新算術無し。独立source/仮定合成reviewは別receiptで保存する。

**定理。** 複素P4のPoisson tensorがrank四で、actual Dが被約既約非正規quinticなら、その点は採録済み閉成分の有限和に入る。

## 導手支持の全入口と全rank

既受理 `NONPLANE_CONDUCTOR_POINTWISE_KNOWN_COVERAGE.md` の全支持入口を使う。非正規integral quinticの導手にはintegral曲面Sがあり、全正規化/非split/特異/vertexを保った既分類により、そのSのdegree一から六の有限支持型が尽くされる。ここでSをsmooth又はsplitと追加仮定しない。

Sの一般点上でPoisson rankは零又は二である。rank四は退化因子上では起こらない。よって次の三つの全pointwise theoremで全支持が尽くされる。

1. Sがnonplaneなら `NONPLANE_CONDUCTOR_POINTWISE_KNOWN_COVERAGE.md`。全rank零/二、degree二–六の全型はT37/T38/T39/T40又は空。quadric coneはFu-only全方向を戻したV2 theoremだけを使い、撤回V1を復活させない。
2. Sがplaneで一般rank二なら `PLANE_RANKTWO_CONDUCTOR_WHOLE_POINTWISE_COVERAGE.md` と `PLANE_RANKTWO_CONDUCTOR_WHOLE_COMPOSITION_SEPARATE_REVIEW.md`。全R/全flags/全root partition、m二・三・四及びm五以上補集合を保持した既知閉和又は空。
3. Sがplaneでrank零なら `PLANE_RANKZERO_CONDUCTOR_WHOLE_POINTWISE_COVERAGE.md`。全coHiggs/global V/scalar/zero symbolを全位数で尽くし、全higher sourceが既知閉和又は空。

少なくとも一つのSを選べば十分で、導手の全既約曲面を同時にnormal formへ固定する必要は無い。同じtensorに複数Sがあっても既知和への所属結論は変わらない。

各theoremは全source/actual inverse/全zero complementのpointwise結論である。G非優越、fixed支持のparameter上界又は一点のorbitからこの結論を推測しない。有限和の閉性と全支持入口を合わせることで、全被約既約非正規D枝が被覆される。

この定理で新成分を加えない。被約既約normal Dの既証明とは別に扱い、可約及び非被約D全支持は本定理の射程外として残す。従って全E・総成分上界・分類主定理はまだ未証明である。
