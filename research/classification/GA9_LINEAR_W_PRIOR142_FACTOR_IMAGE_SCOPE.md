# Ga9 prior142 synthesis：factor multiplication の閉条件の正確な射程

2026-11-08、root手読解による用語限定。`GA9_LINEAR_W_PRIOR142_WHOLE_COMPONENT_EXCLUSIONS.md`（SHA256 `7c1079b25eaedcc27699d974937036e6a94abb4b74ba69bba69a1ea5bbb60fd6`）、その142-row projectionと17-input manifestは全て不変に保つ。

冒頭の「reduced1+4 multiplication image は閉条件」という表現は、以下の **full M14** に限定して読む。

```
M14 = image(P(Sym¹ V*) × P(Sym⁴ V*) → P(Sym⁵ V*)),
       ([L],[G4]) ↦ [L*G4].
```

射影空間の積は proper だから M14 は closed。ここで G4 の reducible／nonreduced、L と G4 の共通因子、全 factor collision を含める。「被約かつ既約quartic」という open 部分自体を closed と主張しない。

Ga9 dense orbit O の actual divisor は被約 H と既約 F4 の積なので、その閉包 Y の非零Pfaffian open 上では D∈M14、essential rank≤3 という二つの閉条件だけを用いる。Pfaffian零も含める場合は対応 affine cone で零を加える。

同dimension二十三の generic factor-type 比較は、仮の包含が C=Y を強制し dense open Ga9 orbit と accepted generic source が交わる、という等号の議論である。reducedness 又は irreducibility が boundary で保存されることは使わない。高dimension成分の除外で使う他の所定factor-degree/multiplicity像も、全 splitting／zero-edge／factor collision を含む proper multiplication image である。

Ga9 の quartic singularity は E6（normalization vanishing (4,3,0)、局所 leading h³+y⁴、Puiseux pair (3,4)）。A6ではない。旧103のA3+A2共通tacnode-tangent成分の全閉包からの非包含は、凍結済み `GA9_LINEAR_W_MARKED_PAIR_SEPARATE_REVIEW.md` の全marked DVR/P1 normalization 論証を使い、singularity adjacencyだけで除かない。historical142 projectionにはこの旧103が一行だけ含まれ、重複加算は無い。

この限定で prior142 の全行理由と結論は変わらない。新成分性・局所支持／smoothness／正確な成分数は別の whole-cover 入力とroot採録の義務である。
