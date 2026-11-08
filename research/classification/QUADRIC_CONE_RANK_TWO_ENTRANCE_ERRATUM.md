# Quadric-cone rank二whole排除稿：非零Cの元source入口を訂正

2026-11-08、producerによる明示訂正。旧note／manifest bytesは変更しない。
独立reviewer quartic_reviewが保存PENCIL全33 directionsとの不一致を発見した。
**旧 `SS_BASE_REDUCTION`、`NIL_BASE_REDUCTION` とそれらを使うfinite入口／whole排除稿は、
元全sourceへの橋が未成立なのでwhole theoremとして採用しない。**
full33から追加項が落ちるというnormalizer／元Jacobiの証明を、その旧稿は提示していなかった。
最終pole／valuation議論の形式が整合していても、この欠落を補わない。

## 1. 省略した元directions

全PENCIL reportのfamily indexはparticularを0、directionsを1,…,33とする。
非零C semisimpleでは次の項を旧六係数表示が落とした。

* p26：B_bcへ4a。
* p27：B_bcへ2b。
* p30：B_acへ2c。
* p32：B_abへ−4c。

非零C nilpotentでは

* p26：B_bcへ−8b。
* p27：B_acへ4b、B_bcへ−2c。
* p29：B_acへ4c。
* p30：B_abへ−2c。

これらは全Fu-only coefficientsと連動するgenuine元source方向であり、
入口で零と置いてよいものではない。
C零caseにはこれらの追加base項が無く、全六係数の `ZERO_C_ACTUAL_FACTOR` はこの誤りの影響を受けない。
そのC零枝は独立reviewで別に全source確認済みである。

## 2. 正しい全B式

A=p1,d=p2,D=p8,r=p7,s=p14,t=p20とする。
SSではxi=p26,zeta=p27,eta=p30,nu=p32、kappa=p31、w=1+2kappaで

```
B_ab=A−rb+sa−ab−4nu*c,
B_ac=d−rc+ta−w*ac+(w−2)b²+2eta*c,
B_bc=D−sc+tb−bc+4xi*a+2zeta*b.
```

nilではxi=p26,eta=p27,zeta=p29,omega=p30、kappa=p33で

```
B_ab=A−rb+sa+2kappa*ac+(1−2kappa)b²−2omega*c,
B_ac=d−rc+ta+2bc+4eta*b+4zeta*c,
B_bc=D−sc+tb+c²−8xi*b−2eta*c.
```

両型で全Yを保持する。T1は既表示のもの、T0も全三constantを残す。
fullJacB／Lie_TB又はactual affine normalizerをこの完全式へ適用し、
旧canonical式が正当化される全unit／zero分岐と、新たな例外枝を区別する必要がある。

## 3. 証拠と受理状況

対象旧whole manifestは
`inputs/quadric_cone_rank_two_irreducible_quintic_exclusion_manifest_20261108.json`、
SHA cb8d1edcc1e5163000555e4edf16c97ead7fd74f24cc49f32582ff81ad8e7f7f。
これは歴史的producer稿として不変保存し、受理済み全空定理とは扱わない。
元PENCIL／triangular／独立全33算術証明書、既global G≤24には変更が無い。
C零の独立receiptは `QUADRIC_CONE_ZERO_C_WHOLE_SEPARATE_REVIEW.md` とそのmanifestを参照。

この修復のため新operator／profile／CASは実行していない。
以下のcomplete sourceから別V2純証明を作り、再独立監査後にだけwhole結論を回復する。
それまでは非零CのE帰属／空性は未解決のまま。新component数の変動は零。
