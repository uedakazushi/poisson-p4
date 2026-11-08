# Smooth cubic scroll 四源 cover の独立棄却と、有効な部分算術

2026-11-08、PRIMARY E。独立 reviewer exceptional_coverage。
**全四源による被覆は FAIL。四表示源それぞれの算術は PASS。**
producer の不変 bytes は変更せず、typed な訂正と追加全源の確認を別義務とする。

## 全被覆を棄却した正確な理由

`CUBIC_SCROLL_ZERO_FOUR_WHOLE_SOURCES.md` §3 の g11 の恒等式は

```
g11=-2R4+4R2-R3
```

であり、表示された `-2R4+2R2-R3` ではない。
同じ restriction 上で g10+g11 は恒等零。
従って元24式から独立に回収できるのは R4-2R2 であって、
R4 と R2=fu を別々に零と置くことは出来ない。
全24式はこの restriction で正確に八式
`R1,R3,R4-2R2,R5,R6,R7,R8,R9` と同値。
この両向き定数線形同値を独立 Fraction で検査した。

具体的には、全 q∈C 上の一 parameter 族

```
(r0,...,r11)=(2,2,0,q,0,1,q,0,1,0,2,0)
```

は元全24式を満たす。元30×12 complete source から実 tensor を直接作ると、
全350 projective Jacobi 式が厳密に零、actual signed Euler-Pfaffian は非零五次。
全150係数と全126係数の一致を、保存 source rows とも別に確認した。
f=u=1 なので、f=u=0 を要する L/S、u=0 を要する F、f=0 を要する U の
どの閉源にも入らない。これは rank二だけに現れる無関係な補集合ではない。

## 有効な四源の算術

独立 checker は標準ライブラリ Fraction と既独立 reviewer の疎多項式基本演算だけを使い、
producer・SymPy・CAS を import しない。
四つの全 outer-product map を別に記述し、元30×12 source から全 tensor を再構成した。
各 map の全 parameter 上で、元24式、全350式、全126 actualD、
表示された全因子積と全 homogeneous coordinate degree は全て一致する。
したがって四表示源の既約五次不存在は有効だが、全 smooth scroll 枝へ広げられない。

全2640確認、0.0566秒、RSS14432KiB。追加rank・minor・operator・点列は無し。
最初の誤った g11 恒等式を実際に棄却した checker/job の exact bytes は
`inputs/cubic_scroll_four_review_first_failure_20261108.json` に保持した。
最終 negative receipt は
`reports/CLASSIFICATION_CUBIC_SCROLL_FOUR_SOURCES_SEPARATE_REVIEW_20261108.json`。

次は正しい八式の全 support を覆う追加 closed source を別に証明する必要がある。
三行三列 rank一 source への還元案は有望だが、本 receipt の受理済み結論に混ぜない。
