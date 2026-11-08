# 新五成分catalogue：全tensor・divisor・一般orbit metadataの別作者監査

2026-11-08、別作者 `/root/quartic_sources`。対象は147/1156 cutの `research/COMPONENTS.json` 末尾五recordと `FIVE_CANDIDATES_FULL_COMPONENTS_AND_ONE_ORBIT_FRONTIER.md`（bc27300d）。新しいPoisson計算、operator、point、rank、jet、profileは無い。catalogue自体は変更していない。元immutable sourceの印刷式と独立証明へ、全十括弧及びmetadataを手で照合した。

## 1. 全十括弧とsigned actualD

全recordのcoordinate順を固定し、未印刷括弧は零として全十を比較した。五sourceは順に、Ga9 linear-w canonical四倍、nodeの全R6、cusp-flexの全R6、cusp接線・割線の全R6、constant DOUBLE forced canonicalそのものである。

|record末尾の識別子|signed Euler-Pfaffian D|factor/multiplicity|本質変数|
|---|---|---|---:|
|GA9_LINEAR_W_E6_CUSP_TANGENT_23|−4h(4h³w−8h²y²+y⁴)|(1,4)/(1,1)|3|
|NODE_FLEX_TANGENT_THROUGH_FLEX_23|−4hx(hy²−x(x−h)²)|(1,1,3)/(1,1,1)|3|
|CUSP_FLEX_SECANT_THROUGH_FLEX_23|−6hx(hy²−(x+h)³)|(1,1,3)/(1,1,1)|3|
|CUSP_CUSP_TANGENT_CUSP_SECANT_23|−9hx(h²y−(h−x)³)|(1,1,3)/(1,1,1)|3|
|CONSTANT_MODULAR_DOUBLE_LINE_CUBIC_23|h²(x²y+xh²+h³)|(1,3)/(2,1)|3|

Pxy=0である三plane源のPfaffianは `-Ux*Vy+Vx*Uy=-det(U,V)` なので、第一二も−4/−6である。初回catalogueの+4/+6はdivisorとして非零scalar同値で成分証明に影響しなかったが、signed式との整合をrootへ通知し、rootが現在の−4/−6へ訂正した。旧producer/projection bytesは変更していない。

Ga9の四次はwについてprimitive一次で既約、その唯一のplane特異点は(h,y)=(0,0)のE6 cusp、hは接線。node cubicは[h:x:y]=[1:1:0]にordinary nodeを持ち、hは変曲接線、xはその変曲点を通る通常接線。初回日本語nameの曖昧さもrootが現在のこの表記へ修正した。第二のcubicはcuspで、hは変曲接線、xは変曲点を通るsecantで、両線ともcuspを避ける。第三はhがcusp接線、xがcusp secantである。DOUBLEのcubicはyについてprimitive一次、唯一のplane cuspはA2でhはcusp secant。

全Dは非零なのでgeneric rank四。Ga9と三plane源はreduced、DOUBLEだけnonreducedというmetadataは正しい。三plane源でR6はPfに寄与しないので、全magnetic係数を保持した式でこの結論が成立する。

## 2. 次元・一般orbit・商

五つとも全source／軌道の飽和閉像は23次元である。これは合成証明bc27300dと各元whole-source theoremの入力を用い、新たに成分性の局所計算をしたとはしない。

Ga9とDOUBLEは、それぞれ元全sourceが一つのactual23軌道の閉包へ戻る全normalizerを持つ。保存Euler24／actual51 rank48からorbit23が確認され、商次元零のmetadataと整合する。

第一二については新独立全源lemma `TWO_LINES_M3_WHOLE_TWO_POISSON_FIELDS_AND_GENERIC_ORBIT.md`（4d2817d9）の全R6場を使う。

```
node:
 Z=(8x-4h)∂z-8y∂w, X=-2(x-h)∂z+2y∂w, Z+4X=4h∂z;
cusp:
 Z=(8x+2h)∂z+8y∂w, X=-2(x+h)∂z-2y∂w, Z+4X=-6h∂z.
```

strict h-Casimir liftのBVとclosed-log Hamiltonianから二場は真Poisson field。全R6で独立modEuler、commuting square-zero matricesなのでGa²を与え、全orbit≤22。保存M47のnonzero minorの開下界と合成して一般orbit22、閉像23との差から商一が正しい。

第三は `CUSP_TANGENT_SECANT_ALL_LAMBDA_ORBITS_AND_GA2_STABILIZER.md`（31af0b4f）を元全式と照合した。全finite lambdaのstabilizerは厳密に

```
(z,w) -> (z+a*h,w+3a*x+d*h)
```

であり、Ga²、orbit22。相異なるlambdaは同値にならず、全A6閉像23のgeneric商一というmetadataが一致する。I零境界へ一般orbit一定を外挿していない。

## 3. 主張していない値と全成分性の射程

全五の `generic_transverse_algebra=null` は適切である。第一二／第三の固定点D/Mから一般scheme厚みを推定していない。DOUBLEの真一normal二次余核零も、bc27300dが明記する通り成分性の前提には使っていない。smoothness、local unique、nil指数の主張は追加されていない。

`meets_SNC_log_symplectic_locus=false` も適切である。各閉像はactualDのessential≤3というclosed determinantal条件内にあり、非零SNC五独立超平面のessential五点を含めない。

全成分性と相違は、既知142＋六閉候補の全支持有限被覆、全prior142陰性、候補相互非包含により、包含するambientcomponentがこの五候補自身以外へ入れないというglobal reduced-support論証である。各source式・次元・orbit metadataはこの論証と整合し、局所obstructionの未完を隠していない。

監査時点のcoverage147–148／残TRIPLE一候補はbc27300dという固定cutでは正確である。後続TRIPLE→Ga9constantの陽性実族が別作者受理された後だけ、最終exact147／global exhaustionへ更新出来る。その更新は本監査に先取りせず、integratorが行う。

**監査結論：rootによるsignedD二箇所とnode nameの修正後、新五recordのtensor/D/rank/dimension/orbit/商metadataに追加defectは見つからない。** exact record projectionを別inputへ保存し、current catalogueのcoverageやproof-navigationの後続更新と区別する。
