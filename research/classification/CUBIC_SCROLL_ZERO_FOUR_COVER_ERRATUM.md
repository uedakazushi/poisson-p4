# Smooth scroll 旧四源 whole-cover の撤回：g11係数誤りと全補集合

2026-11-08、PRIMARY E。別作者 Fraction review が発見した producer の数学的誤りを明示する。
旧 pure note、四源評価 report、manifest の frozen bytes は変更しない。
**旧四源が全十二 Poisson source を覆うという結論、及びそれだけによる whole empty の結論は撤回する。**
独立受理や global claim の登録は行われておらず、既存成分数も変更しない。

旧 `CUBIC_SCROLL_ZERO_FOUR_WHOLE_SOURCES.md` §3 で、
元 g11 の `−a*u+b*u+4f*u` は `(δ+4f)u` であるのに、
誤って `(δ+2f)u` と展開した。
その結果 `g11=−2R4+2R2−R3` と記したが、正しくは
`g11=−2R4+4R2−R3`。
g10+g11は恒等零であり、fu=0 は得られない。
後続の f/u 分岐からの四 Segre cone被覆はこの誤りに依存するため無効である。

具体的な全一parameter反例は

`(r0,…,r11)=(2,2,0,Q,0,1,Q,0,1,0,2,0)`。

f=u=1なので旧四sourceのいずれにも属さないが、元二十四式は全零である。
別作者は元 full350／actualDも独立に確認する。
この例から rank四、既約actualD又は新componentを主張しない。

旧四個別 map の元24／全350 zero と全126 actualPf factor identityは、
それぞれの全 map 上では正しい（別作者も個別に検証中）。
失われたのは全 source を尽くす cover であり、個別 family の算術値を変更していない。
旧 report の `whole_smooth_scroll_rankzero_irreducible_quintic_branch_empty=true` は、
四個別 sourceに限定した producer booleanであり、全十二sourceの証明として使ってはならない。

修正された全被約支持は、旧 L6 と **3×3 rank≤1 Segre cone五**の和である。
正しい八元式・元24との両方向・全 outer-product map は
`CUBIC_SCROLL_ZERO_TWO_WHOLE_SOURCES_V2.md` に別 identityで保存する。
V2の全actualPf／全係数検証と独立reviewが終わるまで、旧四factor結果を whole E closureへ昇格しない。
