# 固定された R434 の一回だけの QQ 因数分解

2026-11-04。既存の pure4 determinant は `a^4*C^3*R434`。R は434項、全次数20、
各次数は(a,e,d)=(11,16,6)以下である。元sourceは一切変えない。

入力のQQ分母を揃え、signed contentを取り、lex先頭係数が正のprimitive整数多項式に
する。その一つだけをFLINTのexact multivariate QQ factorへ渡す。全因子・重複度と
constantを保存し、constant×全factor冪を逐次掛け戻して元434係数と厳密一致させる。
APIのstatusと多項式の零性を区別する。成功はQQ因子の再乗算を証明し、別の因数分解を
試す指示ではない。Sauxには触れない。

予測では非定数因子数は20以下、一因子の支持は全次数20のdense bound1771以下。
総factor支持を過大に35420とし、入力・逐次product・output-sized scratch・一般reserveを
含むaccessible支持を150000以下に制限する。factor library内部の支持数はinstrument
できないため、この不確実性を明記し、全processを60秒・2GBで制限する。予測はlibrary
内部が絶対に小さいという証明ではない。serialized資料は合計32MiB以下にする。

因数分解が予算内で成功した場合だけ、既存p21を各返却因子で割る固定有限検査をする。
全次数から不可能ならその証明を記録し、それ以外は literal商・remainderの掛戻しを
保存する。これは新しい一般gcd/factor searchではない。全pがRを割らない既存結果と、
proper factor共有の有無を混同しない。

因子の零集合は、新しい欠損成分とは認定しない。p/Delta5/h22はいずれもchart単位に
限る。既存direct P5の領域はD(Delta5*Saux)でRを除外しないが、この因数分解は旧P4
chartの補集合を正確に整理する別資料である。全FULLU cubic57、pure残余rank、既知
curve以外の補集合解決は結論しない。
