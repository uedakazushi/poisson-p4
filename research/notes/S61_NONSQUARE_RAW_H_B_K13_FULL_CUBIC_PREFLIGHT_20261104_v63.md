# 非空Ka13上のfull cubic source入力preflight

2026-11-04。actual Bbad=Ka13であってもfull cubic57のrank欠損とは限らない。ここでは保存済みのPOWER-free 8単元checkpointからactual17×24入力のsource・型・support・高さだけを読む。新しい係数評価、点、素数、Gaussian、minor、GBは実行しない。

入力はCUBIC_ZW25_UNITS/checkpoint.jsonのremaining rows/columnsで、408成分、354非零、numerator84253項・1113monomial、次数(13,23,6)である。全分母は非零整数cと1・a・Bの積である。既存8単元minorは−a^4*F*B²*C²/4608で、全FULLU上で単元である。今回のKa13にはFULLUの明示逆元があるので、8段階と各分母はKa13全体で合法である。POWER207 cofactorと未知parameter逆元は閉包に含まれない。

三つのsource parameter a,t,sに関する前の312個の認証済みnodeを利用する。新たな1113monomialは、既知のcoupled monomialまたは冪node高々3個の積として表せる。有限boxの整数indexだけの動的計画法でこのrecipeを作り、Ka13係数は評価しない。用いる既知nodeは272個である。各値の分母のLCMと係数の絶対値上界、g13のcompanion l1 normから、将来の408成分の整数高さ・raw JSON量を保守的に上から押さえる。

将来の係数入力では全cのLCMをcLCMとして、M_ij=p_ij/(c*u)、u=1,a,Bに対してC_ij=(cLCM/c)*(aB/u)*p_ijをKa13上で直接計算すればよい。C=cLCM*aB*Mは単元scalar倍なのでrankを変えない。aB/uをparameter多項式に展開して84253項を増やさず、低次の既知Ka13値として掛ける。全monomial値やcofactorを重複serializeせず、行ごとの408値とsource indexを保存する。

この17×24のrank17は、既受理のprojected criterion・product41・Cdelta unit lineからfull cubic57を含意する。rankが小さくてもfull57の否定にはならない。特にTheta曲線上でだけ証明されたheavy rank6や10×43 exact criterionを、このKa13へ持ち込んではならない。

係数を保存できた後の別手法として、Ka13-linear写像Ka13^24→Ka13^17をQQ上の221×312行列に直す案がある。分母とg13のleading coefficientのp単元性を別途確認し、許可された素数で221の非零係数minorを得れば、QQ上で全射、従ってKa13全体で全射を証明できる。この論理はKa13が体・既約であることを仮定せず、13個の成分を一度に保持する。今回その素数計算も行わない。

再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_full_cubic_preflight_20261104_v63.py --output work/K13-cubic-preflight-fresh`。実際の係数生成には別の一回予算が必要である。

preflightは0.433秒。既知node分母LCMは1826bit、係数の絶対値は2^61以下、g13のcompanion l1 normは2^5以下だった。高々三因子のrecipeにより、numerator評価の共通分母は9399bit以下・分子9951bit以下。cLCM*aBで分母を払った408値でも、分母13013bit以下・分子13804bit以下、raw JSONの保守的上界23905536 bytesである。matrix出力5304個とmonomial cache14469個のQQスロットは150000項guard内にあるが、全cacheを重複serializeしてはならない。

既存P80_CUBIC_GENERIC/source17x24-input.txtは曲線代入前の**普遍的(a,e,d)入力**を保存している。このファイルの全84253係数・408分母型を元checkpointから組み立て直したASCIIとbyte単位で照合した。次段はこの短いflat入力を読めばよく、以前のGaussianの全historyを複製する必要がない。P80の曲線式や特殊化値を借用するものではない。
