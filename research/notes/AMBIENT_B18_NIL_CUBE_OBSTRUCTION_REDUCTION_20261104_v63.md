# 元 nil 三乗の二段 obstruction と honest Tor₂ source gate

2026-11-04。対象は共通 regular G19 を除いた元 B18 の零次数環 D0。D=k[[ell]]、ell=5(sA−rA) とし、rA-adic 深い T-tail と混同しない。N=Nil(D0)、K=ker(D0→E4)、P=N/K である。既存 transfer512、product44、upper7、actual DVR rank11 の独立監査を依存として束縛する。新しい source/jet/積係数を計算せず、三乗零・三乗非零の判定をしない。

## 1. 真の二段 obstruction

既受理 N³⊂tors_D K と NK²=0 を使う。P は D-free7 である。乗法

    β : K ⊗D Sym²D P → N³,   (k,p,q) ↦ k lift(p) lift(q)

は canonical である。lift を K の元だけ変えると、差は NK² に入って零となる。H=imβ=KN² と置く。三次積を H で割ると

    γ : Sym³D P → N³/H

も canonical になる。lift の差を含む全項は少なくとも一つの K と二つの N を持ち、H に入る。従って

    N³=0  ⇔  β=0 かつ γ=0。

ここでは D0-algebra の分裂も R0-module の分裂も使わない。D-free P の基底を選べば七つの lift が存在するが、その source を現在の rank 証拠から構成したとはしない。

μ_D K=L−4≤13 なので、β の image は13·28=364個以下の D-generator product、γ の image は84個以下で張られる。合計448は積の generator 数の上界であり、N³ の長さの上界でも、これらの積の実計算でもない。β=0 だけでも N⁵=N²N³⊂N²K=0。γ=0 だけなら N³=KN² より N⁴=KN³=K²N²=0。両 obstruction は元指数の上界を改善する正当な別 gate である。

W N²=0 を用いれば β は K/(K∩W) を通り、二つの P と γ は P/image(W→P) を通る。これは W=N、W=H0、あるいは W⊂K を仮定しない。**K/JK は β の正当な入力ではない。** 既知 JK²=0 が正当化するのは K² の pairing だけであり、JKN²=0 は未証明である。

## 2. ambient source の J-adic order と strictness の差

J=IA+IB=(zA0,…,zA6,rA²,rA sA,sA²)。全 actual cycle z=(a36,c8) に対して b0=Qeven a∈J、B=Σβi ci∈J、νA∈J である。product44 の全項から、source vector が J^p と J^q に入れば積は J^(p+q+1)R44 に入る。cube44 は J²R44 に入る。

これは F^qN=image(kerV∩J^qR44) に関する filtration であり、F^qN=J^qN とはしていない。従って cube の ambient order から J N³=0 や cubezero を結論できない。N² modulo JN の隠れた部分は Tor₂(C′,R0/J) の edge image、N³ modulo J²N の隠れた部分は Tor₂(C′,R0/J²) の edge image に入る。この記述も image の消滅や honest lift を計算したものではない。

## 3. Tor₂ connecting map の exact source 指定

A=R0/J の free resolution は Koszul(z7) と (r,s)² の Hilbert–Burch resolution の tensor である。最初の rank は (1,10,44,112)。本 packet は d1,d2,d3 全 sparse polynomial 係数を保存し、d1d2=d2d3=0 を全係数で確認する。R0/J² の31基底も保存するが、31-jetからの Artin–Rees strictness は一切主張しない。

元 normalizer 比較の **ell 切断前** V: R0^56→R0³、H1=K、H0=C を使う。source は accepted chain の literal Q48 の非 even 三成分と −U8 であり、ell=0 の V56 cache を代用しない。Tor₂(C,A) の cycle を c∈R0³⊗P2 に持ち上げ、

    (I3⊗d2)c = (V⊗I10)h,   h∈R0^56⊗P1

を **全 formal source 上で**満たす h を取る。そのとき

    δ₂[c] = [(I56⊗d1)h] ∈ K/JK。

constraint は30×692、output は56×560である。h の選択差は kerV⊗P1 なので output の差は JK に入る。c の V-boundary 変更では h を d2b だけ変え、d1d2=0で output は同じ。d3-boundary は d2d3=0で同じである。これが欠けている connecting image の正確な source recipe である。radical 比較 V: R44→R³ についても同じ構成は30×572となり、N/JN の edge map を与える。

この packet は constraint kernel の basis や actual(c,h) lift を求めていない。mod-J、mod-J² の homology だけで全 lift を回復したとはしない。Tor₂→K/JK を完成させても、それだけで節1の β が零にはならない。

## 4. 最小の次 source 要求

448-product route には、actual P の七基底 lift と K の高々十三 D-generator の **全次数 source identity** が必要である。現在の rank11、非標準 D-module 分解、Tnil の長さ・生成数、有限長 H0 の存在はこれらを代用しない。もう一つの route は product44/cube44 の全 cycle image を actual780-boundaryへ入れる exact symbolic factorization である。

新 source/jet を自動で追加せず、まずこの不足した lift/connecting map から、既存 cache と重なる部分だけを選ぶ。元 nil 三乗は有限原点障害のまま未決であり、全 nilideal の生成元、torsion exponent、正確な全長、全 nil 指数は主張しない。
