# HとD12の非除算を調べる固定軸四次商環：実行前source

2026-11-05。準備のみ。行列式・商環値は実行しない。
新しいbeta点やphysical sourceを使わず、既存のe=d=0という係数制限を用いる。
元M12の144軸成分は独立M12_AXIS_ONE_COEFFICIENT peerで全てsource一致を確認済み。
determinantの独立確認が一係数だけだったことと、入力144セルの全source確認を区別する。

同peerの `own-integer-C12-normalization.json` は12×12の整数多項式C12と定数LCDを保存する。
全726係数、最大a次数8、既存整数係数最大110 bit。元native入力との全144辞書の一致を
準備で確認するが、旧67原子・214演算や572 sourceを再評価しない。
正確な正規化は

```
D12(a,0,0)=a^5*det(C12(a))/LCD,
LCD mod32003=30563 !=0.
```

全67literal atomのgeneric e,d指数は非負で、分母はa,Bのみの単項式。
同sourceでB(a,0,0)=-40なので、制限はQQ[a,a^-1]への合法な環写像である。
これはd=0をFULLUの点と呼ぶことではない。
新しいkappa3・anchorや未知polynomial分母を加えない。

94項Hからこの軸に残る係数はそのまま

```
h_axis=-1025280+35712000a-1148800000a²+26260000000a³-96000000000a⁴.
```

既存prime32003で先頭係数は7157、定数係数は30819で、ともに非零。
従って次数4を保持したmonic像hを使い、S=Fp32003[a]/(h)は有限free次元4。
aはSの単元である。Sを体やreduced環と仮定しない。

実行時に144成分をSへ減らし、12×12 determinantを4096 subset DPで求める。
row kの状態は選択列maskのdeterminantであり、最後の行によるLaplace符号を全て保存する。
全24576 ring積、各ring積は高々16+12=28 scalar積なので688128 scalar積以下。
入力726係数の減少・a冪・LC/LCD逆のcheck・最終a5/LCD正規化を含め750000 scalar積未満。
全field演算一時値30 bit。既存ROMの係数110 bit、LCD519 bitとそのmod変換は別scopeである。
全4096×4 DP、144×4成分、powersとguard/control/failure state込みowned28000未満、
RAW1MiB、45秒soft/60秒hard/2GB。JSONはstreamで書き、DP全体の二重copyを作らない。

余りが非零ならHがD12を六因子環で割らないことが従う。
実際、Hは完全反転積adFBCJと互いに素なので、そこでの除算関係があれば
a,Bの分母を払ったprimitive整数分子にもHが割り込む。
Gaussによりquotientは32003で整であり、この合法な軸制限とLC保持還元でも割り切れるはず。
ところがaとLCDはSの単元なので、計算した非零余りと矛盾する。
余り零ならこの固定軸検査はinconclusiveであり、H|D12は主張しない。

成功し独立peerを得た後には、既受理H-squarefreeとunit-coprimeを使って
V(H)∩FULLU∩D(D12)が非空と結論できる。これは選択二tail十分条件のgapであり、
実bad Fittingやrank低下ではない。今回の準備には値実行の許可は含めない。
