# 混合norm係数を旧e-monic16代数へ移す入力preflight

2026-11-04。新しいdeterminant、solve、GB、点、素数、univariate GCDは実行しない。既存norm pencilの全三つの混合係数について共通G除算を確認し、最小の商の乗法行列だけを構成する。

## 1. Gを除ける範囲

旧stripped norm hと追加norm Kにはh=G*h1、K=G*k1があり、三つのstripped mixed係数MiにもMi=G*miがある。三つのliteral積を元のQQ係数で確認する。sourceの全unit stripは独立受理済みpencil reviewへ束縛する。

actual BbadでG!=0ならh1=k1=m1=m2=m3=0が必要である。これはD(G)上だけの必要条件であり、Gをglobalに逆にしていない。G=0側は既存Q/L例外の扱いに残し、この新行列から除去したとみなさない。

## 2. 最小の混合商

三つの商の項数は110,143,168、e次数は19,21,23、a次数は9,11,12である。各e最高係数は非零QQ定数。最小のm1をそのQQ単元でmonic化し、既存h1のe-monic16式Hで割る。

基底1,e,...,e15上で全16本の

```
m1_monic*e^j=H*q_j+sum_i M_ij(a)*e^i
```

を構成・再乗算する。divisionで逆にするのはHの最高係数1だけで、未知a多項式のleading coefficientは反転しない。出力16×16行列はQQ[a]上にあり、detやsolveは計算しない。

## 3. 次の計算の正確な目標

将来この乗法の非零norm delta1(a)と必要なcofactorを得れば、D(G)上のactual旧残余ではdelta1=0が必要になる。既受理の117次delta0とのQQ[a] GCDがa単元なら、この**G非零部分**を除去できる。G0のQ/L例外は引き続き別に扱う。非定数GCDは必要a集合を狭めるだけで実在を与えない。

今回normが非零か、delta0とのGCDが単元かも未計算である。入力matrixの実項数、degree assignmentによるdet次数上界、整数分母払い後のheightと保存量を記録し、次回の予算判断へ渡す。detだけの保存量と、16cofactor・元Bézout全体の保存量は区別する。以前の34.72MB guard失敗を踏まえ、detだけのboundを次回全artifact保証へ流用しない。

## 入力計算の確定結果と提案予算

source/input計算は0.5863秒で完了した。行列は合計2624項、最大一成分18項、最大a次数17。全16除算を再乗算し、処理中の最大supportは324項だった。certificateは1243343 bytesである。

numeric分母払い後の最大整数係数は318bit。実entry次数行列のassignment上界からdetのa次数は高々152、support153、係数5190bit、detだけの保守的RAW上界は254439 bytesとなる。16cofactorだけの粗い上界は6410240 bytesだが、Ve・Qsum・monic/original Bézoutの重複を含む全bundleはこれで保証されない。

次に予算判断を求める候補は**determinant一個だけ**のexact計算と、保存a117とのQQ[a] GCDである。sourceのmonic multiplication modelからnormは共通h1=m1=0の必要条件であり、この集合論的必要coverのために巨大なsolve/Bézout bundleを出力する必要はない。なおm1とh1のgeneric互いに素性は今回示していないので、detが零ならそれを保存して停止する。他のmiや別点・素数への自動切替は提案しない。今回そのdet/GCDを実行したわけではない。
