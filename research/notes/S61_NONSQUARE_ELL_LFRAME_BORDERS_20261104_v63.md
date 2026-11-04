# 非平方 cubic11：新 L frame の structured reduction

凍結済み L-frame preflight の10x26二次倍行317項と genuine16x31の
50085項だけを読む。一回限り60秒・2GB・32MiBの計算であり、旧 dense19
elimination、追加点・標数・GB・四次行列の計算は行わない。

まず二次倍行の2x6 pair block を消し、8x20 block を消す。
genuine16行をその10行で reduce してから、既存と同じ genuine 行
`[1,2,8,10,11]` の heavy5 を消す。ここで残る11x16を完成すれば別に保存する。
残る genuine11 のうち元 relative 行6を最後に置き、他の10行を消して
その六つの残余係数だけを検査する。初めの25行の generic 独立性は既存の
rank10/非零34minor と可逆な座標変更・九つの monomial quotient による。
新しい点で pivot を探すことはない。

core 多項式の gcd を消去前に使い、不要な共通因子を掛けない。
各行の全係数に共通な多項式因子は厳密な division を確認して除く。
これらはすべて QQ(a,e,d) 上の可逆演算で、未知因子を全域単元とはしない。
全行演算の u,v と正規化 g,c を transcript に保存する。
一つの checkpoint が完成するごとに行列・pivot・残る列を保存する。
中断した transcript の末尾は未完成 step である。

全六残余が零なら canonical11 の generic rank は10、非零なら11。
途中の resource/term guard はどちらの結論も与えない。全域57、full239、
新しい全域 kernel 方向の数には直接結論しない。

実行：

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_nonsquare_ell_Lframe_borders_20261104_v63.py
```
