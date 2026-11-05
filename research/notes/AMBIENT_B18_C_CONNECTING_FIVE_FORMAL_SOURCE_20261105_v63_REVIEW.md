# 元五 head の形式 connecting SOURCE：独立検証

SOURCE 検証は PASS した。値は未計算である。
旧 Told の最小五 head、旧56列の行22,23,29,30,55から、真正 C[ell] への
形式接続写像を、実 source 参照・符号・型・正則除算・境界降下まで独立に確認した。
報告 SHA256 は
`23028007c65f9f17db4637e9af563b2e481ee8b380cb06c8abdbb951bdb1ff6c`。

producer は import せず、元 chain の Q48 と負の U8 の三非定数成分から
FULL V の全168参照を復元した。ell=0 の既存 V56 の参照と符号が全て一致する。
必要な元 scalar DAG は269 nodeであり、未使用の残りの chain を展開しない。
旧単位分解では mbar の選択24列の行列式が
−397267578125/50096498540544 で、nbar は選択行へのその単位逆の埋込みである。
型付き参照を検査し、既受理の mbar nbar=I24、d1 mbar=Vbar を用いた。

従って sigma=(I−nbar mbar)E は、受理済み自然写像の literal cycle と同じである。
sigma の最終 ELL0 係数を QQ[[z7,r]] から ell に依存しない形式級数として持ち上げる。
同じ original105 の FULL clone と ELL0 clone は別物であり、最後の座標はそれぞれ
s=r+ell/5 と s=r である。両者の Jacobian の定数行列式は非零である。
FULL m,n への置換は行わない。

Vbar sigma=0 なので、FULL V とこの持上げの積の15成分は元の形式正則環で
ell の倍数である。そこで先に RegDivell を行う。有限 Artin 商に入ってから
ell を割る操作ではない。得た eta は ell[eta]=0 を満たす真正 C[ell] の元である。
独立検証済みの D 線形 q7、ker E7 の自由行1,2,5,6、既存 Jordan 基底変換を
この順に適用する17個の constructor の型も確認した。

境界降下には、元 FULL chain の実恒等式640本 UP=QB と、完全な128本の
normalizer syzygy を用いる。mod ell の境界を FULL の元 X と係数で持ち上げると、
ell 非依存の代表との差は ell t である。RegDivell 後の差は V t であり C で零となる。
これにより全768関係、および Schur/Euler 後の825関係について降下する。
homology の単純な tensor や higher Tor の省略には依存しない。

原 FIVE_FORMAL の「旧第八生成元は f63」とする説明は広過ぎた。
独立検証では後続の不変な PLANE_ORDER 補足による限定訂正を採用する。
実 U9 の honest C0 基底の末尾は gA=node505938 で、scalar9 もその参照である。
OLD640 の Q0 multiplier の末尾は f63=node7790 である。
この二つは全七座標の実 telescope f63=gA+Σ e_i zAi で結ばれるが、同じ関数とは扱わない。
元の168参照、sigma、RegDiv の operator は訂正前後で変わらない。

この接続写像は Rbar 線形だが、途中の T7 は D のみの商である。
T7 の任意ベクトルに単純な R0 作用を置いてはいない。
正しい真正 target 構造を通じて値は既存の M7=Told/(z7,r²)Told を経由し、
三条件を満たす十七個の未知 head 係数へ帰着する。それらは依然未計算である。

再現用 checker:

```sh
.venv/bin/python code/ambient_B18_C_connecting_five_plane_source_review_20261105_v63.py
```

同じ実行で PLANE_ORDER の独立検証も保存した。記録は0.56秒、所有数上界53,136、
RAW入力13,379,320 bytesで、60秒・2 GiB・150,000 slots・RAW32 MiB以内である。
凍結出力先への再実行を拒否する。再検証は既存成果物を保存した別作業コピーで行う。
元 H0、K torsion、nil 三乗、接続写像の実数値はこの SOURCE 検証からは得ていない。
