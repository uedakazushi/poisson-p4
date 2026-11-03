# J=0: 保存行列の単位列正規化と一次 direct lift の厳密な障害

## 保存した最良状態

前の `S61_NONSQUARE_J0_DIXON_UNITS_20261103_v63` で得た39単位消去後の18×52行列に対し、**指定済み単位の共通列因子だけ** を除いた。環は引き続き

\[
 A=\mathbf Q[a,e,1/(aF B_s C(3e+2))]
\]

である。J=0 の代入、F,Bs,C の定義、57行基底への接続は直前の packet の通りで、未知の因子やJは逆転しない。

実際の one-witness API は各行を known-unit primitive にするが、固定引数の経路では column clean を行わない。残っていた共通列因子は `a` が51個、`3e+2` が33個であった。全 entry について割り算を逆乗算して検査した。これは列ごとの単位対角変換なので cokernel を保存し、新しい単位 pivot を増やした操作ではない。

その後の非零 QQ scalar による行・列対角変換は、以前の60秒 cap を超える87.1秒で正常に完了した。新しい120秒/2.2GBの bounded budget の範囲内であり、全 scalar 操作を逆変換して確認した。最良の保存状態は

`work/nonsquare-J0-pid-v63/column-units-v1/primitive.txt`

で、18×52、876,541,698 bytes、SHA256 は

`039ee56e033a9d576aae8e271483f34e10dd111e87f391eb740251dc90c9489a`

である。原 raw 状態、成功した各段階の hash、列因子の指数表、QQ scalar、native log と receipt は保持した。公開 report には再生成可能な大行列を重複保存しない。

## 一次 direct lift を排除する小行列式証明

計算量を抑えるため、単位消去の際に他の行を変えない特別な係数を探す方針を調べた。これは、固定した20×54整数係数行列 M に対して

\[
 Mq=g(a,e)e_r,\quad
 q_j=\alpha_j+\beta_j a+\gamma_j e\in\mathbf Q[a,e],
\]

を解くことである。\(g\) には次数や単位因子の上限を課さない。しかし、この **一次係数 ansatz については、すべての20行 r で非零解が存在しない**。

証明は整数係数の有限次元線形代数だけで与えられる。10個の固定整数点 \((a_s,e_s)\) で、r 以外の19行の零条件を評価する。この190×162整数行列の列は \((\alpha,\beta,\gamma)\) に対応する。各 r について、162×162小行列式が素数32003を法として非零であることを確認した。従ってその **整数小行列式自体が非零** であり、QQ 上の評価条件は単射である。恒等式を満たす q はこれらの必要条件も満たすので、q=0 しかない。

native producer と異なる Python 実装で、固定20×54行列の全係数から全32,400評価行列 entry を再構築した。さらに別の NumPy modular elimination により、各 r の実際の162行と非零行列式を保存・再検算した。これは mod p の cokernel 零性から QQ 上の零性を推論する議論ではない。固定した **整数線形写像の非零小行列式** を使っている。

この否定結果は、固定した20×54 frame 上の total degree≤1 の **多項式** q に限る。高次の係数、既知単位の分母を持つ係数、別の row frame に対する結論ではない。前の通常の row-unit witness は他の行に非零の像を持つので矛盾しない。

## J=0 上の generic full rank

同じ整数係数20×54行列を \((a,e)=(1,2)\) で評価すると、最初の20列の20×20小行列式が32003を法として非零である。従って QQ 上のこの点でも非零である。ここで \(a,F,B_s,C,3e+2\) はすべて非零なので、それまでの QQ および A 上の可逆変換を戻せる。最初の37 pivot と合わせて57行基底の full cubic rank57 を得る。特に J=0 上の **generic rank57** は厳密に確認された。

これは欠損集合が proper closed subset であることを与えるが、**欠損集合が有限であることはまだ与えない**。未知の曲線上で rank が落ちる可能性を除くには、依然として QQ(a)[e] と QQ(e)[a] の両方での unit/cokernel 証明、あるいは同等の determinantal support 証明が必要である。

## 否定的な基底改善試験

全多項式係数の Gram 行列を用いた整数行 LLL も一度だけ試した。\(\det U=\pm1\)、全係数 Gram の合同式、全多項式行の変換後 Gram を厳密に検査したが、最大係数 bit 数は1308から1296、bit 数総和は424,895,290から421,206,788への改善に留まった。続く QQ primitive 操作にも新しい行 content は無かった。この枝は最良18×52状態には使用せず、改善が小さいので反復しない。

## 再現

`polynomial_known_unit_columns_20261103_v63.cpp` は `input output known-units trace` を引数とし、直前 packet の raw18×52から列単位を除く。続いて `polynomial_QQ_scalar_primitive_20261103_v63.cpp input output factors 1` が最良状態を作る。source と trace の hash は report に記録した。

一次障害の全小行列式は、元の大行列を再生成せず軽量に確認できる。

```sh
python code/s61_nonsquare_J0_direct_lift_replay_20261103_v63.py \
  reports/S61_NONSQUARE_J0_PRECONDITIONING_20261103_v63/certificate.json.gz
```

元の20×54行列が存在する場合、native `s61_nonsquare_J0_direct_lift_20261103_v63.cpp` と `s61_nonsquare_J0_preconditioning_20261103_v63.py --output <新規directory>` は、その係数からの独立な評価再構築を含む検算を行う。

次の未解決問題は18×52 cokernel の支持の有限性である。今回の一次 direct lift の不可能性と係数成長を踏まえ、高次の sparse coefficient identity や、両座標関数体における座標 annihilator を狙う必要がある。quartic coverage、ambient nilpotence、一般 affine exhaustion の完了を主張しない。
