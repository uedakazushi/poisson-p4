# 自由 M0：141生成の有限加群と九つの作用

2026年11月2日。指定有限体点 F_32003 における actual S/a³ の圧縮表示を保存する。
元135の高次 full standard basis を繰り返すのではなく、既存の有限商の正常形を使用した。
characteristic zero や全 S の結論ではない。

## 1. 小さい逆変換

R-adic preparation の141横断標準単項式 B を用い、
A3=F_32003[u0,u1,u2]/(u)^3 とする。
A3^141 の標準列は、B と base 次数0/1/2の10単項式の積、計1410列。
保存された正常形写像 C の像は1303次元である。
その1284列は相異なる単位列であり、残り19行だけの可逆 minor を解けば右逆を作れる。
base 二次の846列だけの場合も、737単位列を除いた rank24 の小さい minor で十分。

この右逆を用いて、S/a² の22本の一次関係に二次 tail を付け直す。
補正なしの一次関係をそのまま S/a³ の関係とみなしてはいけない。
二次85関係のうち一次倍60次元を除く25本を加え、計47本の A3 関係を得る。
これらと base 一次倍の F_32003 行空間は rank107である。
1410−107=1303 と C の像次元が一致するので、得た関係は C の全 kernel である。

## 2. 作用と有限代数の整合

既存の保存 full basis に対し、B の normal 八変数倍1128列と、t 倍141列を query した。
各正常形を上の右逆で A3^141 の列に戻し、normal 八変数と t の九作用を構成する。
生成時の native query、19/24 minor と全作用検査は116.27秒だった。

正式 verifier は保存 query の入力単項式と出力、C による全作用列の一致を確認し、

* 47関係×九作用の安定423件、
* 九作用の全36組×141生成元の commutator5076件、
* t² の全141列の零、
* B の単項式評価が選んだ cyclic 生成元に戻ること

を検査する。従って明示的な141生成・47関係の A3 加群上に、
正しい可換作用と t²=0 を持つ有限代数が構成されている。
この代数を actual S/a³ と同定する依存は、既存の元135に由来する full basis と正常形である。
今回はその巨大な基底を再生成せず、source SHA256、再構成手順、query を保存する。
原135式に対する別実装の代数同定を新たに完了した、という主張ではない。

## 3. 次の有限系

B の distinct monomial border は537個だが、最小 monomial border は89個だけである。
base 次数3の tail は、一つの base 三次単項式ごとに89×141=12549係数へ絞れる。
base 三次単項式は10個で、最高 base 次数の線形化には同じ transverse 作用素が現れる。

次に必要なのは、この89個の monic normal 関係を元135へ持ち上げる tracked membership。
それを使い、変形された border relations、元式の残余、および critical relations を
141生成上で処理する。既存 A3 代数の抽象変形だけでは元 S/a⁴ は決まらないので、
この source の拘束を省略しない。
J10 による cutoff は m^(r+8)⊂a^r であるため、r=4では総次数11までの有限情報で足りる。
全11変数の巨大な第三層標準基底を再試行するのではなく、
横断八変数側の89個の source certificate と、141生成の小加群に分ける方針である。
base3 compatibility の解、全 R 加群の認識、QQ への移行はまだ得ていない。

## 4. 保存と再現

データは research/wip/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.tar.gz に圧縮保存する。
47関係、九作用の全係数、小逆行列、1269 query の入出力、探索ソースとログを含む。
74MB の元 full basis は同梱せず、前 preparation archive の再構成順序を参照する。

    .venv/bin/python code/s61_minimal_free_flag_r_adic_actions_20261102.py --output reports/S61_MINIMAL_FREE_FLAG_R_ADIC_ACTIONS_20261102.json

この verifier は有限加群と作用の恒等式を再検査するが、raw135 の展開や native basis を再実行しない。
特に t²∈a³S から全 S の t²=0、t 冪零、relative 飽和除去を推論しない。
