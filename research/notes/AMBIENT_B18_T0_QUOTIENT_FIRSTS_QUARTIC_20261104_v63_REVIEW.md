# 正則商の一次微分・D0一次関係の独立監査

対象は `AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63` と、先行する `FIRSTJET_SENSITIVITY_PREFLIGHT` の凍結証明書である。著者workerを実行せず、新しいsource solve・jet・点・prime・GBを作らない。前段で独立受理されたCA/CB source frame、regular quotient membership、C0 free9表示を使用する。

## 保存source4の独立解釈

独自の疎二変数有理数多項式環で、次数4を超える項だけを捨てる。単位逆元には有限geometric seriesを使い、著者の15成分配列・convolution tableを使わない。実際の元source DAGを、保存source3、および保存source4へそれぞれ代入する。source3のframeと変数・式順が既存cacheと同一であること、追加が二族それぞれの285個の四次係数だけであることを固定する。全57式×15係数を両族で再評価し、保存残差と補正も元のJ57を用いる有限恒等式として照合する。新しい方程式解法は実行しない。

各族のgと四つのBbarの全15係数を元出力から再評価する。正則商の存在と原点値0は前段の全次数のsource-image membershipで既に証明済みなので、`Bbar_4=g_3 H_1` から商の一次微分が厳密に決まる。四次五係数すべてでこの恒等式を確認し、endpoint二係数だけの一致では済ませない。元のregular-function nodeや原点値は変更しない。

CAのH64微分は `(3/8,3/4)`、CBのG64微分は `(1,-1/2)`、残る六商は一次まで零である。しかし二族のr,sを同一視しない。保存source2の共通physical九変数への微分行列から、full9 adapterの81恒等式と非零detを確認する。これによりCBの微分をCAの全九座標へ戻すと、全四つの差 `G_CB,i-H_CA,i` が一次まで厳密に零になる。

## 全640感度の独立復元

疎dual numberでC0表示の全640関係を再評価する。旧CB四商、新CA四商は、互いに異なる八つの正orderの形式germとして保持する。既知のcoordinate quotientは、既に独立受理されたsource2/prefix恒等式に従う。八つの未知微分を先に零と置く操作は行わない。

感度が残るのは正確に列604,613,622,631で、それぞれg_A basis成分の `G_CB,i-H_CA,i` だけである。独立復元した既知部分・八感度を先行preflightと照合し、実際に求めたfull9微分を代入する。全640一次columnのrankは独立消去でも38となる。先行preflightの固定span34、許容未知効果を含む拡大span40とactual rank38を区別する。

## Visible classとsocleの範囲

g_A*tauのbasis classに七zero座標を掛けた一次類は関係spanへ入り、r/sを掛けた二類は独立に残る。全関係の原点が零だから、これは `m D0/m² D0` の厳密な非零判定である。

したがって、D0のresidue fiberへ非零のg_A basis leading classを持つ任意のactual T0 elementはsocleになれない。そのようなelementを `c e_g+m F9` と書けば、r倍の一次類は `c r e_g` であり、追加項はm²へ入るため、上の非零類を相殺できない。sでも同様である。これはraw g_A*tau自体がactual kernel elementであるとの主張ではない。

一方、m D0内に隠れたfinite-support socleは排除していない。T0の全annihilator、cyclicity、元交叉環の深さは未決定のままである。自動的な高次jetへの増加は行わない。

専用reportに現在の全依存hash、独立実行の時間・メモリを保存する。
