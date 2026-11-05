# native 75→72→605 制限と fine 正値上界の独立 SOURCE 査読

対象は `NATIVE_SHARED72_INPUT_SOURCE` (`4f018950…`) と `NATIVE_SHARED72_FINE_POSITIVE` (`072c20b0…`) である。独立報告は `reports/AMBIENT_B18_NATIVE_SHARED72_SOURCE_FINE_20261105_v63_REVIEW/report.json`、SHA-256 `0384d2623fab486eb0dbd1a4df239bd04e6a7ff778f7464c454ec324c541683e`。**制限・数値正値上界の SOURCE 査読は PASS**。符号付き native 応答や605個の upper は計算していない。

同じ canonical G8 の点値は `(1/2,25/8,0,0,0,5,−125/2,0)τAτB`、点 monic は `τA²=τB²=0` である。元の Koszul 符号と free4 の恒等列から、48×128 の点 source 行列全体を独立に作り、元 MU 入力の128列すべてと一致させた。固定23補列の指定23行では `diag(−I7/2,I16)` となり、行列式は −1/128。各 free 列の指定行成分をこの対角補列で消す方法で105列の kernel 基底を独立に導き、全48成分の零と free I105 を確認した。rank search や生産者の明示比率式の import は使っていない。

元619個の源の添字の一意性を確認し、563/609を源の段階で差し引いた全617関係について、literal full48 S0 との積が零、105座標からの full128 再構成が元の源と一致することを検査した。旧部分との M0X0 相殺から native kernel を推測した証明ではない。free16成分がないことも各関係で確認した。

非 opaque 605関係には固定105基底の72列・120座標が現れる。この72列を展開すると元75列を使い、元605の raw 支持に追加される4/8/12はいずれも**点 S0 が非零**である。したがって75列の個々の補助値は known partial-N による形式的 extension であり、未提供の真の N 高位値を得たとは扱わない。実際の射影された Γ 三次＋D 非 zz 二次という意味を持つのは、固定72 kernel 列への制限後の値、およびそこからの605定数写像である。

primitive の J18、N の source/physical 順序、最初八 source 行の prune、G8・monic・完全な Γ56/D は、先行の独立数値査読 `f9581e3b…` と文字通り同じ入力へ結んだ。そのため数値検査済み68本の R 行を利用できる。全 monic 対の積の支持消滅は今回も独立に確認した。単一 monic の全 trace carry と ordinary rr 一回という規約は保持した。

生産者の逐次正値加算とは異なり、源・mask・base の順序を逆にし、経路の未約分分母の LCM を第一パスで求め、第二パスで絶対値分子を共通分母へ移して集計した。全72列の S fingerprint、計4896個の Γ/D 細分 D/U、両 channel の Γ/D/合計、および605関係の両 channel の D/U がすべて厳密に一致した。有理数の比だけの比較ではなく、分母の可除性証明を含む D/U 自体の一致である。

72列のうち6列では、絶対値を取った S0 の上界は非零だった。これらでも、独立に確認済みの full48 kernel 恒等式によって**点係数だけ**を厳密な零へ置き換えた。first、second、r²ℓ、monic の高次項を零にしていない。未提供 N2 mixed/N3 は形式的に残る。D 二次では全 S0 による最高項が消え、Γ 純 r²ℓ では同じ kernel 消去と J18 における z 因子の純平面への不寄与を用いる。個別 raw4/8/12 を零とする論法ではない。

全 factor を作る前の正値クラス上界を再構成し、fine242ビット、605写像254ビットを得た。別法 native 共通上界341/1024と、その605写像354/1063も独立に再構成した。将来の whole worker は fine 上界へ戻してから OLD/FIRST/外側/Q を組み立てる必要があり、coarse354をそのままQへ掛ける許可ではない。

演算前 whole owned 上界109,449、観測保持上界104,467、演算中 committed354・演算前チェック354ビット、実行0.277秒で、hard150,000・384/1153・32MiB・60秒・2GB内で完了した。最初の独立査読では通常 pair のループ変数が primitive 入力名を隠し、fine照合後・coarse再構成前に停止した。失敗時のコードと SOURCE prefix を `_REVIEW_ATTEMPT1` に保存し、三行の変数名だけを訂正した。生産者・primary・入力・上界は変更していない。

この結果は SOURCE と正値上界に限る。OLD605の形式消去証明 `36a453a1…`、全 FIRST/Psi/comparison carry、個々の native92外側項を維持する。実605応答、全605零、HF・長さ・元H0への結論はなく、whole worker/custody と ROOT の新しい実行許可は別に必要である。
