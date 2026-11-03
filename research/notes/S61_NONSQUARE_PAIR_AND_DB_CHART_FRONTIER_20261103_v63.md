# Pair / DB chart の再開可能な bounded frontier（2026-11-03）

全体の QQ pair gate は未解決。actual pair の座標は `A=x0=W=0`、B 自由であり、`x2=v,x4=z,x1=3t0,x3=-7t0/2+u`、`t0=(30PaTc*v-6P Vc*u-4aKa*z)/(60P²)` とする。元の804項から A=0、B=-10aPJz、Q1=Q2prime=0、および literal2×6係数行列を再検証した。

Pair第0行には3個/10項/次数2の QQ cofactor による厳密恒等式 `sum(row0_j*c_j)=a*P^4` がある。safe append+unit pivot の実際の出力は1×6、466項、次数12、11,300bytes、分子/分母最大33/45bits。QQ/既知unit正規化後は次数11、51bit整数係数となる。これは pair rank≥1 のみを示す。D1/D2、さらにD3/D4の有限 primitive cofactor 探索でも次の witness は得られなかった。D3/D4 は a,P のみを target とし、係数次数の許す全 monomial target（最大14/15）も調べた。これは探索不成功であり、QQ非所属や全次数での不可能性の証明ではない。

DB!=0 の別 chart は、split frame で B_r=-20aDB を用いて r を消去する。actual Q0,Q3,D の係数を a²B_r² で polynomial にし、3×10行列を保存。2個の安全な known-unit pivot 後は1×8、593項、次数12、11,820bytes、係数最大35bits。DB はこの chart でだけ逆にする。残りの scalar ideal の gcd は1だが、これだけでは localized unit ideal を意味しない。D1/D2の有限探索も不成功。

次の実際の制約到達点は別々に保存した。

- 以前の補助3×19QQmodule計算: 300.014秒、RSS427424KiB、timeout、GB完了なし。既存 frozen v56_five_rank_gate_budget300（archive SHA dee071b4821cbfc2818ebb02077881faf6d00c87c99bcc577876aef900fd9294）を参照。
- Pair1×6 scalar の4変数 QQ liftstd: 60.013秒、RSS521232KiB、timeout、GB完了なし。work/.../pair-scalar-QQ-v66。
- DB chart1×8 scalar の3変数 QQ sequential saturation: 59.964秒、RSS43336KiB、最初の sat(a) で timeout、QQ stage 完了なし。work/.../five-DBchart-scalar-QQ-v69。先行0.033秒の initializer エラーは INVALID_rational_exponent_initializer 以下に保持し、修正版は coefficient-first serialization を使う。エラー run の exit0 は成功ではない。

同じ32003だけを使った軽量 diagnostic は両方成功した。Cloudの `sat` は **idealそのものを返す** APIであり `[1]` は使用しない。Pair は a,P の順の saturation 後に ideal(1)、DB chart は FULLU と DB の saturation 後に ideal(1) となった。各0.033/0.065秒。この結果は準備した行列の有限体上の unit ideal 証明であり、QQへの昇格ではない。既知 unit のみを使用し、未知因子を逆にしない。

補助5-rankのDB0全体の厳密証明と、四線形変数/binaryDの構造縮約は別 packet に保存・reviewされる。後者によりpair gateが成立すれば五二次式のrank5も従うため、DB!=0 scalar idealを別途解決する必要はなくなる。main constructive pivot34 はこの packet で変更しない。

軽量再現は `source /workspace/cloud-setup/activate.sh; python code/s61_nonsquare_pair_gate_frontier_20261103_v63.py`。重い失敗計算を再実行せず、actual804→pair座標、全12係数、10項QQ unit identity、expanded scalar row の全byteを約4.6秒で検証する。GB再試行、追加prime、無期限cofactor次数増加を次の既定動作とはしない。自然なcore座標による新しい構造的縮約を継続中である。
