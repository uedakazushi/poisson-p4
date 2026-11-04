# 4×4ノルムのexact行列式・随伴・元minorからのmembership

既存Mのuniform clearing s=D*a^10による整数行列N=sMを使う。
固定24 permutationの行列式Rと、adj(N)の第1列vだけを計算し、全4式Nv=Re0を
検査する。元のモニック除算恒等式からV=Σvj*d^j、Q=Σvj*qjとして
s*g*V−s*H*Q=R。既存H=αf+βgを代入すると
R=(−sQα)f+(sV−sQβ)g。このliteralな全係数membershipを保存する。

全native恒等式を別の標準Python整数分母払い積和でも確認する。個別係数評価や
有限体からの推論は使わない。Rは非零で、QQ contentのみを外した後、gcd(R,U)を
反復してU-supported因子を全て除く。各除去でRの商とUの商を保存・再構成する。
最後のhとUのcoprimalityを別実装exact GCDで確認する。F/B/Cなどの既約性は
仮定せず、未知因子の因数分解は行わない。

hが非定数なら、B-rankbad点の(a,e)射影はh=0という必要曲線に含まれるだけで、
その曲線に実際の階数欠損点が存在するとは主張しない。hが定数ならRはFULLUで
単元なのでpair idealが単位になりBの全域rank3が従う。この場合分けは出力の
厳密なhに基づく。full239など他の未解決条件は別である。
