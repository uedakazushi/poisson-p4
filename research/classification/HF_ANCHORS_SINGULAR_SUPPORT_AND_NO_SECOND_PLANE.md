# 二HF anchorにはconcurrent m三を使える第二planeがない

2026-11-08、PRIMARY E、純証明。対象は固定した二anchorだけ。
後続finite graphを避けられる既全concurrent m三の直接適用入口を先に確認する。

Aのprimitive Dは `25c*u⁴−4b(3a*u+5v²)²−20a*u³v−50u²v³`。
そのc偏微分は25u⁴、b偏微分は−4(3au+5v²)²。
SingDではu=v=0が従い、このplane上で全偏微分は零なので、reduced singular supportは唯一のS={u=v=0}。
Sに沿うgeneric normal位数は二であり、第二rank二plane／m三concurrent sourceはない。

Bでは `D=−b(3au+2v²)²−6ac*u²v+2a*v⁴−12c²u³+4cu*v³`。
b偏微分から3au+2v²=0。
u=0ではv=0となりSを得る。
u非零ではa=−2v²/(3u)、c偏微分からc=v³/(3u²)、他の全偏微分はこれらの代入で零となる。
従ってSingDのreduced supportはSと、次のirreducible cubic cone S'の和である。

```
(a,c,u,v)=(−2s*t²,t³,3s³,3s²t),  bはvertex自由座標。
```

これはP3のtwisted cubicをb方向にconeにしたsurfaceでありplaneではない。
S'∩Sはa=u=v=0のline、S'以外のsecond surfaceはない。
S'の一般点ではΠ_uv=−2u²が非零、D=0なのでPoisson rankは二である。
既cubic-cone G非優越はpointwise parentを認識していないため、この情報だけでEを閉じない。
Bにも第二planeはなく、Sはnormal位数二のままである。

従ってどちらのanchorにも、既wholeconcurrent m三plane theoremへの直接橋はない。
これは既24/25親への包含を排除したものではない。
唯一normal character−2と完全IFT graphの先験次数4／3を使う別有限経路を適用する。
