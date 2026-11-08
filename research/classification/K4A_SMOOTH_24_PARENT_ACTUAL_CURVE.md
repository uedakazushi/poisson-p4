# K4a の全一軌道は smooth24 親の境界：完全有限切片からの実曲線

2026-11-08、PRIMARY E。全350式の全次数十を確認した producer PASS。
独立 Fraction 再生は未完。新しい成分は計数しない。
既知 MPS24 のどれか、又は未登録24成分かという親の識別は別に残す。

## 完全切片の結果

`K4A_COMPLETE_TORUS_SLICE.md` の事前C*次数境界とactual action局所積を使う。
同じ D101・G48・既N23・既101minorだけから四つの固定係数 X2,...,X5 を解いた。
全input/output weight対応、true102切片、全101選択式の次数十までの零、
全350式の全次数0,...,10を確認した結果、**全て厳密に零**である。
0.225秒、RSS50 MiB。追加point、rank、minor、normal補正次数を選んでいない。

従って原点のformal切片はsmooth A1、actual48作用方向とのétale productからraw局所schemeはsmooth49。
Euler24とscale1を除いた P0 は **smooth unique projective24** 点である。
以前のD101/G48によるtangent24/orbit23 limitationを、完全有限モデルで実24親へ進めた。

## 整数だけで書ける literal actual family

保存graphのparameter t を `t=174e` と実再parameter化し、

```
P(e)=P0+e*A1+e²*A2+e³*A3+e⁴*A4+e⁵*A5
```

と書く。座標 `(c,a,b,u,v)`、P0 は同じ18倍原始整数K4a tensor：

```
P0ca=4ab, P0cb=27a²-2b², P0cu=-6av+4bu, P0cv=-8bv+24au,
P0ab=18av, P0au=18u², P0av=36uv, P0bu=18uv,
P0bv=18u²+36v², P0uv=0.
```

以下に各 A_i の全非零 entry を示す。省略entryは零、反対称性を補う。

```
A1ca=-72bc, A1cb=-420ac, A1cu=37a²-6b²,
A1cv=12ab-36cu, A1ab=-261a²+54b², A1au=36bu, A1bu=174au;

A2cb=216c², A2cu=-126ca, A2ab=-324ca, A2au=27a²,
A2av=-216cu+72ab, A2bu=108cu+18ab, A2bv=285a²;

A3cu=216c², A3ab=432c², A3au=-180ca, A3av=216cb,
A3bu=-36cb, A3bv=1044ca, A3uv=-93a²+18b²;

A4au=324c², A4bv=-1404c², A4uv=612ca;

A5uv=-1080c².
```

これは形式jetでなく、全eについて元350式を満たす斉次二次のactual polynomial tensorである。
旧literal pencil P0+tN が自己Jacobi非零だった事実と矛盾しない。
完全equivariant graphの必須二次から五次係数を同時に残した実族である。

## 全非零 parameter は一つの24軌道

normalized作用 `s^-1 diag(s²,s,s,1,1)_*` のもとで

```
P(e) -> P(e/s).
```

各A_jのcharacterは-jなので全非零eは同じactual projective orbitに入る。
e=0で全actual作用48列にcurve derivativeを加えるとrank49。
これは既Nの非group性（functional値8/87、再parameter後16）だけから分かり、新minorを計算しない。
rankの非消失開性より、一般eの群×curve像はraw49次元。
e非零ではcurve方向自身が上のtorus作用方向なので、同じraw49は群作用像だけで実現する。
従って非零orbitはprojective24次元。その閉包はP0を含むunique24成分そのものである。

wholeK4a sourceは既受理全normalizer theoremによりP0の23軌道へ入るから、
**全K4a sourceがこのactual24成分の境界に入る**。
親を既知catalogueへ認識すること、又は全prior126と区別して新成分として採録することは、
別作者の完全切片再生と別の幾何監査を経て統合担当が行う。
