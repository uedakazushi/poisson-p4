# 六quadの全22 syzygy加群の独立監査

2026-11-03。著者のSymPy producerをimport・再実行せず、元の六quadの保存係数を
別のFraction sparse辞書で読み、元8座標から七light座標への変換を再構成する。
保存したL、ell、Tz、Twについて、二つのZW制限式、全六quadのbasis変更、
liftした全syzygy恒等式をQQ上で照合する。168列行列、新点、有限体は用いない。

独立性に使う二つの単元恒等式も全係数で確認する。

- `ell_x1=-16*C*J/243`。
- `Delta=ell|ZW-(2J/1215)L`はx1係数零で、
  `3*Delta_x2-(7/2)*Delta_x4=-80*a*d²/81`。

既受理のZW pair定理から六quadは定数上独立である。15 Koszulと
`W*(ell Z)-Z*(ell W)`のlight倍との重なりはell倍だけなので、x1を除いた
六lightを選べば旧21方向は独立。新syzygyが旧21方向に入るならZW制限の
q1/q3成分からell|ZWとLの比例が必要となり、第二単元恒等式に反する。
従って全FULLUで22方向は点別独立である。

ここで主張するkernel等号は**加群として**のものである。22方向が作る
subbundleをR168から割るとlocally free rank146となる。既知generic rank146に
より、この商からR210への写像はgenericに単射。integralな底上のlocally free
加群はtorsion-freeだから加群kernelは零であり、22本が全kernelのbasisとなる。
像加群もlocally free rank146である。ただし特殊fiberでこの埋込みが単射とは
限らず、全点ambient rank146、cokernel flatness、全quartic全射は未証明である。

再現は以下の通り。既存証明書を上書きしない出力先を指定する。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/s61_nonsquare_universal_quad22_review_20261103_v63.py \
  --output work/universal-quad22-independent-new
```
