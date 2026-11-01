# 追加の M0 無限遠点の孤立被約性

2026年11月1日。固定一次部分での 134 式、135 式の二次斉次化を考える。
`reports/S61_MINIMAL_AXIS_LOCAL_POINT_20261101.json` に指定した有理点は、
全 135 最高式を満たす。M0 正常形では
G=−19x²z/2−4x²w/3+5xz²+16xzw/3−10z²w/3 である。

元の S61 式から最高 Jacobian と斉次化変数 t の列を直接再構成すると、
最高 Jacobian の階数は44、射影 chart の全45変数の Jacobian 階数は45となる。
非零の45次小行列は最初の134行だけから選ばれており、証明書に行番号と厳密な有理行列式を保存した。

従って両方の未飽和斉次化の局所環で最大 ideal の余接空間は零である。
Noether 局所環に Nakayama の補題を適用して最大 ideal は零、局所環は C となる。
特に t=0 であり、この点は t による飽和後の閉包から除かれる。
同じ小行列が非零の開部分にこの結論は延びる。
この一点の証明書だけで対応する族の β への優越性や全無限遠の尽くしを主張しない。

再現：

```sh
.venv/bin/python code/s61_minimal_axis_local_point_20261101.py --output /tmp/minimal-axis-local-point.json
```
