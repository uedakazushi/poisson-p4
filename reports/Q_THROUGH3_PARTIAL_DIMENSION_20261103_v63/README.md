# 新しい二次57式＋真の三次17式：部分基底と上界22

[証明・停止理由・範囲](../../research/notes/Q_THROUGH3_PARTIAL_DIMENSION_20261103_v63.md)。最終状態は `summary.json`。

一回のcached M2計算から、次数3の329多項式と全cofactorを保存した。全329恒等式 A*C=B を別プロセスでも再生した。13列の先頭単項式はsupportが互いに素であり、高さ13の明示的な完全交差を与える。従って特殊ファイバー次元は22以下。元の74式はQQ上の斉次式で32003-integralなので、明示的なV-modelのProjのpropernessを使い、QQおよびC上の全35変数完備局所モデルの次元も22以下となる。

目標12は未達。次数4の継続は3GiB制約下でメモリ不足により終了した。再度の重いGBや四次DAGは実行していない。

## 保存済みの13列をGBなしで検証

未使用のworkファイル名で実行する。

```sh
source /workspace/cloud-setup/activate.sh
python - <<'PYREPLAY'
import gzip
from pathlib import Path
p=Path('reports/Q_THROUGH3_PARTIAL_DIMENSION_20261103_v63')
text=gzip.decompress((p/'input.m2.gz').read_bytes())+gzip.decompress((p/'thirteen-cofactors.m2.gz').read_bytes())
Path('work/q-thirteen-replay.m2').write_bytes(text+b'print "PASS13"; exit 0;\n')
PYREPLAY
M2 --script work/q-thirteen-replay.m2
```

全329列の場合は `thirteen-cofactors.m2.gz` を `stage-3.witness.m2.gz` に置き換える。いずれも74式のGBを計算せず、保存された多項式恒等式を検査する。元の生成元・全変換係数・単項式support・メモリ監視・停止ログを保持している。
