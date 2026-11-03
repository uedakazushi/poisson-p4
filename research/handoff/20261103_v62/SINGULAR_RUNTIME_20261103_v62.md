# 既存 Cloud ライブラリを使う薄い Singular 実行系

2026-11-03。Macaulay2 の QQ module Gröbner 計算で実際に600秒の時間切れが生じたため、別の既存算術 engine を試せるようにした。公開済み環境には既に Singular 4.4.1 の共有ライブラリと `singular-data` があり、console executable と高速演算 module だけが無かった。同じ Debian trixie version の `singular-ui`（13,648 bytes）と `singular-modules`（1,201,912 bytes）だけを追加で取得し、既存 apt metadata の SHA256 と照合して、研究用の ignored `work/` 内へ展開した。既存の環境セットアップ、M2、共有ライブラリ、proxy/CA 設定、managed configuration は変更していない。

実行 wrapper は現在

    /workspace/poisson-p4/work/nonsquare-v62/singular-runtime/run-singular.sh

である。既存 `activate.sh` を読み、`SINGULARPATH` と root/data/procs directory を設定する。実行例は

    work/nonsquare-v62/singular-runtime/run-singular.sh --cpus=0 -q input.sing

で、通常の exact `std`、`lift` に使える。`--cpus=2` の modular routines は worker 間の localhost socket を使う。この desktop sandbox の既定 network namespace では bind が禁止されていたため、その実行だけ `exec_command` の `with_additional_permissions` と `network.enabled=true` を使った。HTTP policy を変更したり proxy を迂回したりしていない。

再取得する公開 package は次の二つである。URL の host は既存の package-manager network policy に含まれる `deb.debian.org`。

- `https://deb.debian.org/debian/pool/main/s/singular/singular-ui_4.4.1+ds-2_amd64.deb`、SHA256 `4e28bb5d607d2d4f2174dbb19ec8e82a33687139a63faccd2dcdcc3776808dff`
- `https://deb.debian.org/debian/pool/main/s/singular/singular-modules_4.4.1+ds-2_amd64.deb`、SHA256 `c156abdfc3370f84232c57b51ba6d78dff68347d8cf52fc25a3106a094b2426c`

`dpkg-deb -x` で個別 prefix へ展開するだけで root install は不要である。wrapper、実際の binary/共有ライブラリ/算術 library source、軽量 fixture の hash を `singular_runtime_manifest_20261103_v62.json` に保存した。

runtime の軽量 fixture は成功した。今回の数学への適用では、raw inhomogeneous `modStd(M,1)` の unit-basis 出力を、それだけで新定理には採用していない。公開 library が exact と呼ぶ flag と、この研究の有限証明書に必要な両方向の module equality の監査を区別した。実際の二つの thin matrix の raw 出力は候補であり、直接 lift/CRT right-inverse は制限時間内に得られていない。

代わりに homogeneous certificate の一般補題を独立に確認した。各元 column を固定次数で斉次化し、同じ target shifts=0 と module order を保つ。候補 G の QQ Gröbner 性と元 module の G への包含を exact に確認し、一つの良い素数で先頭 module の保存と Gbar=元 module bar の両方向を exact に確認すれば、次数ごとの rank の半連続性が逆包含も与える。これを検査する corrected generator は `homogeneous-module-input.py`。最初の prototype は Singular が module 同士の `==` を定義しないことに当たり、matrix cast による比較へ直した。エラーのある prototype log は証明に使わない。corrected V2 は各600秒まで継続した。J=0 は明示的な6,000,000 KiB virtual-memory cap に達して exit14（`no more memory`）、J≠0 は600秒で exit124となり、どちらも候補GBを返さなかった。これは計算資源内で未完了という記録であり、数学的な不可能性や一般的な engine の限界を意味しない。終了後に全worker descendantの消滅を確認した。`--cpus=2` でも modular library は多数の forked task を作りうるので、次の実行は実際のprocess treeとcgroup memoryを監視する。

従って新しい Singular runtime を使えることと、未完了の大きい QQ module 計算が解決したことは別である。J=0 の product 欠損の有限性は、別途、FLINT の四つの exact QQ minor と gcd によって証明した。
