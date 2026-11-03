# Cloud 研究再開時の接続確認 — 2026-11-03

記録時刻（このセッションの確認結果）：UTC `2026-11-03T02:02:41+00:00` / Asia/Tokyo `2026-11-03T11:02:41+09:00`。
開始時・照合時の GitHub `refs/heads/main` はともに
`4c3f623fb698f9adc6077015b25e2325b11072af`。

## 実行経路と GitHub

Cloud runtime skill と networking reference、`/etc/codex/network-policy.json` を確認。
`environment_status` は `observations_current: true`、`observed_phase: running`、
HTTP policy `restricted/enforced`。7 文献ホストは設定済み、VPN は未設定。
認証情報の値は表示・変更していない。公開済み環境のセットアップは再実行していない。

最初の通常 `exec_command` による fetch/curl は、`proxy:8080` への接続失敗で終了した。
`getent hosts proxy` は名前解決に成功した。これはホスト個別の拒否を示すものではない。
同じ inherited proxy / CA trust を維持したまま、次の正規 tool grant を付けると接続できた。

```json
{"sandbox_permissions":"with_additional_permissions", "additional_permissions":{"network":{"enabled":true}}}
```

その経路で実行した確認：

```sh
git ls-remote origin refs/heads/main
# 4c3f623fb698f9adc6077015b25e2325b11072af  refs/heads/main
```

fetch/push にも同じ tool network grant を使う。proxy の解除、直接接続、TLS 検証無効化、
network policy の編集はしていない。この報告は push 成功の証拠ではない。

## 文献サイト

各 URL に対し原則次の GET を実行し、本文は保存せず状態と最終 URL を検査した。
`curl` の終了値0は HTTP 状態と別に記録する。

```sh
curl --max-time 20 --silent --show-error --location --output /dev/null \
  --write-out '%{http_code} %{url_effective} %{content_type}' URL
```

| 設定済みホストの入口 | HTTP | curl 終了値 | 結果 |
|---|---:|---:|---|
| `https://arxiv.org/` | 200 | 0 | 接続成功 |
| `https://export.arxiv.org/` | 200 | 0 | 接続成功 |
| `https://doi.org/` | 200 | 0 | `https://www.doi.org/` へ転送 |
| `https://pmc.ncbi.nlm.nih.gov/` | 200 | 0 | 接続成功 |
| `https://stacks.math.columbia.edu/` | 200 | 0 | 接続成功 |
| `https://www.cambridge.org/` | 200 | 0 | 接続成功 |
| `https://www.dpmms.cam.ac.uk/` | 200 | 0 | 接続成功 |

既存研究ノート・文献表の実参照先も確認した（Cambridge/DOI は上限25秒）。

| 実参照先 | HTTP / 終了値 | 確認範囲 |
|---|---|---|
| `https://arxiv.org/html/2010.08692v2` | 200 / 0 | HTML、MPS 論文本文 |
| `https://stacks.math.columbia.edu/tag/02OF` | 200 / 0 | HTML、proper + finite fibers の参照先 |
| `https://pmc.ncbi.nlm.nih.gov/articles/PMC5818673/` | 200 / 0 | HTML、Pym 総説 |
| `https://www.dpmms.cam.ac.uk/~taf1000/papers/invenq.pdf` | 200 / 0 | `application/pdf`、Fisher 著者公開PDF |
| `https://www.cambridge.org/core/journals/proceedings-of-the-edinburgh-mathematical-society/article/nef-vector-bundles-on-a-quadric-threefold-with-first-chern-class-two/A30A0DFEB6D7BA030E24F4B8D72DD982` | 200 / 0 | HTML、既存ノートの nef bundles 記事 |
| `https://doi.org/10.1112/blms.12339`（転送なし） | 302 / 0 | DOI 解決は成功 |
| 同 DOI（転送あり） | 000 / 56 | `https://onlinelibrary.wiley.com/doi/abs/10.1112/blms.12339` で proxy `CONNECT ... 403` |
| `https://export.arxiv.org/api/query?id_list=2010.08692` | 000 / 28 | 20秒、0 bytes の timeout。入口200だけからAPI可用性を推論しない |

## 限定と再開操作

Wiley の転送先ホストは許可一覧にない。DOI サービス自体の接続成功は、転送先の全文取得成功を
意味しない。必要なら Cloud 環境の正式な設定画面で `onlinelibrary.wiley.com`、既存の直接引用URLの
`londmathsoc.onlinelibrary.wiley.com` を許可する設定変更を利用者レビューにかける。
このセッションでは設定変更や拒否の迂回を行っていない。手元の証拠や許可済み一次資料を用いる研究は継続できる。
arXiv export API は一回の timeout で原因を確定できない。必要時に同じ許可経路で上限時間を明示して再試行する。

これらは接続・コンテンツ種別の確認であり、今回新たに文献の定理・仮定を精査したという主張ではない。
