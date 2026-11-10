# GitHub mainからの手動再現

2026-11-10 / v12。GitHub mainのみ正本。Actions不使用、Library・旧ZIP・未追跡workは入力にしない。生成PDFはGit追跡しない。

## 環境

Python 3、SymPy 1.14.0、LuaLaTeX（ltjsarticle、luatexja-preset、HaranoAji、amsmath、mathtools、fix-cm、geometry、hyperref、xurl）、Poppler（pdfinfo、pdftotext）。既存全分類の再現にはroot requirements.txtおよび各producerが指定するCAS等の依存も必要。新しい接錐producerはSymPyだけで完結し、旧producer、外部データ、ネットワークを使わない。

## v12：非被約接錐の全検算とPDF

```bash
bash scripts/verify_github_source.sh tangent-cone
```

check_isotropic_quadratic_cone.pyは有理数体上で次を構成・照合する：全350×150 Jacobi線形化のrank95、Eulerゲージ・尺度・30接方向のrank55、全255左核×465接方向対の二次写像のrank7、正規化された7式と疎な左核証明書、二次イデアルのradical分解、平方零元とcolon、対数型源の28次minor・augmented rank29、混合ac補正の全350係数。
保存先 research/subvariety_incidence/certificates/isotropic_quadratic_cone.json と完全照合する。任意の出力先を指定した直接実行も可能：

```bash
python3 code/subvariety_incidence/check_isotropic_quadratic_cone.py \
  --output work/isotropic_quadratic_cone.json \
  --check research/subvariety_incidence/certificates/isotropic_quadratic_cone.json
```

続いて新稿 papers/isotropic_log_tangent_cone_ja.tex をLuaLaTeX二回で組版し、参照・Citation未解決、Overfull、空のPDF抽出がないことを確認する。出力はwork/isotropic-tangent-cone-check/以下のみ。通常の組版はrootの.latexmkrcにより latexmk papers/isotropic_log_tangent_cone_ja.tex でも可能。

## 既存incidenceと分類

```bash
bash scripts/verify_github_source.sh incidence
bash scripts/verify_github_source.sh full
```

incidenceは旧11 producer/11証明書と同伴incidence原稿の再生成に加え、新しい接錐検査も実行する。fullは既存依存hash、147カタログ、global structure、四つの分類原稿の組版とincidenceを行う。これを全147の原証明の自動認証と呼ばない。既存有限証明の射程は research/publication/REFEREE_REPORT.md、research/reproduction/README.md と各元証明を参照。

必要なら python scripts/make_structural_bundle.py --output work/current-github-source-bundle で旧分類用の現行一式を生成できる。新接錐結果の最小再現単位は新TeX、単独producer、JSON、verify_isotropic_tangent_cone.shおよび現在の方針・ノートであり、旧bundleへの依存はない。将来bundle generatorを使う場合には新稿も含まれることを別途確認する。

## 今回の実行範囲

新producerを実行して全有理計算PASS、再実行で保存JSONに一致。別の直接SymPy Jacobi実装で二つのテンソルの全350係数も照合した。新6頁原稿を二回組版し、未解決参照・Overfullなし、全6頁を画像で確認した。旧11 producerの一括再実行、147分類全体の再監査、旧incidence全稿の今回の組版、fresh cloneからの完全再現は未実行。

GitHub保存は期待したmain先端に対する非force更新後、remote SHAと変更ファイルを読み戻して確認する。GitHub書込み、数学的検証、PDF生成、Library退避は別々の状態。通常Library保存はせず、実push失敗時のみ緊急退避する。既存Libraryファイルは無断削除しない。原稿に公表価値のある既存証明を失わせず、TeXにないPDF-onlyの数学入力は削除しない。
## 完全等方Log点の明示23次元交叉族（2026-11-10 / v14）

`code/subvariety_incidence/check_isotropic_intersection_dim23.py` はSymPy 1.14.0で独立に4パラメータPoisson構造を構築し、全Jacobi、rank2 Log1112のepsilon族に対する全10二形式係数恒等式、rank4非定数g長鎖源に対する射影線形変換と全六括弧、24次minor 10871635968、Euler gaugeとの独立性を検証する。保存JSONは `research/subvariety_incidence/certificates/isotropic_intersection_dim23.json`。手動入口 `bash scripts/verify_github_source.sh tangent-cone` で従来の接錐・軌道の証明書と同時に再計算・照合する。数値下界23は一般交叉次元不等式からも出るが、両成分への具体的な源レベルの構成を新しく与えた。交叉の正確な次元23/24、全接錐初期イデアル、ambient局所環の被約性は未解決。

